from src.para_seq import UINT_ERR
from src.para_seq.io_manager import *
import pytest

# Tests for CustomErr extensions have been skipped, add them if new bahaviour is added.

# isValidFastaFilePath--------------------------------------------------------------------
def test_isValidFastaFilePathEmpty():
    assert not isValidFastaFilePath("")

def test_isValidFastaFilePathOther():
    assert not isValidFastaFilePath("file.txt")

def test_isValidFastaFilePath():
    assert isValidFastaFilePath("file.fasta")

def test_isValidFastaFilePathShort():
    assert isValidFastaFilePath("file.fa")

def test_isValidFastaFilePathTechnically():
    assert isValidFastaFilePath(".fa")

# isValidDNA------------------------------------------------------------------------------
def test_isValidDNAEmpty():
    assert not isValidDNA("")

def test_isValidDNA():
    assert isValidDNA("ACGTNTGCGCTACACTGNNNNGGGTTTCCCAAA")

def test_isValidDNAMixedCase():
    assert isValidDNA("ACGaTNtaTGttCGccCgTcAggCcAtCaTcGnnNNtNnNaGnGaGaTTTCggCCAAA")

def test_isValidDNAInvalid():
    assert not isValidDNA("ACTFG")

# uint------------------------------------------------------------------------------------
@pytest.mark.parametrize("value", ["", "d", "d3", "2.2", "-1", "2e04"])
def test_uintInvalid(value):
    with pytest.raises(ValueError) as errInfo: uint(value)
    assert str(errInfo.value) == UINT_ERR.format(value)

def test_uint():
    assert uint("12345") == 12345

# parseRawSeq-----------------------------------------------------------------------------
def test_parseRawSeqQuery():
    seqObj = parseRawSeq("ACGT", SeqName.Query)
    assert seqObj.seq  == "ACGT"
    assert seqObj.name == SeqName.Query.value
    assert seqObj.id   == SeqName.Query.id

def test_parseRawSeqTarget():
    seqObj = parseRawSeq("ACGT", SeqName.Target)
    assert seqObj.seq  == "ACGT"
    assert seqObj.name == SeqName.Target.value
    assert seqObj.id   == SeqName.Target.id

def test_parseRawSeqQueryLower():
    seqObj = parseRawSeq("acGt", SeqName.Query)
    assert seqObj.seq  == "ACGT"

@pytest.mark.parametrize("seq", ["", "QWERTY", "GA.TC"])
def test_parseRawSeqInvalid(seq):
    with pytest.raises(InvalidSeqErr) as errInfo: parseRawSeq(seq, SeqName.Query)
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + ": " + "invalid raw DNA string."

# loadFasta-------------------------------------------------------------------------------
def test_loadFastaEmpty(): # This should never happen normally
    with pytest.raises(FileNotFoundError) as errInfo: loadFasta("")
    assert str(errInfo.value) == "[Errno 2] No such file or directory: ''"

def test_loadFastaNonexistent():
    with pytest.raises(FileNotFoundError) as errInfo: loadFasta(".fa")
    assert str(errInfo.value) == "[Errno 2] No such file or directory: '.fa'"

# Always run pytest from the root ParaSeq folder.
TEST_DATA_PATH = ".\\data\\{}.fasta"
def test_loadFasta():
    seqObj = loadFasta(TEST_DATA_PATH.format("good"))[0]
    assert seqObj.seq  == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seqObj.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seqObj.id   == "seq1|accession=XM_123456.7|organism=Homo"

# All these tests extract the SeqRecord instance directly, which should never be done in
# the program. This stands to test how the FASTA file is parsed as a whole, with none of
# the necessary processing that must be done on the specific requested sequences.
def test_loadFastaEmptyFile():
    with pytest.raises(IndexError) as errInfo:
        loadFasta(TEST_DATA_PATH.format("empty"))[0]
    
    assert str(errInfo.value) == "list index out of range"

def test_loadFastaHeader():
    seqObj = loadFasta(TEST_DATA_PATH.format("header"))[0]
    assert seqObj.seq  == ""
    assert seqObj.name == "seq1|accession=NR_000000.1|organism=Escherichia"
    assert seqObj.id   == "seq1|accession=NR_000000.1|organism=Escherichia"

def test_loadFastaCorrupted():
    seqObj = loadFasta(TEST_DATA_PATH.format("bad"))[0]
    assert seqObj.seq  == "QWERTYUIOPOIUYTREERTYUI"
    assert seqObj.name == "seq1|access3456.7|org=an=ism=Homo"
    assert seqObj.id   == "seq1|access3456.7|org=an=ism=Homo"

# parseFastaSeq---------------------------------------------------------------------------
def test_parseFastaSeq():
    path = TEST_DATA_PATH.format("good")
    seqObj = parseFastaSeq(loadFasta(path), 0, path)
    assert seqObj.seq  == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seqObj.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seqObj.id   == "seq1|accession=XM_123456.7|organism=Homo"

# The function doesn't distinguish between no sequence and bad sequence:
@pytest.mark.parametrize("path", ["header", "bad"])
def test_parseFastaSeqInvalid(path):
    path = TEST_DATA_PATH.format(path)
    with pytest.raises(InvalidSeqErr) as errInfo: parseFastaSeq(loadFasta(path), 0, path)
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + f": sequence at position 0 in FASTA\
 file \"{path}\" is not valid DNA."

def test_parseFastaSeqOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq(loadFasta(path), 45, path)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": list index out of range, FASTA \
file \"{path}\" doesn't contain a sequence at position 45."

# parseSeq--------------------------------------------------------------------------------
# Coverage is lower here as this function just calls other functions
def test_parseSeq():
    seqObj = parseSeq(TEST_DATA_PATH.format("good"), 0, SeqName.Target)
    assert seqObj.seq  == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seqObj.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seqObj.id   == "seq1|accession=XM_123456.7|organism=Homo"

def test_parseSeqBadPath():
    with pytest.raises(InvalidSeqErr) as errInfo: parseSeq("good.fsta", 0, SeqName.Target)
    # ^^^ Bad file ext leads to interpreting this as DNA
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + f": invalid raw DNA string, your \"\
{SeqName.Target.value}\" sequence was interpreted as DNA, if you intended to provide a FA\
STA file instead make sure your file path ends with a .fa or .fasta extension."

def test_parseSeqOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo:
        parseSeq(path, 45, SeqName.Target)
    
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": list index out of range, FASTA \
file \"{path}\" doesn't contain a sequence at position 45."

# parseSeqsFromFile----------------------------------------------------------------------
def test_parseSeqsFromFile():
    seq1, seq3 = parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 2)
    assert seq1.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seq3.name == "seq3|accession=NR_000000.1|organism=Escherichia"

def test_parseSeqsFromFileOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo:
        parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 3)
    
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": list index out of range, FASTA \
file \"{path}\" doesn't contain a sequence at position 3."

def test_parseSeqsFromFileDefault():
    seq1, seq2 = parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 0)
    assert seq1.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seq2.name == "seq2|accession=NM_654321.4|organism=Mus"

def test_parseSeqsFromFileIdentical():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(IdenticalSeqsErr) as errInfo: parseSeqsFromFile(path, 1, 1)
    assert str(errInfo.value) == f"Alignment of identical sequences is pointless: tried l\
oading sequences coming from the same file ({path}), with identical positions (1)."

# setupArgParser--------------------------------------------------------------------------
def test_setupArgParser():
    args = setupArgParser().parse_args(('0', '1', "-m", '2', "-mm", '3', "-g", '4', '-tp', '5'))
    assert args.target_seq       == '0'
    assert args.query_seq        == '1'
    assert args.match_score      == 2
    assert args.mismatch_penalty == 3
    assert args.gap_penalty      == 4
    assert args.target_pos       == 5

@pytest.mark.parametrize("args", [
    (),
    ("-m", "2", "-mm", '3', "-g", '4'),
    ('0', "-mm", '3', "-g", '4'),
    ('0', "-m", "2", "-g", '4'),
    ('0', "-m", "2", "-mm", '3'),
    ('0', '1', "-m", "foo", "-mm", '3', "-g", '4'),
    ('0', '1', "-m", '2', "-mm", "foo", "-g", '4'),
    ('0', '1', "-m", '2', "-mm", '3', "-g", "foo"),
])
def test_setupArgParserInvalidOrMissingArgs(args):
    with pytest.raises(SystemExit) as errInfo: setupArgParser().parse_args(args)

# parseInputArgs--------------------------------------------------------------------------
def test_parseInputArgs():
    args = setupArgParser().parse_args(("ACG", "CGT", "-m", '2', "-mm", '3', "-g", '4'))
    target, query, match, mismatch, gap = parseInputArgs(args)
    assert isinstance(target, SeqRecord)
    assert isinstance(query, SeqRecord)
    assert target.seq == "ACG"
    assert query.seq == "CGT"
    assert target.name == SeqName.Target.value
    assert query.name == SeqName.Query.value
    assert target.id == SeqName.Target.id
    assert query.id == SeqName.Query.id
    assert match == 2
    assert mismatch == 3
    assert gap == 4

def test_parseInputArgs():
    args = setupArgParser().parse_args(("ACG", "CGT", "-m", '2', "-mm", '3', "-g", '4'))
    target, query, match, mismatch, gap = parseInputArgs(args)
    assert isinstance(target, SeqRecord)
    assert isinstance(query, SeqRecord)
    assert target.seq == "ACG"
    assert query.seq == "CGT"
    assert target.name == SeqName.Target.value
    assert query.name == SeqName.Query.value
    assert target.id == SeqName.Target.id
    assert query.id == SeqName.Query.id
    assert match == 2
    assert mismatch == 3
    assert gap == 4

def test_parseInputArgs():
    args = setupArgParser().parse_args(("ACG", "CGT", "-m", '2', "-mm", '3', "-g", '4'))
    target, query, match, mismatch, gap = parseInputArgs(args)
    assert isinstance(target, SeqRecord)
    assert isinstance(query, SeqRecord)
    assert target.seq == "ACG"
    assert query.seq == "CGT"
    assert target.name == SeqName.Target.value
    assert query.name == SeqName.Query.value
    assert target.id == SeqName.Target.id
    assert query.id == SeqName.Query.id
    assert match == 2
    assert mismatch == 3
    assert gap == 4

def test_parseInputArgs():
    args = setupArgParser().parse_args(("ACG", "CGT", "-m", '2', "-mm", '3', "-g", '4'))
    target, query, match, mismatch, gap = parseInputArgs(args)
    assert isinstance(target, SeqRecord)
    assert isinstance(query, SeqRecord)
    assert target.seq == "ACG"
    assert query.seq == "CGT"
    assert target.name == SeqName.Target.value
    assert query.name == SeqName.Query.value
    assert target.id == SeqName.Target.id
    assert query.id == SeqName.Query.id
    assert match == 2
    assert mismatch == 3
    assert gap == 4

def test_parseInputArgsSameDNA():
    args = setupArgParser().parse_args(("ACG", "ACG", "-m", '2', "-mm", '3', "-g", '4'))
    with pytest.raises(IdenticalSeqsErr) as errInfo: parseInputArgs(args)
    assert str(errInfo.value) == "Alignment of identical sequences is pointless: tried lo\
ading identical DNA sequences, the 2 provided sequences were interpreted as DNA."

def test_parseInputArgsSameFasta():
    path = TEST_DATA_PATH.format("good")
    args = setupArgParser().parse_args((path, path, "-m", '2', "-mm", '3', "-g", '4'))
    parseInputArgs(args)

def test_parseInputArgsSingleDNA():
    args = setupArgParser().parse_args(("ACG", "-m", '2', "-mm", '3', "-g", '4'))
    with pytest.raises(MissingSeqErr) as errInfo: parseInputArgs(args)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + ": a single DNA sequence was provid\
ed."

def test_parseInputArgsSingleFasta():
    args = setupArgParser().parse_args((TEST_DATA_PATH.format("good"), "-m", '2', "-mm", '3', "-g", '4'))
    target, query, match, mismatch, gap = parseInputArgs(args)
    assert isinstance(target, SeqRecord)
    assert isinstance(query, SeqRecord)
    assert target.seq == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert query.seq == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"
    assert target.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert query.name == "seq2|accession=NM_654321.4|organism=Mus"
    assert target.id == "seq1|accession=XM_123456.7|organism=Homo"
    assert query.id == "seq2|accession=NM_654321.4|organism=Mus"
    assert match == 2
    assert mismatch == 3
    assert gap == 4

# collectAndParseInputArgs----------------------------------------------------------------
def test_collectAndParseInputArgs():
    target, query, match, mismatch, gap = collectAndParseInputArgs((TEST_DATA_PATH.format("good"), "-m", '2', "-mm", '3', "-g", '4'))
    assert isinstance(target, SeqRecord)
    assert isinstance(query, SeqRecord)
    assert target.seq == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert query.seq == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"
    assert target.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert query.name == "seq2|accession=NM_654321.4|organism=Mus"
    assert target.id == "seq1|accession=XM_123456.7|organism=Homo"
    assert query.id == "seq2|accession=NM_654321.4|organism=Mus"
    assert match == 2
    assert mismatch == 3
    assert gap == 4