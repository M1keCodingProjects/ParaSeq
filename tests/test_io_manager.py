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

# uint------------------------------------------------------------------------------------
def test_uintEmpty():
    with pytest.raises(ValueError):
        uint("")

def test_uintAlpha():
    with pytest.raises(ValueError):
        uint('d')

def test_uintAlphaNum():
    with pytest.raises(ValueError):
        uint("d3")

def test_uintDecimal():
    with pytest.raises(ValueError):
        uint("2.2")

def test_uintNeg():
    with pytest.raises(ValueError):
        uint("-1")

def test_uintExp():
    with pytest.raises(ValueError):
        uint("2e04")

def test_uint():
    assert uint("12345") == 12345

# parseRawSeq-----------------------------------------------------------------------------
def test_parseRawSeqEmpty(): # This should never happen normally
    seqObj = parseRawSeq("", SeqName.Query)
    assert seqObj.seq  == ""

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

def test_parseRawSeqInvalid():
    with pytest.raises(ValueError):
        parseRawSeq("QWERTY", SeqName.Query)

def test_parseRawSeqOneInvalid():
    with pytest.raises(ValueError):
        parseRawSeq("GA.TC", SeqName.Query)

def test_parseRawSeqQueryLower():
    seqObj = parseRawSeq("acGt", SeqName.Query)
    assert seqObj.seq  == "ACGT"

# loadFasta-------------------------------------------------------------------------------
def test_loadFastaEmpty(): # This should never happen normally
    with pytest.raises(FileNotFoundError):
        loadFasta("")

def test_loadFastaNonexistent():
    with pytest.raises(FileNotFoundError):
        loadFasta(".fa")

# Always run pytest from the root ParaSeq folder.
TEST_DATA_PATH = ".\\data\\{}.fasta"
def test_loadFasta():
    seqObj = loadFasta(TEST_DATA_PATH.format("good"))[0]
    assert seqObj.seq  == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seqObj.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seqObj.id   == "seq1|accession=XM_123456.7|organism=Homo"

def test_loadFastaEmptyFile():
    with pytest.raises(IndexError):
        loadFasta(TEST_DATA_PATH.format("empty"))[0]

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

def test_parseFastaSeqOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqsErr):
        parseFastaSeq(loadFasta(path), 45, path)

# parseSeq--------------------------------------------------------------------------------
# Coverage is lower here as this function just calls other functions
def test_parseSeq():
    seqObj = parseSeq(TEST_DATA_PATH.format("good"), 0, SeqName.Target)
    assert seqObj.seq  == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seqObj.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seqObj.id   == "seq1|accession=XM_123456.7|organism=Homo"

def test_parseSeqBadPath():
    with pytest.raises(ValueError):
        parseSeq("good.fsta", 0, SeqName.Target)
        # ^^^ Bad file ext leads to interpreting this as DNA

def test_parseSeqOOB():
    with pytest.raises(MissingSeqsErr):
        parseSeq(TEST_DATA_PATH.format("good"), 45, SeqName.Target)

# parseSeqsFromFile----------------------------------------------------------------------
def test_parseSeqsFromFile():
    seq1, seq3 = parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 2)
    assert seq1.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seq3.name == "seq3|accession=NR_000000.1|organism=Escherichia"

def test_parseSeqsFromFileOOB():
    with pytest.raises(MissingSeqsErr):
        parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 3)

def test_parseSeqsFromFileDefault():
    seq1, seq2 = parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 0)
    assert seq1.name == "seq1|accession=XM_123456.7|organism=Homo"
    assert seq2.name == "seq2|accession=NM_654321.4|organism=Mus"

def test_parseSeqsFromFileIdentical():
    with pytest.raises(IdenticalSeqsErr):
        parseSeqsFromFile(TEST_DATA_PATH.format("good"), 1, 1)

# setupArgParser--------------------------------------------------------------------------
def test_setupArgParser():
    args = setupArgParser().parse_args(('0', '1', "-m", '2', "-mm", '3', "-g", '4', '-tp', '5'))
    assert args.target_seq       == '0'
    assert args.query_seq        == '1'
    assert args.match_score      == 2
    assert args.mismatch_penalty == 3
    assert args.gap_penalty      == 4
    assert args.target_pos       == 5

def test_setupArgParserNoArgs():
    with pytest.raises(SystemExit):
        setupArgParser().parse_args(())

def test_setupArgParserMissingArgs():
    with pytest.raises(SystemExit):
        setupArgParser().parse_args(("-m", "2", "-mm", '3', "-g", '4'))
    
    with pytest.raises(SystemExit):
        setupArgParser().parse_args(('0', "-mm", '3', "-g", '4'))
    
    with pytest.raises(SystemExit):
        setupArgParser().parse_args(('0', "-m", "2", "-g", '4'))
    
    with pytest.raises(SystemExit):
        setupArgParser().parse_args(('0', "-m", "2", "-mm", '3'))

def test_setupArgParserInvalidArgs():
    with pytest.raises(SystemExit):
        setupArgParser().parse_args(('0', '1', "-m", "foo", "-mm", '3', "-g", '4'))
    
    with pytest.raises(SystemExit):
        setupArgParser().parse_args(('0', '1', "-m", '2', "-mm", "foo", "-g", '4'))
    
    with pytest.raises(SystemExit):
        setupArgParser().parse_args(('0', '1', "-m", '2', "-mm", '3', "-g", "foo"))

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
    with pytest.raises(IdenticalSeqsErr):
        parseInputArgs(args)

def test_parseInputArgsSameFasta():
    path = TEST_DATA_PATH.format("good")
    args = setupArgParser().parse_args((path, path, "-m", '2', "-mm", '3', "-g", '4'))
    parseInputArgs(args)

def test_parseInputArgsSingleDNA():
    args = setupArgParser().parse_args(("ACG", "-m", '2', "-mm", '3', "-g", '4'))
    with pytest.raises(MissingSeqsErr):
        parseInputArgs(args)

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
