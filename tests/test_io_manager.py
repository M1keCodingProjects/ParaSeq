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

# validateDNA-----------------------------------------------------------------------------
def test_validateDNA():
    assert validateDNA("ACGT") == "ACGT"

def test_validateDNALower():
    assert validateDNA("acGt") == "ACGT"

@pytest.mark.parametrize("seq", ["", "QWERTY", "GA.TC"])
def test_validateDNAInvalid(seq):
    with pytest.raises(InvalidSeqErr) as errInfo: validateDNA(seq)
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + ": ."

# parseFastaSeq---------------------------------------------------------------------------
def test_parseFastaSeqEmpty(): # This should never happen normally
    with pytest.raises(FileNotFoundError) as errInfo: parseFastaSeq("", 0)
    assert str(errInfo.value) == "[Errno 2] No such file or directory: ''"

def test_parseFastaSeqNonexistent():
    with pytest.raises(FileNotFoundError) as errInfo: parseFastaSeq(".fa", 0)
    assert str(errInfo.value) == "[Errno 2] No such file or directory: '.fa'"

# Always run pytest from the root ParaSeq folder.
TEST_DATA_PATH = ".\\data\\{}.fasta"
def test_parseFastaSeqEmptyFile():
    path = TEST_DATA_PATH.format("empty")
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq(path, 0)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": FASTA file \"{path}\" doesn't c\
ontain a sequence at position 0."

def test_parseFastaSeq():
    seq1, seq2 = parseFastaSeq(TEST_DATA_PATH.format("good"), 0, 1)
    assert seq1 == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seq2 == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"

def test_parseFastaSeqSwitchedIndexes():
    seq1, seq2 = parseFastaSeq(TEST_DATA_PATH.format("good"), 2, 0)
    assert seq1 == "GCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCATCGATCGTAGCTAGCTAGCTAA"
    assert seq2 == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"

def test_parseFastaSeqMultiline():
    seq1, seq2 = parseFastaSeq(TEST_DATA_PATH.format("multiline"), 0, 1)
    assert seq1 == "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    assert seq2 == "GGCCTTAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAG"

# The function doesn't distinguish between no sequence and bad sequence:
def test_parseFastaSeqHeaderOnly():
    path = TEST_DATA_PATH.format("header")
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq(path, 0)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": FASTA \
file \"{path}\" doesn't contain a sequence at position 0."

def test_parseFastaSeqInvalid():
    path = TEST_DATA_PATH.format("bad")
    with pytest.raises(InvalidSeqErr) as errInfo: parseFastaSeq(path, 0)
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + f": sequence at position 0 in FASTA\
 file \"{path}\" is not valid DNA."

def test_parseFastaSeqOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq(path, 45)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": FASTA \
file \"{path}\" doesn't contain a sequence at position 45."

# parseSeq--------------------------------------------------------------------------------
# Coverage is lower here as this function just calls other functions
def test_parseSeqRaw():
    assert parseSeq("ACGTGT", 300, SeqName.Target) == "ACGTGT"
    # ^^^ pos is ignored

def test_parseSeq():
    assert parseSeq(TEST_DATA_PATH.format("good"), 0, SeqName.Target) == "ATGCGTACGTAGCTA\
GCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"

def test_parseSeqBadPath():
    with pytest.raises(InvalidSeqErr) as errInfo: parseSeq("good.fsta", 0, SeqName.Target)
    # ^^^ Bad file ext leads to interpreting this as DNA
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + f": Your \"{
        SeqName.Target.value}\" sequence is not valid DNA."

def test_parseSeqOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo:
        parseSeq(path, 45, SeqName.Target)
    
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": FASTA \
file \"{path}\" doesn't contain a sequence at position 45."

# parseSeqsFromFile----------------------------------------------------------------------
def test_parseSeqsFromFile():
    seq1, seq2 = parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 2)
    assert seq1 == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seq2 == "GCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCATCGATCGTAGCTAGCTAGCTAA"

def test_parseSeqsFromFileOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo:
        parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 3)
    
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": FASTA \
file \"{path}\" doesn't contain a sequence at position 3."

def test_parseSeqsFromFileDefault():
    seq1, seq2 = parseSeqsFromFile(TEST_DATA_PATH.format("good"), 0, 0)
    assert seq1 == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seq2 == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"

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
    assert target == "ACG"
    assert query == "CGT"
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
    assert target == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert query == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"
    assert match == 2
    assert mismatch == 3
    assert gap == 4

# collectAndParseInputArgs----------------------------------------------------------------
def test_collectAndParseInputArgs():
    target, query, match, mismatch, gap = collectAndParseInputArgs((TEST_DATA_PATH.format("good"), "-m", '2', "-mm", '3', "-g", '4'))
    assert target == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert query == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"
    assert match == 2
    assert mismatch == 3
    assert gap == 4