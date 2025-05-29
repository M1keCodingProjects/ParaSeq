from src.para_seq import UINT_ERR
from para_seq.input_manager import *
import pytest

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
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq("", 1)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + ": [Errno 2] No such file or directory: '', the provided FASTA file doesn't exist."

def test_parseFastaSeqNonexistent():
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq(".fa", 1)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + ": [Errno 2] No such file or directory: '.fa', the provided FASTA file doesn't exist."

# Always run pytest from the root ParaSeq folder.
TEST_DATA_PATH = ".\\data\\{}.fasta"
def test_parseFastaSeqEmptyFile():
    path = TEST_DATA_PATH.format("empty")
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq(path, 1)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": provided FASTA file \"{path}\" is empty."

def test_parseFastaSeqMalformedFile():
    path = TEST_DATA_PATH.format("malformed")
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq(path, 1)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": {path} is not plain or gzip compressed fasta formatted file, file is malformed."

def test_parseFastaSeq():
    seq1, seq2 = parseFastaSeq(TEST_DATA_PATH.format("good"), 1, 2)
    assert seq1 == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seq2 == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"

# 0-based positions are not allowed and stopped much earlier, so this should never happen:
def test_parseFastaSeq0():
    seq1, _ = parseFastaSeq(TEST_DATA_PATH.format("good"), 0)
    assert seq1 == "GCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCATCGATCGTAGCTAGCTAGCTAA"
    # 0 becomes -1 to correct from the 1-based system, which ends up loading the last seq

def test_parseFastaSeqSwitchedIndexes():
    seq1, seq2 = parseFastaSeq(TEST_DATA_PATH.format("good"), 2, 1)
    assert seq1 == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"
    assert seq2 == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"

def test_parseFastaSeqMultiline():
    seq1, seq2 = parseFastaSeq(TEST_DATA_PATH.format("multiline"), 1, 2)
    assert seq1 == "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    assert seq2 == "GGCCTTAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAGGCTAAG"

# The function doesn't distinguish between no sequence and bad sequence:
def test_parseFastaSeqHeaderOnly():
    path = TEST_DATA_PATH.format("header")
    with pytest.raises(InvalidSeqErr) as errInfo: parseFastaSeq(path, 1)
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + f": sequence at position 1 in FASTA file \"{path}\" is not valid DNA."

def test_parseFastaSeqInvalid():
    path = TEST_DATA_PATH.format("bad")
    with pytest.raises(InvalidSeqErr) as errInfo: parseFastaSeq(path, 1)
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + f": sequence at position 1 in FASTA file \"{path}\" is not valid DNA."

def test_parseFastaSeqOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo: parseFastaSeq(path, 45)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": FASTA file \"{path}\" doesn't contain a sequence at position 45."

# parseSeq--------------------------------------------------------------------------------
# Coverage is lower here as this function just calls other functions
def test_parseSeqRaw():
    assert parseSeq("ACGTGT", 300, SeqName.Target) == "ACGTGT"
    # ^^^ pos is ignored

def test_parseSeq():
    assert parseSeq(TEST_DATA_PATH.format("good"), 1, SeqName.Target) == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"

def test_parseSeqBadPath():
    with pytest.raises(InvalidSeqErr) as errInfo: parseSeq("good.fsta", 1, SeqName.Target)
    # ^^^ Bad file ext leads to interpreting this as DNA
    assert str(errInfo.value) == INVALID_SEQ_PREFIX + f": Your \"{SeqName.Target.value}\" sequence is not valid DNA."

def test_parseSeqOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo:
        parseSeq(path, 45, SeqName.Target)
    
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": FASTA file \"{path}\" doesn't contain a sequence at position 45."

# parseSeqsFromFile----------------------------------------------------------------------
def test_parseSeqsFromFile():
    seq1, seq2 = parseSeqsFromFile(TEST_DATA_PATH.format("good"), 1, 2)
    assert seq1 == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seq2 == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"

def test_parseSeqsFromFileOOB():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(MissingSeqErr) as errInfo:
        parseSeqsFromFile(TEST_DATA_PATH.format("good"), 1, 45)
    
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + f": FASTA file \"{path}\" doesn't contain a sequence at position 45."

def test_parseSeqsFromFileDefault():
    seq1, seq2 = parseSeqsFromFile(TEST_DATA_PATH.format("good"), 1, 1)
    assert seq1 == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert seq2 == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"

def test_parseSeqsFromFileIdentical():
    path = TEST_DATA_PATH.format("good")
    with pytest.raises(IdenticalSeqsErr) as errInfo: parseSeqsFromFile(path, 2, 2)
    assert str(errInfo.value) == f"{IDENTICAL_SEQS_PREFIX}: tried loading sequences coming from the same file ({path}), with identical positions (2)."

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
    with pytest.raises(SystemExit): setupArgParser().parse_args(args)

# parseInputArgs--------------------------------------------------------------------------
def test_parseInputArgs():
    args = setupArgParser().parse_args((
        "ACG", "CGT", "-m", '2', "-mm", '3', "-g", '4', "-o", "foo", "-ma", '5', "-ls", '6'))
    target, query, match, mismatch, gap, outputPath, maxAlignments, maxSeqLen = parseInputArgs(args)
    assert target == "ACG"
    assert query == "CGT"
    assert match == 2
    assert mismatch == 3
    assert gap == 4
    assert outputPath == "foo"
    assert maxAlignments == 5
    assert maxSeqLen == 6

def test_parseInputArgsSameDNA():
    args = setupArgParser().parse_args(("ACG", "ACG", "-m", '2', "-mm", '3', "-g", '4'))
    with pytest.raises(IdenticalSeqsErr) as errInfo: parseInputArgs(args)
    assert str(errInfo.value) == IDENTICAL_SEQS_PREFIX + ": tried loading identical DNA sequences, the 2 provided sequences were interpreted as DNA."

def test_parseInputArgsSameFasta():
    path = TEST_DATA_PATH.format("good")
    args = setupArgParser().parse_args((path, path, "-m", '2', "-mm", '3', "-g", '4'))
    parseInputArgs(args)

def test_parseInputArgsSingleDNA():
    args = setupArgParser().parse_args(("ACG", "-m", '2', "-mm", '3', "-g", '4'))
    with pytest.raises(MissingSeqErr) as errInfo: parseInputArgs(args)
    assert str(errInfo.value) == MISSING_SEQ_PREFIX + ": a single DNA sequence was provided."

def test_parseInputArgsSingleFasta():
    args = setupArgParser().parse_args((TEST_DATA_PATH.format("good"), "-m", '2', "-mm", '3', "-g", '4'))
    target, query, match, mismatch, gap, outputPath, maxAlignments, maxSeqLen = parseInputArgs(args)
    assert target == "ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG"
    assert query == "CGATCGATGCATGCATGCATGCTAGCTAGCATCGTAGCTAGCTAGCTAACGATCGATCGTAGCGTACG"
    assert match == 2
    assert mismatch == 3
    assert gap == 4
    assert outputPath == "./output/output.txt"
    assert maxAlignments == MAX_DISPLAYED_ALIGNMENTS
    assert maxSeqLen == MAX_DISPLAYED_SEQ_LEN

def test_parseInputArgsInvalidOutPath():
    args = setupArgParser().parse_args(("ACG", "-m", '2', "-mm", '3', "-g", '4', "-o", "./output/"))
    with pytest.raises(InvalidFileErr) as errInfo: parseInputArgs(args)
    assert str(errInfo.value) == INVALID_FILE_PREFIX + ": [Errno 13] Permission denied: './output/', \"./output/\" is not a valid output file path."