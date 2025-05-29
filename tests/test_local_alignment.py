from numpy import any, shape, int64
from para_seq.local_alignment import *
import pytest

# getMatrixShape--------------------------------------------------------------------------
def test_getMatrixShape():
    assert getMatrixShape("ACTG", "TTG") == (4, 5)

def test_getMatrixShapeEmptyTarget():
    assert getMatrixShape("", "TTG") == (4, 1)

def test_getMatrixShapeEmptyQuery():
    assert getMatrixShape("ACTG", "") == (1, 5)

def test_getMatrixShapeSameLen():
    assert getMatrixShape("ACTGACTGACTGACTG", "ACTGACTGACTGACTG") == (17, 17)

# createMatrices--------------------------------------------------------------------------
def test_createMatrices():
    scoreMat, scoreMem, dirsMat, dirsMem = createMatrices((3, 4))
    assert isinstance(scoreMat, ndarray)
    assert isinstance(scoreMem, SharedMemory)
    assert isinstance(dirsMat,  ndarray)
    assert isinstance(dirsMem,  SharedMemory)

    assert scoreMem.name == SCORE_MATRIX_SHMEM_NAME
    assert dirsMem.name  == DIRS_MATRIX_SHMEM_NAME

    # Check if zeros-filled:
    assert not any(scoreMat)
    assert not any(dirsMat)

    assert shape(scoreMat) == (3, 4)
    assert shape(dirsMat)  == (3, 4)

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

# This in theory could never happen as empty seqs are stopped before.
def test_createMatricesZeroDims():
    with pytest.raises(ValueError) as errInfo: createMatrices((0, 0))
    assert str(errInfo.value) == "\'size\' must be a positive number different from zero"

# TODO: Testing the retrieval case needs a process

# freeSharedMem---------------------------------------------------------------------------
def test_freeSharedMem():
    mem = SharedMemory(name = "test", create = True, size = 1)

    assert mem.buf is not None
    freeSharedMem(mem, isFreedCompletely = True)

    # Any access to the buffer after closing is no longer possible:
    assert mem.buf is None

    # Trying to attach an unlinked memory would trigger an error:
    with pytest.raises(FileNotFoundError) as errInfo: SharedMemory(name = mem.name)
    assert str(errInfo) == "<ExceptionInfo FileNotFoundError(2, 'The system cannot find t\
he file specified') tblen=2>"

def test_freeSharedMemOnlyClose():
    mem = SharedMemory(name = "test", create = True, size = 1)

    assert mem.buf is not None
    freeSharedMem(mem)

    assert mem.buf is None

    with pytest.raises(FileNotFoundError) as errInfo: SharedMemory(name = mem.name)
    assert str(errInfo) == "<ExceptionInfo FileNotFoundError(2, 'The system cannot find t\
he file specified') tblen=2>"

# createSharedMatrix----------------------------------------------------------------------
def test_createSharedMatrix():
    mat, mem = createSharedMatrix((3, 4), uint32, "matrix", isNew = True)
    assert isinstance(mat, ndarray)
    assert isinstance(mem, SharedMemory)
    assert shape(mat) == (3, 4)
    assert mem.name == "matrix"
    assert not any(mat)

    freeSharedMem(mem)

# computeAntidiagCoords-------------------------------------------------------------------
def test_computeAntidiagCoords():
    rows, cols = 3, 5
    assert computeAntidiagCoords(0, rows, cols).tolist() == [[0, 0]]
    assert computeAntidiagCoords(1, rows, cols).tolist() == [[1, 0], [0, 1]]
    assert computeAntidiagCoords(2, rows, cols).tolist() == [[2, 0], [1, 1], [0, 2]]
    assert computeAntidiagCoords(3, rows, cols).tolist() == [[3, 0], [2, 1], [1, 2]]
    assert computeAntidiagCoords(4, rows, cols).tolist() == [[4, 0], [3, 1], [2, 2]]
    assert computeAntidiagCoords(5, rows, cols).tolist() == [[4, 1], [3, 2]]
    assert computeAntidiagCoords(6, rows, cols).tolist() == [[4, 2]]

def test_computeAntidiagCoordsTall():
    rows, cols = 5, 3
    assert computeAntidiagCoords(0, rows, cols).tolist() == [[0, 0]]
    assert computeAntidiagCoords(1, rows, cols).tolist() == [[1, 0], [0, 1]]
    assert computeAntidiagCoords(2, rows, cols).tolist() == [[2, 0], [1, 1], [0, 2]]
    assert computeAntidiagCoords(3, rows, cols).tolist() == [[2, 1], [1, 2], [0, 3]]
    assert computeAntidiagCoords(4, rows, cols).tolist() == [[2, 2], [1, 3], [0, 4]]
    assert computeAntidiagCoords(5, rows, cols).tolist() == [[2, 3], [1, 4]]
    assert computeAntidiagCoords(6, rows, cols).tolist() == [[2, 4]]

def test_computeAntidiagCoordsSquare():
    rows, cols = 3, 3
    assert computeAntidiagCoords(0, rows, cols).tolist() == [[0, 0]]
    assert computeAntidiagCoords(1, rows, cols).tolist() == [[1, 0], [0, 1]]
    assert computeAntidiagCoords(2, rows, cols).tolist() == [[2, 0], [1, 1], [0, 2]]
    assert computeAntidiagCoords(3, rows, cols).tolist() == [[2, 1], [1, 2]]
    assert computeAntidiagCoords(4, rows, cols).tolist() == [[2, 2]]

# This in theory could never happen as empty seqs are stopped before.
def test_computeAntidiagCoordsEmpty():
    assert computeAntidiagCoords(0, 0, 0).tolist() == []

# This in theory could never happen as antidiags are computed based on seq len.
def test_computeAntidiagCoordsOOB():
    assert computeAntidiagCoords(3, 1, 1).tolist() == []

# fillMatrices----------------------------------------------------------------------------
def test_fillMatrices():
    _, scoreMem, _, dirsMem = createMatrices((4, 7))
    assert fillMatrices(("ATTTCG", "TTT", 2, 2, 1)) == 6

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

# Should be the same exact score as the previous test, since the sequences are the same
def test_fillMatricesTall():
    _, scoreMem, _, dirsMem = createMatrices((7, 4))
    assert fillMatrices(("TTT", "ATTTCG", 2, 2, 1)) == 6

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

def test_fillMatricesSquare():
    _, scoreMem, _, dirsMem = createMatrices((7, 7))
    assert fillMatrices(("TTAAAT", "ATTTCG", 2, 2, 1)) == 4

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

def test_fillMatricesAllZeros():
    _, scoreMem, _, dirsMem = createMatrices((7, 4))
    assert fillMatrices(("TTT", "ATTTCG", 0, 0, 0)) == 0

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

# This cannot happen, negative scores are invalidated way before this point:
def test_fillMatricesNegativeScores():
    _, scoreMem, _, dirsMem = createMatrices((7, 4))
    assert fillMatrices(("TTT", "ATTTCG", -2, -2, -1)) == 9
    # Apparently it works, it just makes no sense

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

def test_fillMatricesIncompatibleSeqs():
    _, scoreMem, _, dirsMem = createMatrices((7, 4))
    assert fillMatrices(("TTT", "AAAAAA", 2, 2, 1)) == 0

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

# This can't even happen, but let's see what score it gives:
def test_fillMatricesIdenticalSeqs():
    _, scoreMem, _, dirsMem = createMatrices((7, 7))
    assert fillMatrices(("AAAAAA", "AAAAAA", 2, 2, 1)) == 12

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

# reconstructAlignments-------------------------------------------------------------------
def test_reconstructAlignments():
    params = ("ATTTCG", "TTT", 2, 2, 1)
    scoreMat, scoreMem, _, dirsMem = createMatrices((4, 7))
    maxScore = fillMatrices(params)
    
    expectedAlignments = [(2, 1, "TTT", "TTT")]
    assert all(map(
        lambda alignment: alignment in expectedAlignments,
        reconstructAlignments(scoreMat, maxScore, params)))

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

def test_reconstructAlignmentsCommutative():
    # Testing this for multiple alignments gets complicated, we can pretty much infer
    # that seq order doesn't matter
    seq1 = "TTTACATATCGGTGTC"
    seq2 = "ACGCG"
    params = (seq1, seq2, 2, 2, 1)
    scoreMat, scoreMem, _, dirsMem = createMatrices((len(params[1]) + 1, len(params[0]) + 1))
    maxScore = fillMatrices(params)
    assert reconstructAlignments(scoreMat, maxScore, params) == [(8, 1, "ATCG-G", "A-CGCG")]
    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)
    
    params = (seq2, seq1, 2, 2, 1)
    scoreMat, scoreMem, _, dirsMem = createMatrices((len(params[1]) + 1, len(params[0]) + 1))
    maxScore = fillMatrices(params)
    assert reconstructAlignments(scoreMat, maxScore, params) == [(1, 8, "A-CGCG", "ATCG-G")]

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

def test_reconstructAlignmentsMany():
    params = ("ATGCGTACGTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGTTAGCATCGATCGATCGTACGTAGCTAGCTAGCTAACG", "AAAATTTAAAAA", 2, 2, 1)
    scoreMat, scoreMem, _, dirsMem = createMatrices((13, len(params[0]) + 1))
    maxScore = fillMatrices(params)

    expectedAlignments = [(43, 4, "ATCGTTA", "AT--TTA"), (43, 4, "ATCGTTAGCA", "AT--TTA--A")]
    assert all(map(
        lambda alignment: alignment in expectedAlignments,
        reconstructAlignments(scoreMat, maxScore, params)))

    freeSharedMem(scoreMem)
    freeSharedMem(dirsMem)

def test_reconstructAlignmentsZeroScore():
    mat, mem = createSharedMatrix((1, 1), uint32, SCORE_MATRIX_SHMEM_NAME, isNew = True)
    assert reconstructAlignments(mat, 0, ()) == []
    freeSharedMem(mem)

# findLocalAlignments---------------------------------------------------------------------
def test_findLocalAlignments(capsys):
    assert findLocalAlignments(
        ("TTTACATATCGGTGTC", "ACGCG", 2, 2, 1)) == (6, [(8, 1, "ATCG-G", "A-CGCG")])
    
    out, err = capsys.readouterr()
    assert out == err == ""

def test_findLocalAlignmentsPrints(capsys):
    assert findLocalAlignments(("TTTACATATCGGTGTC", "ACGCG", 2, 2, 1),
        doLogProgress = True) == (6, [(8, 1, "ATCG-G", "A-CGCG")])
    
    out, err = capsys.readouterr()
    assert err == ""
    assert out == "Filling score and directions matrices...\nReconstructing best local alignments...\n"