## Analysis pipeline module
from numpy           import ndarray, uint8, uint32, dtype, arange, column_stack, argwhere
from para_seq        import DIRS_MATRIX_SHMEM_NAME, SCORE_MATRIX_SHMEM_NAME
from multiprocessing import Pool

from multiprocessing.shared_memory import SharedMemory

type Alignment      = tuple[int, int, str, str]
type AnalysisParams = tuple[str, str, int, int, int]

# The 3 possible backtracking dirs are encoded as single bits of different value, such
# that a single bitflag can hold all combinations:
UP_DIR, DIAG_DIR, LEFT_DIR = 1, 2, 4

# Structuring the values needed by all processes as global consts allows me to set them
# during Pool init, greatly reducing the amount of args I need to pass to each process
# and avoiding the creation of long lists of the same values copied over and over
TARGET_SEQ       = ""
QUERY_SEQ        = ""
MATCH_SCORE      = 0
MISMATCH_PENALTY = 0
GAP_PENALTY      = 0
MATRIX_SHAPE     = (0, 0)
def _setProcessTaskConsts(analysisParams:AnalysisParams) -> None:
    """
    Sets values for unchanging analysis parameters as global constants, also computing
    matrix shape from the provided sequences. Meant as an initializer for pooled
    processes, can also be called whenever to reduce the amount of args that need to be
    passed around.

    Args:
        analysisParams (AnalysisParams): A tuple containing:
        - targetSeq (str) : The target sequence to align.
        - querySeq (str) : The query sequence to align.
        - matchScore (int) : The alignment score bonus for a nucleotide match.
        - mismatchPenalty (int) : The alignment score penalty for a nucleotide mismatch.
        - gapPenalty (int) : The alignment score gap penalty for gap opening and extension.
    """
    global TARGET_SEQ, QUERY_SEQ, MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY, MATRIX_SHAPE
    TARGET_SEQ, QUERY_SEQ, MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY = analysisParams
    MATRIX_SHAPE = getMatrixShape(TARGET_SEQ, QUERY_SEQ)

# Untested, as it would be a very convoluted setup. Sufficient test coverage on the
# process-joining functions should be enough to test this as well.
def computeCellScoreAndDirs(x:int, y:int) -> int:
    """
    **Only works as process task**\n
    Computes alignment score and backtracking directions for cell at the provided
    coordinates.

    Args:
        x (int): The x coordinate of the cell, corresponding to a nucleotide in the target sequence (or a gap if 0).
        y (int): The y coordinate of the cell, corresponding to a nucleotide in the query sequence (or a gap if 0).
    
    **Side effects**
        scoreMatrix: mutates
        dirsMatrix: mutates
        thread safe
    
    Returns:
        int: The computed alignment score for this cell.
    """
    global UP_DIR, DIAG_DIR, LEFT_DIR
    global MATRIX_SHAPE, MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY, QUERY_SEQ, TARGET_SEQ

    # The whole thing is 0-init so we just skip the first row/column cells
    if not y or not x: return 0
    
    scoreMatrix, scoreSharedMem, dirsMatrix, dirsSharedMem = createMatrices(
        MATRIX_SHAPE, isNew = False)
    
    # vvv int casting prevents underflow errors
    insertion  = int(scoreMatrix[y    , x - 1]) - GAP_PENALTY
    deletion   = int(scoreMatrix[y - 1, x    ]) - GAP_PENALTY
    comparison = int(scoreMatrix[y - 1, x - 1]) + (
        MATCH_SCORE if QUERY_SEQ[y - 1] == TARGET_SEQ[x - 1] else -MISMATCH_PENALTY)
    # ^^^ -1 on seq pos is due to the matrix having an extra row/column for gaps.

    scoreMatrix[y, x] = score = max(0, comparison, deletion, insertion)
    if score: dirsMatrix[y, x] = (
        (score == deletion)   * UP_DIR   |
        (score == comparison) * DIAG_DIR |
        (score == insertion)  * LEFT_DIR)

    freeSharedMem(scoreSharedMem)
    freeSharedMem(dirsSharedMem)
    return score

def getMatrixShape(targetSeq:str, querySeq:str) -> tuple[int, int]:
    """
    Computes shape (rows, columns) of the alignment score and directions matrices based on
    the sequences to align.

    Args:
        targetSeq (str): The target sequence, defining the amount of columns.
        querySeq (str): The query sequence, defining the amount of rows.

    Returns:
        tuple:
        - int: The amount of rows.
        - int: The amount of columns.
    """
    return (len(querySeq) + 1, len(targetSeq) + 1)

# Helper method since we need to repeat this bit of code in main and in the processes:
def createMatrices(shape:tuple[int, int], *, isNew = True) -> tuple[ndarray, SharedMemory, ndarray, SharedMemory]:
    """
    Creates or retrieves reference to alignment score and directions matrices with
    provided shape. Remember to call freeSharedMem on the SharedMemory instance at the end
    of every process.
    Args:
        shape (tuple[int, int]): The dimensions (rows, columns) of the matrices.
        isNew (bool, optional): Whether to create (True) the matrix or simply retrieve it (False). Defaults to True.

    Returns:
        tuple:
        -np.ndarray: The created/retrieved alignment score matrix.
        -SharedMemory: The SharedMemory instance tied to the alignment score matrix.
        -np.ndarray: The created/retrieved directions matrix.
        -SharedMemory: The SharedMemory instance tied to the directions matrix.
    """
    # In local alignment scores can never be negative (therefore uint32)
    return (
        *createSharedMatrix(shape, uint32, SCORE_MATRIX_SHMEM_NAME, isNew = isNew),
        *createSharedMatrix(shape, uint8,  DIRS_MATRIX_SHMEM_NAME,  isNew = isNew))

def freeSharedMem(mem:SharedMemory, *, isFreedCompletely = False) -> None:
    """
    Releases, and ultimately destroys if specified, the provided shared memory region.
    Remember to call this method on all live SharedMemory instances at the end of every
    process that uses them.

    Args:
        mem (SharedMemory): The SharedMemory to release and/or destroy.
        isFreedCompletely (bool, optional): Fully frees the shared memory region if True, only do this in the main process after all child processes have finished executing. Defaults to False.
    """
    mem.close()
    if isFreedCompletely: mem.unlink()

def createSharedMatrix(shape:tuple[int, int], itemType:dtype, name:str, *, isNew = False) -> tuple[ndarray, SharedMemory]:
    """
    Creates or retrieves reference to matrix with provided shape and itemType and bound to
    SharedMemory buffer with provided name. Upon creation the matrix will be filled with
    zeros. Remember to call freeSharedMem on the SharedMemory instance at the end of every
    process.

    Args:
        shape (tuple[int, int]): The dimensions (rows, columns) of the matrix.
        itemType (np.dtype): Type of the values in the matrix, compatible with it holding uint values and starting out filled with zeros.
        name (str): Name for the SharedMemory instance whose buffer is bound to the matrix, necessary to be able to retrieve the same matrix in another process.
        isNew (bool, optional): Whether to create (True) the matrix or simply retrieve it (False). Defaults to False.

    Returns:
        tuple: A tuple containing:
        - np.ndarray: The created/retrieved matrix.
        - SharedMemory: The SharedMemory instance tied to the matrix.
    """
    sharedMem  = SharedMemory(
        name   = name,
        create = isNew,
        size   = shape[0] * shape[1] * dtype(itemType).itemsize)
    #   ^^^ When attaching to an existing shared memory block this param is ignored.
    
    matrix = ndarray(shape, dtype = itemType, buffer = sharedMem.buf)
    if isNew: matrix.fill(0) # I can't use zeros because it creates its own buffer.

    return matrix, sharedMem

def computeAntidiagCoords(antidiagId:int, rowsAmt:int, columnsAmt:int) -> ndarray:
    """
    Computes the coordinates of cells belonging to the antidiagonal at the provided index
    in a matrix with the provided amount of rows and columns.

    Args:
        antidiagId (int): The index of the antidiagonal of which to compute the cell coordinates.
        rowsAmt (int): The amount of rows in the matrix.
        columnsAmt (int): The amount of columns in the matrix.

    Returns:
        np.ndarray: 2D-array collection of cell coordinates, in order, belonging to the specified antidiagonal. Each coordinate pair is a row of this 2D array.
    """
    ys = arange(
        max(0, antidiagId - columnsAmt + 1),
        min(antidiagId, rowsAmt - 1) + 1, # The +1 accounts for upper bound exclusion
        dtype = uint32) # Non-neg coords can be uint
    
    # Each i-th antidiag is composed of the cells whose coords sum to i
    xs = antidiagId - ys
    return column_stack((xs, ys))

def fillMatrices(analysisParams:AnalysisParams) -> int:
    """
    Fill aligment score and directions matrices based on the provided analysis parameters.

    Args:
        analysisParams (AnalysisParams): A tuple containing:
        - targetSeq (str) : The target sequence to align.
        - querySeq (str) : The query sequence to align.
        - matchScore (int) : The alignment score bonus for a nucleotide match.
        - mismatchPenalty (int) : The alignment score penalty for a nucleotide mismatch.
        - gapPenalty (int) : The alignment score gap penalty for gap opening and extension.
    
    Returns:
        int: The maximum alignment score found in the score matrix. All the cells with this value are the starting point for the backtracking step.
    """
    maxScore = 0
    # Recomputing this a lot is not a problem since it's a simple operation and it helps
    # isolate the function for testing:
    rowsAmt, columnsAmt = getMatrixShape(*analysisParams[:2])
    # Python automatically spreads initargs into the initializer, so I need to wrap them:
    with Pool(initializer = _setProcessTaskConsts, initargs = (analysisParams,)) as pool:
        for antidiagId in range(rowsAmt + columnsAmt - 1):
            # Each cell in the same antidiag can be computed in parallel:
            antidiag = computeAntidiagCoords(antidiagId, rowsAmt, columnsAmt)
            antidiagMaxScore = max(pool.starmap(computeCellScoreAndDirs, antidiag))
            if maxScore < antidiagMaxScore: maxScore = antidiagMaxScore

    return maxScore

def reconstructAlignments(scoreMatrix:ndarray, maxScore:int, analysisParams:AnalysisParams) -> list[Alignment]:
    """
    Reconstruct all best local alignments based on the filled matrices, the maximum
    alignment score identified and the provided analysis parameters.

    Args:
        scoreMatrix (np.ndarray): The filled alignment score matrix.
        maxScore (int): The maximum alignment score found in the score matrix. All the cells with this value are the starting point for the backtracking step.
        analysisParams (AnalysisParams): A tuple containing:
            - targetSeq (str) : The target sequence to align.
            - querySeq (str) : The query sequence to align.
            - matchScore (int) : The alignment score bonus for a nucleotide match.
            - mismatchPenalty (int) : The alignment score penalty for a nucleotide mismatch.
            - gapPenalty (int) : The alignment score gap penalty for gap opening and extension.

    Returns:
        list[Alignment]: All the optimal local alignments, ignoring exact duplicates.
    """
    if not maxScore: return [] # No point in aligning if maxScore is 0

    bestLocalAlignments :set[Alignment] = set()
    # vvv Python automatically spreads initargs into the initializer, so I need to wrap them:
    with Pool(initializer = _setProcessTaskConsts, initargs = (analysisParams,)) as pool:
        for alignments in pool.starmap(execTraceback, argwhere(scoreMatrix == maxScore)):
            bestLocalAlignments.update(alignments)

    # List conversion is useful for slicing this collection and everything else we want
    # to do in output.
    return list(bestLocalAlignments)

# Untested, as it would be a very convoluted setup. Sufficient test coverage on the
# process-joining functions should be enough to test this as well.
# Coords here are accepted y first to comply with numpy.
def execTraceback(startY:int, startX:int) -> list[Alignment]:
    """
    **Only works as process task**\n
    Executes traceback and reconstructs local alignments, starting from cell at the
    provided coordinates.

    Args:
        y (int): The y coordinate of the local alignment starting cell.
        x (int): The x coordinate of the local alignment starting cell.
    
    **Side effects**
        scoreMatrix: mutates
        dirsMatrix: mutates
        thread safe
    
    Returns:
        list[Alignment]: All the best local alignments.
    """
    global UP_DIR, DIAG_DIR, LEFT_DIR, MATRIX_SHAPE, QUERY_SEQ, TARGET_SEQ

    scoreMatrix, scoreSharedMem, dirsMatrix, dirsSharedMem = createMatrices(
        MATRIX_SHAPE, isNew = False)
    
    bestAlignments :list[Alignment] = []
    
    stack = []
    stack.append((startX, startY, "", ""))
    while stack:
        x, y, targetAlignment, queryAlignment = stack.pop()

        # Local alignment ends at any cell with a value of 0:
        if not scoreMatrix[y, x]:
            if targetAlignment and queryAlignment: # No point in saving empty alignments
                # Coords are shifted by 1 to enter a 1-based system of reference:
                bestAlignments.append((x + 1, y + 1, targetAlignment, queryAlignment))
            
            continue

        # Each fork in the road adds a stack entry that sends the program to the next cell
        # and keeps track of the aligned seqs so far:
        cellDirs = dirsMatrix[y, x]
        if cellDirs & UP_DIR: stack.append((
            x, y - 1, '-' + targetAlignment, QUERY_SEQ[y - 1] + queryAlignment))
        
        if cellDirs & DIAG_DIR: stack.append((
            x - 1, y - 1, TARGET_SEQ[x - 1] + targetAlignment, QUERY_SEQ[y - 1] + queryAlignment))

        if cellDirs & LEFT_DIR: stack.append((
            x - 1, y, TARGET_SEQ[x - 1] + targetAlignment, '-' + queryAlignment))

    freeSharedMem(scoreSharedMem)
    freeSharedMem(dirsSharedMem)
    return bestAlignments

# Contains some prints since it's intended as the main collection of analysis pipeline
# steps, to be called in the main file:
def findLocalAlignments(analysisParams:AnalysisParams, *, doLogProgress = False, doShowMatrices = False) -> tuple[int, list[Alignment]]:
    """
    Find all local alignments starting from the provided analysis parameters.

    Args:
        analysisParams (AnalysisParams): A tuple containing:
            - targetSeq (str) : The target sequence to align.
            - querySeq (str) : The query sequence to align.
            - matchScore (int) : The alignment score bonus for a nucleotide match.
            - mismatchPenalty (int) : The alignment score penalty for a nucleotide mismatch.
            - gapPenalty (int) : The alignment score gap penalty for gap opening and extension.

        doLogProgress (bool, optional): If True prints analysis progress messages to standard output. Defaults to: False.
        doShowMatrices (bool, optional): If True prints the filled score and directions matrices to standard output, useful for debugging. Defaults to: False.
    
    Returns:
        tuple: The maximum alignment score and all the local alignments, ignoring exact duplicates.
    """
    _setProcessTaskConsts(analysisParams)
    scoreMatrix, scoreSharedMem, dirsMatrix, dirsSharedMem = createMatrices(MATRIX_SHAPE)

    if doLogProgress: print("Filling score and directions matrices...")
    maxScore = fillMatrices(analysisParams)

    if doShowMatrices:
        print("score matrix:", scoreMatrix, "directions matrix:", dirsMatrix,
              sep = "\n\n", end = "\n\n")

    if doLogProgress: print("Reconstructing best local alignments...")
    bestLocalAlignments = reconstructAlignments(scoreMatrix, maxScore, analysisParams)

    freeSharedMem(scoreSharedMem, isFreedCompletely = True)
    freeSharedMem(dirsSharedMem,  isFreedCompletely = True)
    return maxScore, bestLocalAlignments

# The main is used here to showcase how to use this file's functions:
def main() -> None:
    targetSeq = "GTCACCGTAGGATCGATGCTTAGCTACGATCGATCGATCGTAGCTAGCTAGCTAGTCGATCGATCGATAGCTAGCTAGCTAGCTAGTACGTAGCTAGCTGATCGATCGTGACGTAGC"
    querySeq  = "TTGACCGTAGGATAGTCGATCGATCGATAGCTAGCTAGCTAGCTAGTACGATAGCTTTCGATAGCTAGCATGCTAGC"
    maxScore, bestLocalAlignments = findLocalAlignments((targetSeq, querySeq, 2, 1, 2))
    print(maxScore)
    print(*bestLocalAlignments, sep = '\n')

if __name__ == "__main__": main()