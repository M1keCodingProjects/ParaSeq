## Output manager module
from para_seq                 import ALIGNMENT_INFO
from para_seq.utils           import ellipsize
from para_seq.local_alignment import Alignment

def displayOutputSummary(maxScore:int, bestLocalAlignments:list[Alignment], maxDisplayedAlignments:int, maxDisplayedSeqLen:int) -> None:
    """
    Prints a summary of the result of the alignment procedure to standard output,
    including the first few alignments up to the provided amount and truncating the
    aligned sequences when longer than the provided "maximum displayed sequence length"
    parameter.

    Args:
        maxScore (int): The maximum alignment score found in the score matrix.
        bestLocalAlignments (list[Alignment]): All the optimal local alignments.
        maxDisplayedAlignments (int): The maximum number of shown alignments in the terminal output.
        maxDisplayedSeqLen (int): The length after which aligned sequences are truncated in the terminal output.
    """
    print("Best local alignment score:", maxScore)
    outputBuf = ""
    for x, y, alignedTarget, alignedQuery in bestLocalAlignments[:maxDisplayedAlignments]:
        outputBuf += ALIGNMENT_INFO.format(x, y,
            ellipsize(alignedTarget, maxDisplayedSeqLen),
            ellipsize(alignedQuery, maxDisplayedSeqLen))
    
    print(outputBuf)

def saveOutput(outputPath:str, maxScore:int, bestLocalAlignments:list[Alignment]) -> None:
    """
    Saves entire result of the alignment procedure to a file at the provided path,
    creating it if it doesn't exist and overwriting it otherwise.

    Args:
        outputPath (str): The path to the output file.
        maxScore (int): The maximum alignment score found in the score matrix.
        bestLocalAlignments (list[Alignment]): All the optimal local alignments.
    """
    outputBuf  = f"Score: {maxScore}\nTotal alignments: {len(bestLocalAlignments)}\n"
    outputBuf += "".join( # ALIGNMENT_INFO already has newlines
        map(lambda alignment: ALIGNMENT_INFO.format(*alignment), bestLocalAlignments))
    
    with open(outputPath, 'w') as fd: fd.write(outputBuf)