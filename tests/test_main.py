import pytest
from src.para_seq.main import *
from src.para_seq import ALIGNMENT_INFO

# main------------------------------------------------------------------------------------
OUT_INTRO = """Starting analysis...
Retrieving sequences...
Filling score and directions matrices...
"""

@pytest.mark.parametrize("args", [
    ("AAAAAAAAAAAAAAAAAA", "TTT", "-m" '2', "-mm", '2', "-g", '1', "-ma", '1'),
    ("AAAAAAAAAAAAAAAAAA", "AAA", "-m" '0', "-mm", '0', "-g", '0', "-ma", '1')])
def test_mainNoAlignments(capsys, args):
    main(args)
    
    out, err = capsys.readouterr()
    assert err == ""
    assert out == OUT_INTRO + "Reconstructing best local alignments...\nNo alignments were found, which might indicate that your sequences' nucleotides are completely different.\n"
    #                          ^^^ This part is separate from OUT_INTRO in case debug mode
    #                              is on, which would print the matrices before this line

def test_main(capsys):
    outputStr  = OUT_INTRO + "Reconstructing best local alignments...\nBest local alignment score: 4\n"
    outputStr += ALIGNMENT_INFO.format("{}", 8, "TT", "TT") + '\n'
    outputStr += "All done! Check the full list of alignments at \"./output/output.txt\".\n"
    # ^^^ I avoid formatting the target alignment start pos because the alignments are
    # collected in a set, which means their insertion order is random and when the first
    # few are selected for output this param (in this example) can be either 0 or 1.

    main(("TTT", "AAAAAAATTCAAA", "-m" '2', "-mm", '2', "-g", '1', "-ma", '1'))
    out, err = capsys.readouterr()
    assert err == ""
    assert out == outputStr.format(1) or out == outputStr.format(2)

    # The file output is tested elsewhere.

# Whole tool tests: these tests precisely check the results of the tool against some local
# alignment problems solved by hand:
def test_example1(capsys):
    main(("ACGGTC", "TGGATCTCCAACG", "-m" '2', "-mm", '2', "-g", '1'), isDebugMode = True)

    matricesOutput = """score matrix:

[[0 0 0 0 0 0 0]
 [0 0 0 0 0 2 1]
 [0 0 0 2 2 1 0]
 [0 0 0 2 4 3 2]
 [0 2 1 1 3 2 1]
 [0 1 0 0 2 5 4]
 [0 0 3 2 1 4 7]
 [0 0 2 1 0 3 6]
 [0 0 2 1 0 2 5]
 [0 0 2 1 0 1 4]
 [0 2 1 0 0 0 3]
 [0 2 1 0 0 0 2]
 [0 1 4 3 2 1 2]
 [0 0 3 6 5 4 3]]

directions matrix:

[[0 0 0 0 0 0 0]
 [0 0 0 0 0 2 4]
 [0 0 0 2 2 5 0]
 [0 0 0 2 2 4 4]
 [0 2 4 1 1 7 7]
 [0 1 0 0 1 2 4]
 [0 0 2 4 5 1 2]
 [0 0 1 7 0 3 1]
 [0 0 2 4 0 1 3]
 [0 0 2 4 0 1 3]
 [0 2 5 0 0 0 1]
 [0 2 4 0 0 0 1]
 [0 1 2 4 4 4 2]
 [0 0 1 2 6 4 4]]
"""
    matricesOutput = matricesOutput.split('\n')

    out, err = capsys.readouterr()
    assert err == ""
    assert out.startswith(OUT_INTRO)
    outLines = out.split('\n')[3:] # Skipping OUT_INTRO
    for i, line in enumerate(matricesOutput):
        assert line == outLines[i], f"failed for line {i}"
    
    pos = len(matricesOutput)
    assert outLines[pos] == "Reconstructing best local alignments..."
    assert outLines[pos + 1] == "Best local alignment score: 7"

    assert ALIGNMENT_INFO.format(3, 2, "GG-TC", "GGATC") in out
    assert out.endswith("All done! Check the full list of alignments at \"./output/output.txt\".\n")

def test_example2(capsys):
    main(("CTG", "CTTGTGCTTGGGACTAAAGACTAAAGCTTGCATG", "-m" '3', "-mm", '3', "-g", '1'), isDebugMode = True)

    matricesOutput = """score matrix:

[[0 0 0 0]
 [0 3 2 1]
 [0 2 6 5]
 [0 1 5 4]
 [0 0 4 8]
 [0 0 3 7]
 [0 0 2 6]
 [0 3 2 5]
 [0 2 6 5]
 [0 1 5 4]
 [0 0 4 8]
 [0 0 3 7]
 [0 0 2 6]
 [0 0 1 5]
 [0 3 2 4]
 [0 2 6 5]
 [0 1 5 4]
 [0 0 4 3]
 [0 0 3 2]
 [0 0 2 6]
 [0 0 1 5]
 [0 3 2 4]
 [0 2 6 5]
 [0 1 5 4]
 [0 0 4 3]
 [0 0 3 2]
 [0 0 2 6]
 [0 3 2 5]
 [0 2 6 5]
 [0 1 5 4]
 [0 0 4 8]
 [0 3 3 7]
 [0 2 2 6]
 [0 1 5 5]
 [0 0 4 8]]

directions matrix:

[[0 0 0 0]
 [0 2 4 4]
 [0 1 2 4]
 [0 1 3 5]
 [0 0 1 2]
 [0 0 3 1]
 [0 0 1 3]
 [0 2 4 1]
 [0 1 2 4]
 [0 1 3 5]
 [0 0 1 2]
 [0 0 1 3]
 [0 0 1 3]
 [0 0 1 1]
 [0 2 4 1]
 [0 1 2 4]
 [0 1 1 5]
 [0 0 1 5]
 [0 0 1 5]
 [0 0 1 2]
 [0 0 1 1]
 [0 2 4 1]
 [0 1 2 4]
 [0 1 1 5]
 [0 0 1 5]
 [0 0 1 5]
 [0 0 1 2]
 [0 2 4 1]
 [0 1 2 4]
 [0 1 3 5]
 [0 0 1 2]
 [0 2 1 1]
 [0 1 1 1]
 [0 1 2 1]
 [0 0 1 2]]
"""
    matricesOutput = matricesOutput.split('\n')

    out, err = capsys.readouterr()
    assert err == ""
    assert out.startswith(OUT_INTRO)

    outLines = out.split('\n')[3:] # Skipping OUT_INTRO
    for i, line in enumerate(matricesOutput):
        assert line == outLines[i], f"failed for line {i}"
    
    pos = len(matricesOutput)
    assert outLines[pos] == "Reconstructing best local alignments..."
    assert outLines[pos + 1] == "Best local alignment score: 8"

    assert ALIGNMENT_INFO.format(1, 31, "C-TG", "CATG") in out
    assert ALIGNMENT_INFO.format(1, 27, "C-TG", "CTTG") in out
    assert ALIGNMENT_INFO.format(1, 27, "CT-G", "CTTG") in out
    assert ALIGNMENT_INFO.format(1,  7, "C-TG", "CTTG") in out
    assert ALIGNMENT_INFO.format(1,  7, "CT-G", "CTTG") in out
    assert ALIGNMENT_INFO.format(1,  1, "C-TG", "CTTG") in out
    assert ALIGNMENT_INFO.format(1,  1, "CT-G", "CTTG") in out
    assert out.endswith("All done! Check the full list of alignments at \"./output/output.txt\".\n")