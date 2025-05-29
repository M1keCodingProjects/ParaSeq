import pytest
from src.para_seq.main import *
from src.para_seq import ALIGNMENT_INFO

# main------------------------------------------------------------------------------------
@pytest.mark.parametrize("args", [
    ("AAAAAAAAAAAAAAAAAA", "TTT", "-m" '2', "-mm", '2', "-g", '1', "-ma", '1'),
    ("AAAAAAAAAAAAAAAAAA", "AAA", "-m" '0', "-mm", '0', "-g", '0', "-ma", '1')])
def test_mainNoAlignments(capsys, args):
    main(args)
    
    out, err = capsys.readouterr()
    assert err == ""
    assert out == """Starting analysis...
Retrieving sequences...
Filling score and directions matrices...
Reconstructing best local alignments...
No alignments were found, which might indicate that your sequences' nucleotides are completely different.
"""

def test_main(capsys):
    outputStr  = """Starting analysis...
Retrieving sequences...
Filling score and directions matrices...
Reconstructing best local alignments...
Best local alignment score: 4
"""
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