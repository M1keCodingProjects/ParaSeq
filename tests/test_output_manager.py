import pytest
from para_seq.output_manager import *
from para_seq import MAX_DISPLAYED_SEQ_LEN, MAX_DISPLAYED_ALIGNMENTS

# displayOutputSummary--------------------------------------------------------------------
def test_displayOutputSummary(capsys):
    # Nothing gets cut off here:
    displayOutputSummary(
        10, [(0, 3, "ATT-C", "A-TTC"), (4, 5, "ACG", "GGTG")],
        MAX_DISPLAYED_ALIGNMENTS, MAX_DISPLAYED_SEQ_LEN)
    
    out, err = capsys.readouterr()
    assert err == ""
    assert out == ("Best local alignment score: 10\n" +
        ALIGNMENT_INFO.format(0, 3, "ATT-C", "A-TTC") +
        ALIGNMENT_INFO.format(4, 5, "ACG", "GGTG") + '\n')

def test_displayOutputSummaryCutoff(capsys):
    displayOutputSummary(
        10, [(0, 3, "ATT-CATT-CATT-C", "A-TATT-CTATT-C-C"), (4, 5, "ACG", "GGTG")], 1, 10)

    out, err = capsys.readouterr()
    assert err == ""
    assert out == ("Best local alignment score: 10\n" +
        ALIGNMENT_INFO.format(0, 3, "ATT-CAT...", "A-TATT-...") + '\n')

# Shouldn't happen, as empty outputs are caught beforehand:
def test_displayOutputEmpty(capsys):
    displayOutputSummary(0, [], 0, 0)
    out, err = capsys.readouterr()
    assert err == ""
    assert out == "Best local alignment score: 0\n\n"

# saveOutput------------------------------------------------------------------------------
# Shouldn't happen, as empty outputs are caught beforehand:
def test_saveOutputEmpty(tmp_path):
    path = tmp_path / "output.txt"
    saveOutput(path, 0, [])
    with open(path) as fd:
        assert fd.read() == "Score: 0\nTotal alignments: 0\n"

def test_saveOutput(tmp_path):
    path = tmp_path / "output.txt"
    saveOutput(path, 10, [
        (0, 3, "ACT", "A--"), (3, 2, "ATTAT", "CGCATAT"), (1, 1, "TT-TTT-T", "-T")])

    with open(path) as fd:
        assert fd.read() == """Score: 10
Total alignments: 3

Target start pos: 0
Query start pos: 3
Target sequence: ACT
Query sequence:  A--

Target start pos: 3
Query start pos: 2
Target sequence: ATTAT
Query sequence:  CGCATAT

Target start pos: 1
Query start pos: 1
Target sequence: TT-TTT-T
Query sequence:  -T
"""

# Shouldn't happen, as invalid file paths are caught beforehand:
def test_saveOutputInvalidPath():
    with pytest.raises(PermissionError) as errInfo: saveOutput("./output/", 0, [])
    assert str(errInfo.value) == "[Errno 13] Permission denied: './output/'"