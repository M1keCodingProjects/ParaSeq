from src.para_seq.main import *

# main------------------------------------------------------------------------------------
def test_main():
    assert main(("ACGT", ".\\data\\good.fasta", "-m" '2', "-mm", '2', "-g", "1")) == None