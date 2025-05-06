from src.para_seq.__init__ import *

# SeqName---------------------------------------------------------------------------------
def test_SeqName():
    assert SeqName.Target.id == 't'
    assert SeqName.Query.id  == 'q'