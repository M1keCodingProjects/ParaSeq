from src.para_seq.utils import *

# CustomErr-------------------------------------------------------------------------------
class CustomErrExt(CustomErr):
    msgPrefix = "prefix"

def test_CustomErr():
    assert str(CustomErrExt("msg", "details")) == "prefix: msg, details."

def test_CustomErrNoDetails():
    assert str(CustomErrExt("msg")) == "prefix: msg."