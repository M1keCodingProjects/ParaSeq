from src.para_seq.utils import *

# CustomErr-------------------------------------------------------------------------------
class CustomErrExt(CustomErr):
    msgPrefix = "prefix"

def test_CustomErr():
    err = CustomErrExt("msg", "details")
    assert str(err) == "prefix: msg, details."
    assert err.msg == "msg"
    assert err.details == "details"

def test_CustomErrNoDetails():
    assert str(CustomErrExt("msg")) == "prefix: msg."