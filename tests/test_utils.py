from src.para_seq.utils import *

# ellipsize-------------------------------------------------------------------------------
def test_ellipsize():
    assert ellipsize("longer than 5", 5) == "lo..."

def test_ellipsizeExact():
    assert ellipsize("exactly 9", 9) == "exactly 9"

def test_ellipsizeShorter():
    assert ellipsize("foo", 100) == "foo"

def test_ellipsizeInvalidSize():
    assert ellipsize("", -1) == "..."
    assert ellipsize("foo", -10) == "..."

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