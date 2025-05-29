## Generic utilities
def ellipsize(s:str, size:int) -> str:
    """
    Truncates to size and adds ellipsis (...) to the end of the string, if needed.

    Args:
        s (str): The provided string to truncate.
        size (int): The max accepted size for the provided string for which it won't be truncated.

    Returns:
        str: The possibly truncated and ellipsized string.
    """
    return s[:size - 3] + "..." if len(s) > size else s

class CustomErr(Exception):
    """General custom error template class."""
    msgPrefix = ""
    def __init__(self, msg:str, details = "") -> None:
        """
        Create a CustomErr object with a specific message and, optionally, more details
        about the error.

        Args:
            msg (str): The error message.
            details (str, optional): More details about the error. Defaults to "".
        """
        super().__init__(f"{self.msgPrefix}: {msg}" + f", {details}" * bool(details) + '.')
        # ^^^ details are added only if present.
        self.msg     = msg
        self.details = details