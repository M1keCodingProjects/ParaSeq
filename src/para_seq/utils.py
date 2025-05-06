## Generic utilities
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