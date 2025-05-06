## Main application file, run this if starting the project manually from an editor.
from para_seq.io_manager import collectAndParseInputArgs

def main(args :tuple[str, ...]|None = None) -> None:
    """
    Main application entry point.
    
    Args:
        args (tuple[str, ...] | None): The input arguments, if passed manually for testing purposes. Defaults to: None.
    """
    targetSeq, querySeq, matchScore, mismatchPenalty, gapPenalty = collectAndParseInputArgs(args)

if __name__ == "__main__": main()