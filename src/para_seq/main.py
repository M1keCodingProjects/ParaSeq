## Main application file, run this if starting the project manually from an editor.
from para_seq.input_manager   import setupArgParser, parseInputArgs
from para_seq.local_alignment import findLocalAlignments
from para_seq.output_manager  import displayOutputSummary, saveOutput

def main(args :tuple[str, ...]|None = None, *, isDebugMode = False) -> None:
    """
    Main application entry point.
    
    Args:
        args (tuple[str, ...] | None): The input arguments, if passed manually for testing purposes. Defaults to: None.
        isDebugMode (bool, optional): If True, the method prints additional information for debugging purposes. Defaults to: False.
    """
    print("Starting analysis...")
    args = setupArgParser().parse_args(args)

    print("Retrieving sequences...")
    *analysisParams, outputPath, shownAlignments, maxSeqLen = parseInputArgs(args)

    maxScore, bestLocalAlignments = findLocalAlignments(
        analysisParams, doLogProgress = True, doShowMatrices = isDebugMode)
    
    if not bestLocalAlignments:
        print("No alignments were found, which might indicate that your sequences' \
nucleotides are completely different.")
        return
    
    displayOutputSummary(maxScore, bestLocalAlignments, shownAlignments, maxSeqLen)
    saveOutput(outputPath, maxScore, bestLocalAlignments)

    print(f"All done! Check the full list of alignments at \"{outputPath}\".")

# Why here? Because I want to be able to test main and catch specific errors. Meanwhile
# a user running the project doesn't want their terminal polluted with the whole stack
# trace, but rather for the program to exit cleanly and explain what happened.
if __name__ == "__main__":
    try: main()
    except BaseException as err:
        msg = "The analysis was interrupted due to an error"
        msg += f", exit code: {err.code}." if isinstance(err, SystemExit) else f":\n{err}"
        print(msg)