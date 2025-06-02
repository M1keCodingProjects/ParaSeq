## Input manager module
from os             import SEEK_END
from pyfastx        import Fasta
from para_seq       import *
from argparse       import ArgumentParser, Namespace
from para_seq.utils import CustomErr, ellipsize

# Custom errors:
class IdenticalSeqsErr(CustomErr):
    """Error class for attempted alignment of identical sequences."""
    msgPrefix = IDENTICAL_SEQS_PREFIX

class MissingSeqErr(CustomErr):
    """Error class for missing sequences."""
    msgPrefix = MISSING_SEQ_PREFIX

class InvalidSeqErr(CustomErr):
    """Error class for invalid sequences."""
    msgPrefix = INVALID_SEQ_PREFIX

class InvalidFileErr(CustomErr):
    """Error class for invalid FASTA or output (.txt) files."""
    msgPrefix = INVALID_FILE_PREFIX

# This method can be considered enough to discern between a file and a raw DNA seq.
def isValidFastaFilePath(filePath:str) -> bool:
    """
    Simply checks if the provided string ends with a valid FASTA file extension.

    Args:
        filePath (str): The string to check.

    Returns:
        bool: If the string can be considered a valid FASTA file path or not.
    """
    return filePath.lower().endswith((".fa", ".fasta"))

def isValidDNA(seq:str) -> bool:
    """
    Checks case-insensitively that the provided string is non-empty and made up entirely
    of valid DNA nucleotide characters in the (mostly) unambiguous "ACGTN" alphabet.

    Args:
        seq (str): The provided sequence string.

    Returns:
        bool: If the sequence can be considered valid DNA or not.
    """
    # Adding the lowercase variants to the set instead of calling seq.upper() is perhaps
    # premature and needless optimization, but let's go for it:
    return seq and set(seq) <= {'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'}
    # ^^^ checking the seq itself avoids the set computation for empty strs, which are
    # obviously invalid but also should never happen given how the program works.

# Type casting function passed to some ArgumentParser args
def uint(value:str) -> int:
    """
    Type casting function from string to non-negative integer (unsigned, uint).

    Args:
        value (str): The string representation of a uint value.

    Raises:
        ValueError: When the provided string does not represent a uint value, or
        if said value is written in scientific notation or anything else besides
        fully explicit decimal notation.

    Returns:
        int: The converted value.
    """
    # This immediately excludes the "-" sign and any decimal or scientific notation
    if not value.isdigit(): raise ValueError(UINT_ERR.format(value))
    
    return int(value)

type DNA = str # Valid DNA, all the characters belong to the ACGTN set.
def validateDNA(seq:str) -> DNA:
    """
    Parse raw DNA sequence string into valid DNA string with upper-case nucleotides.

    Args:
        seq (str): The raw DNA sequence string.

    Raises:
        InvalidSeqErr: If the sequence string is not valid DNA (only characters in the ACGTN set).

    Returns:
        DNA: The valid DNA sequence.
    """
    if not isValidDNA(seq): raise InvalidSeqErr("")
    # The specifics of the error are filled in later, which isn't ideal in general but
    # the error class does a good job of explaining what more or less happened.

    return seq.upper()

# Untested as it requires a mock file descriptor which is a whole mess, I verified
# manually that the newline is being added.
# Add a trailing newline if needed: this solves a bug in pyfastx
def addTrailingNewline(filePath:str) -> None:
    """
    Adds trailing empty line to the (assumed FASTA) file at the provided path, if missing.

    Args:
        filePath (str): The provided file path.

    Raises:
        MissingSeqErr: When the provided FASTA file is empty.
        InvalidFileErr: When the provided file path doesn't lead to an existent FASTA file.
    """
    try:
        # Reading in bytes makes the process more streamlined:
        with open(filePath, 'rb+') as f:
            f.seek(0, SEEK_END) # Move cursor to EOF
            if f.tell() == 0:   # If we're still at the start the file is empty
                raise MissingSeqErr(f"provided FASTA file \"{filePath}\" is empty")

            # Go back 1 byte
            f.seek(-1, SEEK_END)
            if f.read(1) != b'\n': f.write(b'\n')
    
    except FileNotFoundError as err:
        raise InvalidFileErr(err, "the provided FASTA file doesn't exist")

# Untested as it's hard to isolate (I'd have to create a Fasta instance)
def _getValidSeqFromCollection(collection:Fasta, pos:int, filePath:str) -> DNA:
    """
    Retrieve and validate the sequence at a given position in the provided collection, if
    present.

    Args:
        collection (Fasta): The collection containing the desired sequence.
        pos (int): The 1-based position of the desired sequence in the collection.
        filePath (str): The FASTA file path the collection was loaded from, for error reporting reasons.

    Raises:
        InvalidSeqErr: If the sequence string is not valid DNA (only characters in the ACGTN set).
        MissingSeqErr: If the provided collection does not contain a sequence at the desired position.
    
    Returns:
        DNA: The valid DNA sequence.
    """
    try: return validateDNA(collection[pos - 1].seq)
    except InvalidSeqErr: raise InvalidSeqErr(
        f"sequence at position {pos} in FASTA file \"{filePath}\" is not valid DNA")
    except IndexError: raise MissingSeqErr(
        f"FASTA file \"{filePath}\" doesn't contain a sequence at position {pos}")

# Obtain seq object from FASTA file path, tests on this also cover the func above
def parseFastaSeq(filePath:str, seq1Pos:int, seq2Pos :int = None) -> tuple[DNA, DNA]:
    """
    Attempts to load a FASTA file at the provided path, then validates, parses and
    retrieves the sequence(s) at the desired position(s).

    Args:
        filePath (str): The provided FASTA file path.
        seq1Pos (int): 1-based position of the desired sequence.
        seq2Pos (int, optional): 1-based position of a second desired sequence, if provided. Defaults to: None.
    
    Raises:
        InvalidFileErr: When the provided FASTA file is malformed.
    
    Returns:
        tuple: The first, and optionally second, sequence(s) at the desired position(s).
    """
    addTrailingNewline(filePath)

    try: seqs = Fasta(filePath)
    except RuntimeError as err: raise InvalidFileErr(err, "file is malformed")

    seq1 = _getValidSeqFromCollection(seqs, seq1Pos, filePath)
    seq2 = "" if seq2Pos is None else _getValidSeqFromCollection(seqs, seq2Pos, filePath)
    return seq1, seq2

# A method that checks which type of seq info was given and acts accordingly:
def parseSeq(seqOrFilePath:str, pos:int, name:SeqName) -> DNA:
    """
    Obtains valid DNA sequence from raw DNA sequence string or FASTA file path.

    Args:
        seqOrFilePath (str): The raw DNA sequence string or FASTA file path.
        pos (int): The 1-based position of the sequence in the FASTA file, if provided.
        name (SeqName): Sequence metadata merely for error reporting purposes, used if the sequence is provided as raw DNA.

    Raises:
        InvalidSeqErr: If the provided string could not be interpreted as a FASTA file but is also not valid DNA (only characters in the ACGTN set).

    Returns:
        DNA: The valid DNA sequence.
    """
    if isValidFastaFilePath(seqOrFilePath):
        return parseFastaSeq(seqOrFilePath, pos)[0] # <-- We only take one seq here
    
    try: return validateDNA(seqOrFilePath)
    except InvalidSeqErr:
        raise InvalidSeqErr(f"Your \"{name}\" sequence is not valid DNA")
    
    # ^^^ Not as dumb as it looks, I didn't want to pass useless info to the validateDNA
    # function, the error stays the same but gets enriched with the seq name/path here.

# Helper to isolate the "same file" case a bit more:
def parseSeqsFromFile(filePath:str, targetPos:int, queryPos:int) -> tuple[DNA, DNA]:
    """
    Obtain valid DNA target and query sequences, taken from the same FASTA file.

    Args:
        filePath (str): The provided FASTA file path.
        targetPos (int): The 1-based position of the target sequence in the FASTA file.
        queryPos (int): The 1-based position of the query sequence in the FASTA file.

    Raises:
        IdenticalSeqsErr: If the sequences to load are the same.

    Returns:
        tuple[DNA, DNA]: Respectively the target and query sequences.
    """
    # It makes no sense to align a sequence with itself:
    if queryPos == targetPos:
        if targetPos != 1: raise IdenticalSeqsErr(
            f"tried loading sequences coming from the same file ({filePath})",
            f"with identical positions ({targetPos})")

        # If both positions are at 1 then they haven't been set (or the user feels very
        # funny), in this "one file" case the target and query seqs will be respectively
        # the first 2 seqs of the file:
        queryPos = 2
    
    return parseFastaSeq(filePath, targetPos, queryPos)

def setupArgParser() -> ArgumentParser:
    """
    Setup an argparse.ArgumentParser instance and obtain the values for the input
    arguments needed for the analysis.

    Returns:
        Namespace: Object containing the arguments and their values as properties:
            - .target_seq (str): Target DNA sequence or FASTA file path.
            - .query_seq (str): Query DNA sequence or FASTA file path, if provided.
            - .target_pos (int): 1-based target DNA sequence position in FASTA file.
            - .query_pos (int): 1-based query DNA sequence position in FASTA file.
            - .match_score (int): Match score.
            - .mismatch_penalty (int): Mismatch penalty.
            - .gap_penalty (int): Constant gap penalty.
            - .output_path (str): Path to output file.
            - .max_alignments_shown (int): Maximum number of alignments shown before terminal output is cut off.
            - .longest_sequence_shown (int): Maximum aligned sequence length before output is truncated.
        
    All positions, scores and penalties are non-negative.
    """
    parser = ArgumentParser(prog = PACKAGE_NAME, description = PACKAGE_DESCR)

    # The user can specify a single seq arg, as long as it's a FASTA file path
    parser.add_argument("target_seq", type = str, help = TARGET_SEQ_HELP)
    parser.add_argument("query_seq",
        type = str, nargs = '?', default = "", help = QUERY_SEQ_HELP)
    # ^^^ The use of _ here is necessary, as for some reason it doesn't convert
    # - to _ in pos args..

    # When a file containing the desired seq(s) is provided, the user can specify
    # which specific sequence of that file he wants. Positions are not allowed to be
    # negative and are interpreted as 1-based which is more user-friendly.
    parser.add_argument("--target-pos", "-tp",
        type = uint, default = 1, help = TARGET_POS_HELP)
    
    parser.add_argument("--query-pos", "-qp",
        type = uint, default = 1, help = QUERY_POS_HELP)

    # As is customary for simple DNA alignment, all the relevant scores are
    # (non-negative) integers:
    parser.add_argument("--match-score", "-m",
        type = uint, required = True, help = MATCH_HELP)
    
    parser.add_argument("--mismatch-penalty", "-mm",
        type = uint, required = True, help = MISMATCH_HELP)
    
    parser.add_argument("--gap-penalty", "-g",
        type = uint, required = True, help = GAP_HELP)

    # Output customization:
    parser.add_argument("--output-path", "-o",
        type = str, default = "./output/output.txt", help = OUT_PATH_HELP)  

    parser.add_argument("--max-alignments-shown", "-ma",
        type = uint, default = MAX_DISPLAYED_ALIGNMENTS, help = MAX_ALIGN_HELP)
    
    parser.add_argument("--longest-sequence-shown", "-ls",
        type = uint, default = MAX_DISPLAYED_SEQ_LEN, help = MAX_SEQ_LEN_HELP)

    return parser

def parseInputArgs(args:Namespace) -> tuple[DNA, DNA, int, int, int, str, int, int]:
    """
    Parse all the CLI input arguments passed by the user and necessary for
    the analysis.

    Args:
        args (Namespace): The Namespace object containing the arguments and their values.

    Raises:
        InvalidFileErr: If the output file path argument is invalid.
        MissingSeqsErr: When not enough information was provided to retrieve the two sequences.
        IdenticalSeqsErr: If the sequences to load are the same.

    Returns:
        tuple:
        - DNA: The target sequence.
        - DNA: The query sequence.
        - int: The match score.
        - int: The mismatch penalty.
        - int: The gap penalty.
        - str: The path to the output file.
        - int: The maximum number of shown alignments in the terminal output.
        - int: The length after which aligned sequences are truncated in the terminal output.
    """
    if not args.target_pos: raise MissingSeqErr(
        "Target sequence position cannot be 0", "provide a positive and valid position")
    
    if not args.query_pos: raise MissingSeqErr(
        "Query sequence position cannot be 0", "provide a positive and valid position")

    try:
        # This is much more straightforward and all-encompassing than a bunch of os checks
        with open(args.output_path, 'w'): pass
    
    except Exception as err: raise InvalidFileErr(
        err, f"\"{args.output_path}\" is not a valid output file path")

    # It's impossible for this check to fail when query doesn't exist AND the 2 seqs are
    # the same, as target must exists:
    areTheSame = args.query_seq == args.target_seq
    if args.query_seq and not areTheSame:
        return (
            parseSeq(args.target_seq, args.target_pos, SeqName.Target),
            parseSeq(args.query_seq,  args.query_pos,  SeqName.Query),
            args.match_score, args.mismatch_penalty, args.gap_penalty,
            args.output_path, args.max_alignments_shown, args.longest_sequence_shown)

    # If only the target seq exists that's fine, as long as it's a file from which to load
    # both seqs:
    if not isValidFastaFilePath(args.target_seq):
        # This means that 2 raw DNA seqs are the same:
        if areTheSame: raise IdenticalSeqsErr(
            "tried loading identical DNA sequences",
            "the 2 provided sequences were interpreted as DNA")

        # This means that we only have a single raw DNA seq:
        raise MissingSeqErr("a single DNA sequence was provided")

    # the following ignores the case where the 2 paths are the same, and handles the
    # case where the query seq doesn't exist:
    return (*parseSeqsFromFile(args.target_seq, args.target_pos, args.query_pos),
        args.match_score, args.mismatch_penalty, args.gap_penalty, args.output_path,
        args.max_alignments_shown, args.longest_sequence_shown)

# The main is used here to showcase how to use this file's functions:
if __name__ == "__main__":
    args = ('data\\good.fasta', 'GGC', "-m", '0', "-mm", '0', "-g", '1')
    args = setupArgParser().parse_args(args)
    seq1, seq2, *rest = parseInputArgs(args)
    print(ellipsize(seq1, MAX_DISPLAYED_SEQ_LEN), ellipsize(seq2, MAX_DISPLAYED_SEQ_LEN), *rest)