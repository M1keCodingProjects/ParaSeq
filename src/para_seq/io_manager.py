## Input-Output manager module
from Bio            import SeqIO
from para_seq.utils import CustomErr
from para_seq       import SeqName
from argparse       import ArgumentParser, Namespace
from Bio.SeqRecord  import SeqRecord, Seq

# All the consts I need here...
from para_seq import GAP_HELP, MATCH_HELP, MISMATCH_HELP, PACKAGE_DESCR, PACKAGE_NAME, QUERY_POS_HELP, QUERY_SEQ_HELP, TARGET_POS_HELP, TARGET_SEQ_HELP, UINT_ERR, INVALID_SEQ_PREFIX, MISSING_SEQ_PREFIX

# Custom errors:
class IdenticalSeqsErr(CustomErr):
    """Error class for attempted alignment of identical sequences."""
    msgPrefix = "Alignment of identical sequences is pointless"

class MissingSeqErr(CustomErr):
    """Error class for missing sequences."""
    msgPrefix = MISSING_SEQ_PREFIX

class InvalidSeqErr(CustomErr):
    """Error class for invalid sequences."""
    msgPrefix = INVALID_SEQ_PREFIX

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
    Checks case-insensitively if the provided string is made up entirely of valid DNA
    nucleotide characters in the (mostly) unambiguous "ACGTN" alphabet.

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

# Obtain seq object from raw DNA input:
def parseRawSeq(seq:str, name:SeqName) -> SeqRecord:
    """
    Parses raw DNA sequence string into SeqRecord object.

    Args:
        seq (str): The raw DNA sequence string.
        name (SeqName): Sequence metadata, merely for output purposes.

    Raises:
        InvalidSeqErr: if the sequence string is not valid DNA (only characters in the ACGTN set).
    
    Returns:
        SeqRecord: The parsed DNA sequence, with all the nucleotides set to uppercase.
    """
    # TODO: speed up with numpy
    if not isValidDNA(seq): raise InvalidSeqErr("invalid raw DNA string")

    return SeqRecord(seq = Seq(seq.upper()), id = name.id, name = name.value)

# Obtain seq object from FASTA file path
# Loading the entire file once can improve performance:
def loadFasta(filePath:str) -> list[SeqRecord]:
    """
    Attempts to load a FASTA file at the provided path with SeqIO.

    Args:
        filePath (str): The provided FASTA file path.

    Returns:
        list[SeqRecord]: The collection of sequences in the FASTA file as
        SeqRecord objects.
    """
    return list(SeqIO.parse(filePath, format = "fasta"))

# Even if this was just a wrapper around a single indexing operation it would still be
# worth it for the consistent error handling it provides. This is also where sequence
# validation and uppercasing happens, so it's essential 
def parseFastaSeq(seqs:list[SeqRecord], pos:int, filePath:str) -> SeqRecord:
    """
    Retrieves and parses SeqRecord object at specified position in list, if available and
    valid DNA.

    Args:
        seqs (list[SeqRecord]): The list of SeqRecord objects.
        pos (int): The position at which to retrieve the desired sequence.
        filePath (str): File path of the original FASTA file the sequences were obtained from, for error reporting purposes.

    Raises:
        MissingSeqsErr: If the requested sequence position was invalid.
        InvalidSeqErr: if the sequence string is not valid DNA (only characters in the ACGTN set).

    Returns:
        SeqRecord: The sequence at the specified position.
    """
    try: seqObj = seqs[pos]
    except IndexError as err: raise MissingSeqErr(str(err),
        f"FASTA file \"{filePath}\" doesn't contain a sequence at position {pos}")
    
    if not isValidDNA(str(seqObj.seq)): raise InvalidSeqErr(
        f"sequence at position {pos} in FASTA file \"{filePath}\" is not valid DNA")

    # It makes no sense to use MutableSeq just to do this and never touch the seq again,
    # so I reassign instead:
    seqObj.seq = seqObj.seq.upper()
    return seqObj

# A method that checks which type of seq info was given and acts accordingly:
def parseSeq(seqOrFilePath:str, pos:int, name:SeqName) -> SeqRecord:
    """
    Obtains SeqRecord object from raw DNA sequence string or FASTA file path.

    Args:
        seqOrFilePath (str): The raw DNA sequence string or FASTA file path.
        pos (int): The position of the sequence in the FASTA file, if provided, as a non-negative integer.
        name (SeqName): Sequence metadata merely for output purposes, used if the sequence is provided as raw DNA.

    Raises:
        MissingSeqsErr: If the provided string was interpreted as a FASTA file but the requested sequence position was invalid.
        InvalidSeqErr: If the provided string could not be interpreted as a FASTA file but is also not valid DNA (only characters in the ACGTN set).

    Returns:
        SeqRecord: The parsed DNA sequence. Availability of certain metadata depends on how the sequence was retrieved.
    """
    if isValidFastaFilePath(seqOrFilePath):
        return parseFastaSeq(loadFasta(seqOrFilePath), pos, seqOrFilePath)
    
    try: return parseRawSeq(seqOrFilePath, name)
    except InvalidSeqErr as err: raise InvalidSeqErr(err.msg,
        f"your \"{name}\" sequence was interpreted as DNA, if you intended to provide a F\
ASTA file instead make sure your file path ends with a .fa or .fasta extension")

# Helper to isolate the "same file" case a bit more:
def parseSeqsFromFile(filePath:str, targetPos:int, queryPos:int) -> tuple[SeqRecord, SeqRecord]:
    """
    Obtain SeqRecord objects for both the target and query sequence, taken from the same
    FASTA file.

    Args:
        filePath (str): The provided FASTA file path.
        targetPos (int): The position of the target sequence in the FASTA file, as a non-negative integer.
        queryPos (int): The position of the query sequence in the FASTA file, as a non-negative integer.

    Raises:
        IdenticalSeqsErr: If the sequences to load are the same.

    Returns:
        tuple[SeqRecord, SeqRecord]: Respectively the target and query sequences.
    """
    # It makes no sense to align a sequence with itself:
    if queryPos == targetPos:
        if targetPos: raise IdenticalSeqsErr(
            f"tried loading sequences coming from the same file ({filePath})",
            f"with identical positions ({targetPos})")

        # If both positions are at 0 then they haven't been set (or the user feels very
        # funny), in this "one file" case the target and query seqs will be respectively
        # the first 2 seqs of the file:
        queryPos = 1
    
    seqs = loadFasta(filePath)
    return parseFastaSeq(seqs, targetPos, filePath), parseFastaSeq(seqs, queryPos, filePath)

def setupArgParser() -> ArgumentParser:
    """
    Setup an argparse.ArgumentParser instance and the input arguments needed for the
    analysis.

    Returns:
        Namespace: Object containing the arguments and their values as properties:
            - .target_seq (str): Target DNA sequence or FASTA file path.
            - .query_seq (str): Query DNA sequence or FASTA file path, if provided.
            - .target_pos (int): Target DNA sequence position in FASTA file.
            - .query_pos (int): Query DNA sequence position in FASTA file.
            - .match_score (int): Match score.
            - .mismatch_penalty (int): Mismatch penalty.
            - .gap_penalty (int): Linear gap penalty.
        
    All positions, scores and penalties are non-negative.
    """
    parser = ArgumentParser(prog = PACKAGE_NAME, description = PACKAGE_DESCR)

    # The user can specify a single seq arg, as long as it's a FASTA file path
    parser.add_argument("target_seq", type = str, help = TARGET_SEQ_HELP)
    parser.add_argument(
        "query_seq", type = str, nargs = '?', default = "", help = QUERY_SEQ_HELP)
    # ^^^ The use of _ here is necessary, as for some reason it doesn't convert
    # - to _ in pos args..

    # When a file containing the desired seq(s) is provided, the user can specify
    # which specific sequence of that file he wants. Negative positions are not
    # allowed.
    parser.add_argument(
        "--target-pos", "-tp", type = uint, help = TARGET_POS_HELP, default = 0)
    parser.add_argument(
        "--query-pos", "-qp", type = uint, help = QUERY_POS_HELP, default = 0)

    # As is customary for simple DNA alignment, all the relevant scores are
    # (non-negative) integers:
    parser.add_argument(
        "--match-score", "-m", type = uint, required = True, help = MATCH_HELP)
    parser.add_argument(
        "--mismatch-penalty", "-mm", type = uint, required = True, help = MISMATCH_HELP)
    parser.add_argument(
        "--gap-penalty", "-g", type = uint, required = True, help = GAP_HELP)

    return parser

def parseInputArgs(args:Namespace) -> tuple[SeqRecord, SeqRecord, int, int, int]:
    """
    Parse all the CLI input arguments passed by the user and necessary for
    the analysis.

    Args:
        args (Namespace): The Namespace object containing the arguments and their values.

    Raises:
        MissingSeqsErr: When not enough information was provided to retrieve the two sequences.
        IdenticalSeqsErr: If the sequences to load are the same.

    Returns:
        tuple:
        - SeqRecord: The target sequence.
        - SeqRecord: The query sequence.
        - int: The match score.
        - int: The mismatch penalty.
        - int: The gap penalty.
    """
    # It's impossible for this check to fail when query doesn't exist AND the 2 seqs are
    # the same, as target must exists:
    areTheSame = args.query_seq == args.target_seq
    if args.query_seq and not areTheSame:
        return (
            parseSeq(args.target_seq, args.target_pos, SeqName.Target),
            parseSeq(args.query_seq,  args.query_pos,  SeqName.Query),
            args.match_score, args.mismatch_penalty, args.gap_penalty)

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
        args.match_score, args.mismatch_penalty, args.gap_penalty)

# This function has been split into the 2 above for ease of testing, this down here is
# what one should export:
def collectAndParseInputArgs(args :tuple[str, ...]|None = None) -> tuple[SeqRecord, SeqRecord, int, int, int]:
    """
    Collect and parse all the CLI input arguments passed by the user and necessary for
    the analysis.
    
    Args:
        args (tuple[str, ...] | None): The input arguments, if passed manually for testing purposes. Defaults to: None.

    Returns:
        tuple:
        - SeqRecord: The target sequence.
        - SeqRecord: The query sequence.
        - int: The match score.
        - int: The mismatch penalty.
        - int: The gap penalty.
    """
    # Choosing to wait up to this point to call parse_args is convenient for testing purposes.
    return parseInputArgs(setupArgParser().parse_args(args))

if __name__ == "__main__":
    print(collectAndParseInputArgs(('0', '1', "-m", '0', "-mm", '0', "-g", '1')))