## A place to add package-level constants and constructs
from enum import StrEnum

# -Constants section-
MAX_DISPLAYED_SEQ_LEN        = 50
MAX_DISPLAYED_ALIGNMENTS     = 10
SCORE_MATRIX_SHMEM_NAME = "score_matrix"
DIRS_MATRIX_SHMEM_NAME  = "directions_matrix"

# -Strings section-
# Package description and documentation:
PACKAGE_NAME     = "ParaSeq"
PACKAGE_DESCR    = "A parallel Python implementation of the Smith-Waterman local sequence alignment algorithm."
TARGET_SEQ_HELP  = "Target DNA sequence or FASTA file path"
QUERY_SEQ_HELP   = "Query DNA sequence or FASTA file path"
TARGET_POS_HELP  = "Index (non-negative) of the target sequence in the FASTA file, if provided"
QUERY_POS_HELP   = "Index (non-negative) of the query sequence in the FASTA file, if provided"
MATCH_HELP       = "Match score, as a non-negative integer"
MISMATCH_HELP    = "Mismatch penalty, as a non-negative integer"
GAP_HELP         = "Constant gap penalty, as a non-negative integer"
OUT_PATH_HELP    = "Path to the output file the tool will create. Warning: will override if existing"
MAX_ALIGN_HELP   = "Maximum number of alignments shown in the terminal as output"
MAX_SEQ_LEN_HELP = "Maximum length for aligned chunks before sequences are truncated"

# Output:
ALIGNMENT_INFO = """
Target start pos: {}
Query start pos: {}
Target sequence: {}
Query sequence:  {}
"""

# Error messages:
UINT_ERR = "Expected a non-negative integer, got \"{}\"."
IDENTICAL_SEQS_PREFIX = "Alignment of identical sequences is pointless"
INVALID_SEQ_PREFIX    = "The provided sequence is not valid DNA as it contains characters outside of the ACGTN set"
MISSING_SEQ_PREFIX    = "Please provide at least 1 FASTA file path or 2 DNA sequences or FASTA file paths"
INVALID_FILE_PREFIX   = "The provided path or the corresponding file cannot be used"

# -Classes section-
class SeqName(StrEnum):
    """Enum type for the possible names of a DNA sequence."""
    Query  = "query-sequence"
    Target = "target-sequence"

    @property
    def id(self) -> str:
        """The sequence identifier, built from the 1st character of the name."""
        return self.value[0]