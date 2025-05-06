## A place to add package-level constants and constructs
from enum import StrEnum

# -Strings section-
# Package description and documentation:
PACKAGE_NAME  = "ParaSeq"
PACKAGE_DESCR = "A parallel Python implementation of the Smith-Waterman local\
sequence alignment algorithm."

TARGET_SEQ_HELP = "Target DNA sequence or FASTA file path"
QUERY_SEQ_HELP  = "Query DNA sequence or FASTA file path"

TARGET_POS_HELP = "Index (non-negative) of the target sequence in the FASTA file,\
if provided"
QUERY_POS_HELP  = "Index (non-negative) of the query sequence in the FASTA file,\
if provided"

MATCH_HELP    = "Match score, as a non-negative integer"
MISMATCH_HELP = "Mismatch penalty, as a non-negative integer"
GAP_HELP      = "Linear gap penalty, as a non-negative integer"

# -Classes section-
class SeqName(StrEnum):
    """Enum type for the possible names of a DNA sequence."""
    Query  = "query-sequence"
    Target = "target-sequence"

    @property
    def id(self) -> str:
        """The sequence identifier, built from the 1st character of the name."""
        return self.value[0]