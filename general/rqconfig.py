__author__ = 'letovesnoi'

import platform
import os
import sys

LOGGER_DEFAULT_NAME = 'rnaQUAST'

GENERAL_LOCATION = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

rnaOUAST_LOCATION = os.path.abspath(os.path.join(GENERAL_LOCATION, '..'))

BUSCO_LOCATION = os.path.abspath(os.path.join(rnaOUAST_LOCATION, 'others_libs', 'Busco'))

SUPPORTED_PYTHON_VERSIONS = ['2.5', '2.6', '2.7']

debug = False

# DEFINE MAGIC CONSTANTS:
# MULTITHREADING CONSTANTS:
DEFAULT_MAX_THREADS = 1

NUM_ONE_THREAD_TRANSCRIPTS = 100

# REPORT CONSTANTS:
PRECISION = 3

SORT_INDEX_LEN = 10000

# COVERAGE CONSTANTS:
class well_fully_coverage_thresholds():
    """thresholds for well/fully coverages"""

    def __init__(self, well_thresholds, fully_thresholds):
        self.well_exon_threshold = well_thresholds
        self.fully_exon_threshold = fully_thresholds

        self.well_isoform_threshold = well_thresholds
        self.fully_isoform_threshold = fully_thresholds

        self.well_block_threshold = well_thresholds
        self.fully_block_threshold = fully_thresholds

        self.well_transcript_threshold = well_thresholds
        self.fully_transcript_threshold = fully_thresholds


# written in the mapped coverages report isoforms:
COV_ISOFORMS_THRESHOLD = 0.5


TRANSCRIPT_LENS = [500, 1000]


GENEMARK_ST_GENE_LENS = [500, 1000]


# ALIGNMENTS UTILS CONSTANTS:
# ERR_CROSS_UNION must be less than MIN_SPLIT_ALIGNMENT_THRESHOLD:
class alignment_thresholds():
    """This threshold are used in best_alignment_set algorithms"""
    def __init__(self):
        self.ERR_CROSS_QUERY_UNION = 20

        self.UNION_PENALTY = 20

        self.ERR_CROSS_QUERY_FAKE_BLAT = 20
        self.ERR_SPACE_QUERY_FAKE_BLAT = 300

        self.ERR_CROSS_TARGET_FAKE_BLAT = 20
        self.ERR_SPACE_TARGET_FAKE_BLAT = 1500000 # 750000 blat max intron size

        self.QUERY_GAP_THRESHOLD = 10

        self.MIN_SPLIT_ALIGNMENT_THRESHOLD = 50

        self.LOW_COMPLEXITY_LEN_THRESHOLD = 10


# FUSION MISASSEMBLE CONSTANTS:
FUSION_ERR_THRESHOLD = 1


if platform.system() == 'Darwin':
    platform_name = 'macosx'
else:
    if platform.architecture()[0] == '64bit':
        platform_name = 'linux_64'
    else:
        platform_name = 'linux_32'