__author__ = 'letovesnoi'

from datetime import datetime

import TranscriptsCoverage


class AssemblyCorrectnessMetrics():
    """Class of annotation coverage metrics by aligned transcripts"""

    def __init__(self, args):
        # TRANSCRIPTS COVERAGE DICTIONARY:
        # self.transcripts_coverage_dict = {}

        # AVERAGE TRANSCRIPTS COVERAGE;
        self.transcripts_coverage = None
        if (args.gtf is not None or args.gene_db is not None) and \
                        args.alignment is not None and args.reference is not None and args.transcripts is not None:
            self.transcripts_coverage = TranscriptsCoverage.TranscriptsCoverage()


    def update_assembly_correctness_metrics(self, aligned_transcript, transcript_coverage, well_fully_coverage_thresholds):
        start_time = datetime.now()
        # get aligned transcript coverage:
        # self.transcripts_coverage_dict[aligned_transcript.alignment.query_fragment.name] = transcript_coverage

        # update coverage of alignments:
        if transcript_coverage.id_mapped_isoform is not None:
            self.transcripts_coverage.update_transcripts_coverage\
                (transcript_coverage, transcript_coverage.id_mapped_isoform,
                 aligned_transcript.alignment.query_fragment.name, well_fully_coverage_thresholds)

        elapsed_time = datetime.now() - start_time

        return elapsed_time


    def get_assembly_correctness_metrics(self, simple_metrics, logger):
        if self.transcripts_coverage is not None:
            self.transcripts_coverage.get_transcripts_coverage(simple_metrics, logger)