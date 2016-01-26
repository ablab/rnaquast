__author__ = 'letovesnoi'

from datetime import datetime

from quast23.libs import fastaparser

import TranscriptsCoverage


class AssemblyCorrectnessMetrics():
    """Class of annotation coverage metrics by aligned transcripts"""

    def __init__(self):
        # TRANSCRIPTS COVERAGE DICTIONARY:
        # self.transcripts_coverage_dict = {}

        # AVERAGE TRANSCRIPTS COVERAGE;
        self.transcripts_coverage = TranscriptsCoverage.TranscriptsCoverage()


    def update_assembly_correctness_metrics(self, aligned_transcript, transcript_coverage, well_fully_coverage_thresholds):
        start_time = datetime.now()
        # get aligned transcript coverage:
        # self.transcripts_coverage_dict[aligned_transcript.alignment.query_fragment.name] = transcript_coverage

        # update coverage of alignments:
        if transcript_coverage.id_mapped_isoform != None:
            self.transcripts_coverage.update_transcripts_coverage\
                (transcript_coverage, transcript_coverage.id_mapped_isoform,
                 aligned_transcript.alignment.query_fragment.name, well_fully_coverage_thresholds)

        elapsed_time = datetime.now() - start_time

        return elapsed_time


    def get_assembly_correctness_metrics(self, basic_metrics, simple_metrics, logger):
        self.transcripts_coverage.get_transcripts_coverage(simple_metrics, logger)


    def print_unannotated_transcripts(self, transcripts_dict, path_fa_unannotated, logger):
        logger.info('    Getting Unannotated transcripts report...')

        unannotated_transcripts = []
        for id_transcript in self.transcripts_coverage.ids_unannotated_transcripts:
            unannotated_transcripts.append((id_transcript, transcripts_dict[id_transcript]))

        fastaparser.write_fasta(path_fa_unannotated, unannotated_transcripts)

        logger.info('      saved to {}'.format(path_fa_unannotated))
