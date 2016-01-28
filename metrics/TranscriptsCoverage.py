__author__ = 'lenk'

from quast23.libs import fastaparser


class TranscriptsCoverage(object):
    """ASSEMBLY CORRECTNESS METRICS: coverage of aligned transcripts over all scaffolds/chromosomes/patches"""

    def __init__(self):
        # DISTRIBUTIONS FOR TRANSCRIPTS ALIGNMENTS:
        # average percents of covered bases over all aligned transcripts and distributions:
        self.avg_covered_fraction_aligned_part = 0.0
        self.covered_fraction_aligned_part_distribution = {}

        self.avg_covered_fraction_whole_transcript = 0.0
        self.covered_fraction_whole_transcript_distribution = {}

        self.avg_covered_fraction_block_in_t = 0.0
        self.avg_covered_fraction_block_in_t_distribution = {}

        self.avg_percentage_well_covered_blocks_in_t = 0.0
        self.percentage_well_covered_blocks_in_t_distribution = {}

        self.avg_percentage_fully_covered_blocks_in_t = 0.0
        self.percentage_fully_covered_blocks_in_t_distribution = {}

        # FOR ALL TRANSCRIPTS ALIGNMENTS:
        self.matched_len = 0
        self.unmatched_len = 0

        self.ids_annotated_transcripts = set()
        self.num_annotated_transcripts = 0

        self.ids_unannotated_transcripts = set()
        self.num_unannotated_transcripts = 0
        self.percentage_unannotated_transcripts = 0.0

        # ids, number and percentage of of transcripts well-covered (30% over aligned part) by specific annotated isoform:
        self.ids_well_covered_transcripts = set()
        self.num_well_covered_transcripts = 0
        self.percentage_well_covered_transcripts = 0.0
        # ids, number and percentage of transcripts fully-covered (90% over aligned part) by specific annotated isoform:
        self.ids_fully_covered_transcripts = set()
        self.num_fully_covered_transcripts = 0
        self.percentage_fully_covered_transcripts = 0.0

        # summary count of covered blocks by annotated exons over all transcripts:
        self.num_well_covered_blocks = 0
        self.num_fully_covered_blocks = 0
        # percentages of well/fully-covered blocks over all aligned transcripts:
        self.percentage_well_covered_blocks = 0.0
        self.percentage_fully_covered_blocks = 0.0

        # average percentage of matched bases of block over all aligned blocks:
        self.avg_covered_fraction_block = 0.0
        # number of all blocks from all transcript alignments that have x% matched to an database isoform:
        self.block_fraction_matched_distribution = {}

        # summary number of blocks in aligned transcripts:
        self.tot_blocks_num = 0


    # update coverage of aligned transcripts by coverage of aligned blocks in one aligned transcript by specific isoform:
    def update_transcripts_coverage(self, transcript_coverage, id_specific_isoform, id_transcript, well_fully_coverage_thresholds):
        # DISTRIBUTIONS FOR TRANSCRIPTS ALIGNMENTS:
        self.avg_covered_fraction_aligned_part += transcript_coverage.covered_fraction_aligned_part[id_specific_isoform]
        if transcript_coverage.covered_fraction_aligned_part[id_specific_isoform] not in self.covered_fraction_aligned_part_distribution:
            self.covered_fraction_aligned_part_distribution[transcript_coverage.covered_fraction_aligned_part[id_specific_isoform]] = 0
        self.covered_fraction_aligned_part_distribution[transcript_coverage.covered_fraction_aligned_part[id_specific_isoform]] += 1

        self.avg_covered_fraction_whole_transcript += transcript_coverage.covered_fraction_whole_transcript[id_specific_isoform]
        if transcript_coverage.covered_fraction_whole_transcript[id_specific_isoform] not in self.covered_fraction_whole_transcript_distribution:
            self.covered_fraction_whole_transcript_distribution[transcript_coverage.covered_fraction_whole_transcript[id_specific_isoform]] = 0
        self.covered_fraction_whole_transcript_distribution[transcript_coverage.covered_fraction_whole_transcript[id_specific_isoform]] += 1

        self.avg_covered_fraction_block_in_t += transcript_coverage.avg_covered_fraction_block[id_specific_isoform]
        if transcript_coverage.avg_covered_fraction_block[id_specific_isoform] not in self.avg_covered_fraction_block_in_t_distribution:
            self.avg_covered_fraction_block_in_t_distribution[transcript_coverage.avg_covered_fraction_block[id_specific_isoform]] = 0
        self.avg_covered_fraction_block_in_t_distribution[transcript_coverage.avg_covered_fraction_block[id_specific_isoform]] += 1

        self.avg_percentage_well_covered_blocks_in_t += transcript_coverage.percentage_well_covered_blocks[id_specific_isoform]
        if transcript_coverage.percentage_well_covered_blocks[id_specific_isoform] not in self.percentage_well_covered_blocks_in_t_distribution:
            self.percentage_well_covered_blocks_in_t_distribution[transcript_coverage.percentage_well_covered_blocks[id_specific_isoform]] = 0
        self.percentage_well_covered_blocks_in_t_distribution[transcript_coverage.percentage_well_covered_blocks[id_specific_isoform]] += 1

        self.avg_percentage_fully_covered_blocks_in_t += transcript_coverage.percentage_fully_covered_blocks[id_specific_isoform]
        if transcript_coverage.percentage_fully_covered_blocks[id_specific_isoform] not in self.percentage_fully_covered_blocks_in_t_distribution:
            self.percentage_fully_covered_blocks_in_t_distribution[transcript_coverage.percentage_fully_covered_blocks[id_specific_isoform]] = 0
        self.percentage_fully_covered_blocks_in_t_distribution[transcript_coverage.percentage_fully_covered_blocks[id_specific_isoform]] += 1

        # FOR ALL TRANSCRIPTS ALIGNMENTS:
        for i_block in range(transcript_coverage.blocks_num):
            block_fraction_matched = transcript_coverage.covered_fraction_blocks[id_specific_isoform][i_block]
            if block_fraction_matched not in self.block_fraction_matched_distribution:
                self.block_fraction_matched_distribution[block_fraction_matched] = 0
            self.block_fraction_matched_distribution[block_fraction_matched] += 1

        self.avg_covered_fraction_block += sum(transcript_coverage.covered_fraction_blocks[id_specific_isoform])
        self.tot_blocks_num += transcript_coverage.blocks_num

        if transcript_coverage.covered_bases[id_specific_isoform] >= well_fully_coverage_thresholds.well_transcript_threshold * transcript_coverage.transcript_len:
            self.ids_well_covered_transcripts.add(id_transcript)
        if transcript_coverage.covered_bases[id_specific_isoform] >= well_fully_coverage_thresholds.fully_transcript_threshold * transcript_coverage.transcript_len:
            self.ids_fully_covered_transcripts.add(id_transcript)

        self.ids_annotated_transcripts.add(id_transcript)

        self.num_well_covered_blocks += transcript_coverage.num_well_covered_blocks[id_specific_isoform]
        self.num_fully_covered_blocks += transcript_coverage.num_fully_covered_blocks[id_specific_isoform]

        self.matched_len += transcript_coverage.covered_bases[id_specific_isoform]


    def get_transcripts_coverage(self, simple_transcripts_metrics, logger):
        logger.info('  Getting SPECIFICITY metrics...')

        self.num_annotated_transcripts = len(self.ids_annotated_transcripts)
        self.num_well_covered_transcripts = len(self.ids_well_covered_transcripts)
        self.num_fully_covered_transcripts = len(self.ids_fully_covered_transcripts)

        self.unmatched_len = simple_transcripts_metrics.tot_len - self.matched_len

        # DISTRIBUTIONS FOR TRANSCRIPTS ALIGNMENT:
        if simple_transcripts_metrics.num_alignments != 0:
            self.avg_covered_fraction_aligned_part /= simple_transcripts_metrics.num_alignments
            self.avg_covered_fraction_whole_transcript /= simple_transcripts_metrics.num_alignments
            self.avg_covered_fraction_block_in_t /= simple_transcripts_metrics.num_alignments
            self.avg_percentage_well_covered_blocks_in_t /= simple_transcripts_metrics.num_alignments
            self.avg_percentage_fully_covered_blocks_in_t /= simple_transcripts_metrics.num_alignments

        # FOR ALL TRANSCRIPTS:
        if simple_transcripts_metrics.num_non_misassembled != 0:
            self.ids_unannotated_transcripts = simple_transcripts_metrics.ids_non_misassembled.difference(self.ids_annotated_transcripts)
            self.num_unannotated_transcripts = simple_transcripts_metrics.num_non_misassembled - self.num_annotated_transcripts
            self.percentage_unannotated_transcripts = self.num_unannotated_transcripts * 1.0 / simple_transcripts_metrics.num_non_misassembled
            self.percentage_well_covered_transcripts = self.num_well_covered_transcripts * 1.0 / simple_transcripts_metrics.num_non_misassembled
            self.percentage_fully_covered_transcripts = self.num_fully_covered_transcripts * 1.0 / simple_transcripts_metrics.num_non_misassembled

        # FOR ALL TRANSCRIPTS ALIGNMENTS:
        if self.tot_blocks_num != 0:
            self.percentage_well_covered_blocks = self.num_well_covered_blocks * 1.0 / self.tot_blocks_num
            self.percentage_fully_covered_blocks = self.num_fully_covered_blocks * 1.0 / self.tot_blocks_num
            self.avg_covered_fraction_block /= self.tot_blocks_num

        logger.info('  Done.')


    def print_unannotated_transcripts(self, transcripts_dict, path_fa_unannotated, logger):
        logger.info('    Getting Unannotated transcripts report...')

        unannotated_transcripts = []
        for id_transcript in self.ids_unannotated_transcripts:
            unannotated_transcripts.append((id_transcript, transcripts_dict[id_transcript]))

        fastaparser.write_fasta(path_fa_unannotated, unannotated_transcripts)

        logger.info('      saved to {}'.format(path_fa_unannotated))
