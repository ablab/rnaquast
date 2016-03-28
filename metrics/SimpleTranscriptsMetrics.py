__author__ = 'letovesnoi'

from datetime import datetime

from quast23.libs import N50
from quast23.libs import fastaparser


class SimpleTranscriptsMetrics():
    """Class of simple assembled transcripts metrics with alignment"""

    def __init__(self):
        # number of non misassembled aligned transcripts (best union is single union):
        self.num_non_misassembled = 0
        # ids of non misassembled aligned transcripts:
        self.ids_non_misassembled = set()

        # MISASSEMBLIES BY BLAT:
        # number of misassembled aligned transcripts (best union is multiple union):
        self.num_misassembled_by_blat = 0
        # ids of misassembled aligned transcripts:
        self.ids_misassembled_by_blat = set()

        # MISASSEMBLIES BY BLAT:
        # number of misassembled aligned transcripts (best union is multiple union):
        self.num_misassembled_by_blast = 0
        # ids of misassembled aligned transcripts:
        self.ids_misassembled_by_blast = set()

        # MISASSEMBLIES BY BLAT INTERSECTION MISASSEMBLIES BY BLASTN:
        self.ids_misassembled_together = set()

        # number of aligned transcripts:
        self.num_aligned = 0
        # ids of aligned transcripts:
        self.ids_aligned = set()
        # number of unaligned transcripts:
        self.num_unaligned = 0
        # percentage of unaligned transcripts:
        self.percent_unaligned = 0.0
        # ids of unaligned transcripts:
        self.ids_unaligned = set()

        # number of transcripts aligned to one place:
        self.num_unique_aligned = 0
        self.ids_unique_aligned = set()
        # number of transcripts aligned to more then one place:
        self.num_mul_aligned = 0
        # ids of multiple aligned transcripts:
        self.ids_mul_aligned = set()

        # MAYBE ONLY FOR TESTING:
        # ids and number of transcripts with query fragment aligned to different place in genome:
        # over all aligned transcripts:
        # self.num_mul_equal_query_aligned = 0
        # self.ids_mul_equal_query_aligned = set()
        # over misassembled transcripts:
        # self.num_mul_equal_query_aligned_mis = 0
        # self.ids_mul_equal_query_aligned_mis = set()

        # number of transcript alignments in PSL-file (not psl lines, but best mapped alignments):
        self.num_alignments = 0

        # ALL METRICS OVER ALIGNED NON MISASSEMBLED TRANSCRIPTS:
        # length of aligned transcript and number of alignments having this length over all alignments:
        self.len_distribution = {}
        # len:
        self.min_len = 0
        self.max_len = 0
        self.avg_len = 0.0
        self.tot_len = 0

        # reference (from begin to end of alignment) length of transcript  and number of transcripts having this length:
        self.ref_len_distribution = {}
        self.min_ref_len = 0
        self.max_ref_len = 0
        self.avg_ref_len = 0.0
        self.tot_ref_len = 0

        # alignment (summary length of aligned blocks) length of transcript and number of transcripts having this length:
        self.alignment_len_distribution = {}
        self.min_alignment_len = 0
        self.max_alignment_len = 0
        self.avg_alignment_len = 0.0
        self.tot_alignment_len = 0

        # fraction of transcript alignment length and number of transcripts having this fraction:
        self.fraction_distribution = {}
        # fraction len of transcripts:
        self.min_fraction = 0
        self.max_fraction = 0
        self.avg_fraction = 0.0
        self.tot_fraction = 0

        # number of mismatches in transcript and number of transcripts having this mismatches over all transcripts:
        self.mismatch_num_distribution = {}
        self.avg_mismatch_num = 0
        self.tot_mismatch_num = 0

        # target gap content:
        # number of target gap in transcript and number of transcripts having this gaps over all transcripts:
        self.tgap_num_distribution = {}
        self.min_tgap_num = 0
        self.max_tgap_num = 0
        self.avg_tgap_num = 0.0
        self.tot_tgap_num = 0
        # length of target gap and number of gaps having this length over all transcripts:
        self.tgap_len_distribution = {}
        self.min_tgap_len = 0
        self.max_tgap_len = 0
        self.avg_tgap_len = 0.0
        self.tot_tgap_len = 0

        # query gap content:
        # number of query gap in transcript and number of transcripts having this gaps over all transcripts:
        self.qgap_num_distribution = {}
        self.avg_qgap_num = 0.0
        self.tot_qgap_num = 0
        # length of query gap and number of gaps having this length over all transcripts:
        self.qgap_len_distribution = {}
        self.avg_qgap_len = 0.0
        self.tot_qgap_len = 0

        # blocks statistics:
        # blocks number distribution:
        self.blocks_num_distribution = {}
        # min number of block:
        self.min_blocks_num = 0
        # max number of block:
        self.max_blocks_num = 0
        # average number of blocks in aligned transcript:
        self.avg_blocks_num = 0.0
        # total number of blocks:
        self.tot_blocks_num = 0

        # blocks len distribution:
        self.blocks_len_distribution = {}
        # min length of block:
        self.min_block_len = 0
        # max length of block:
        self.max_block_len = 0
        # average length of block:
        self.avg_block_len = 0.0

        # NA50:
        self.na50 = 0
        self.alignment_lengths = []

        # ALL METRICS OVER ALIGNED TRANSCRIPTS:
        # number of alignments of transcript and count of transcripts having this count of alignments over all transcripts:
        self.mult_alignment_distribution = {}

        # FRACTION OF GENOME MAPPED:
        # dictionary of length of chromosomes:
        self.dict_len_chr = {}
        # summary length of all chromosomes:
        self.sum_len_chr = 0

        # percent of coverage of reference by all aligned transcripts:
        self.fraction_genome_mapped = 0.0

        # duplication ratio (average over all covered positions in chromosomes):
        # self.avg_duplication_ratio = 0.0
        # number of covered positions at least ones:
        # self.num_covered_pos = 0
        # number of aligned transcripts covering this position:
        # self.num_covering_pos_distribution = {}
        # for strand in strands:
        #     self.num_covering_pos_distribution[strand] = {}
        #     for name_seq in ids_chrs:
        #         self.num_covering_pos_distribution[strand][name_seq] = {}

        # INITIALIZE METRICS WITH ALIGNMENT AND ANNOTATION:
        # if args_annotation != None:
            # with alignment and annotation:
            # number of aligned transcripts len out isoforms without introns length:
            # self.num_outliers_len = 0
            # percent of transcripts having length out isoform range:
            # self.percent_outliers_len = 0.0
            # number of aligned transcripts having alignment length less then min isoform without introns length or more then max isoform without introns length:
            # self.num_outliers_alignment_len = 0
            # percent of aligned transcripts having alignment length out isoform range:
            # self.percent_outliers_alignment_len = 0.0
            # number and percentage of aligned transcripts with length more then max isoform with introns length or less then min isoform with introns length:
            # self.num_outliers_ref_len = 0
            # self.percentage_outliers_ref_len = 0.0
            # number and percentage of aligned transcripts with number of blocks more then max number of exons per isoform or less then min:
            # self.num_outliers_blocks_num = 0
            # self.percent_outliers_blocks_num = 0.0
            # number and percentage of blocks with length more then max length of exon and less then min length of exon:
            # self.num_outliers_blocks_len = 0
            # self.percent_outliers_blocks_len = 0.0
            # number and percentage of target gaps with length more then max length of intron and less then min length of intron:
            # self.num_outliers_tgap_len = 0
            # self.percent_outliers_tgap_len = 0
            # number and percentage of transcripts with number of target gaps more then max number of introns per isoform and less then min:
            # self.num_outliers_tgap_num = 0
            # self.percent_outliers_tgap_num = 0


    # update metrics with alignment, without annotation by best mapped aligned transcript
    # (can be only assembled transcript (single union alignment -- not misassembled):
    def update_metrics_by_best_mapped_transcript(self, aligned_transcript):
        start_time = datetime.now()

        self.ids_non_misassembled.add(aligned_transcript.alignment.query_fragment.name)

        self.tot_ref_len += aligned_transcript.reference_len
        if aligned_transcript.reference_len not in self.ref_len_distribution:
            self.ref_len_distribution[aligned_transcript.reference_len] = 0
        self.ref_len_distribution[aligned_transcript.reference_len] += 1

        self.tot_len += aligned_transcript.alignment.query_fragment.size
        if aligned_transcript.alignment.query_fragment.size not in self.len_distribution:
            self.len_distribution[aligned_transcript.alignment.query_fragment.size] = 0
        self.len_distribution[aligned_transcript.alignment.query_fragment.size] += 1

        self.tot_alignment_len += aligned_transcript.alignment_len
        if aligned_transcript.alignment_len not in self.alignment_len_distribution:
            self.alignment_len_distribution[aligned_transcript.alignment_len] = 0
        self.alignment_len_distribution[aligned_transcript.alignment_len] += 1

        self.tot_fraction += aligned_transcript.fraction
        if aligned_transcript.fraction not in self.fraction_distribution:
            self.fraction_distribution[aligned_transcript.fraction] = 0
        self.fraction_distribution[aligned_transcript.fraction] += 1

        self.tot_blocks_num += aligned_transcript.alignment.blocks_num
        if aligned_transcript.alignment.blocks_num not in self.blocks_num_distribution:
            self.blocks_num_distribution[aligned_transcript.alignment.blocks_num] = 0
        self.blocks_num_distribution[aligned_transcript.alignment.blocks_num] += 1

        self.alignment_lengths.append(aligned_transcript.alignment_len)

        for i_block in range(aligned_transcript.alignment.blocks_num):
            if aligned_transcript.alignment.blocks_sizes[i_block] not in self.blocks_len_distribution:
                self.blocks_len_distribution[aligned_transcript.alignment.blocks_sizes[i_block]] = 0
            self.blocks_len_distribution[aligned_transcript.alignment.blocks_sizes[i_block]] += 1

            if i_block + 1 < aligned_transcript.alignment.blocks_num:
                # length of target gap and number of gaps having this length over all transcripts:
                t_base_insert = aligned_transcript.alignment.target_fragment.starts[i_block + 1] - \
                                (aligned_transcript.alignment.target_fragment.starts[i_block] +
                                 aligned_transcript.alignment.blocks_sizes[i_block] - 1) - 1
                if t_base_insert not in self.tgap_len_distribution:
                    self.tgap_len_distribution[t_base_insert] = 0
                self.tgap_len_distribution[t_base_insert] += 1

                # length of query gap and number of gaps having this length over all transcripts:
                q_base_insert = aligned_transcript.alignment.query_fragment.starts[i_block + 1] - \
                                (aligned_transcript.alignment.query_fragment.starts[i_block] +
                                 aligned_transcript.alignment.blocks_sizes[i_block] - 1) - 1
                if q_base_insert not in self.qgap_len_distribution:
                    self.qgap_len_distribution[q_base_insert] = 0
                self.qgap_len_distribution[q_base_insert] += 1

        # mismatches:
        if aligned_transcript.alignment.mismatches not in self.mismatch_num_distribution:
            self.mismatch_num_distribution[aligned_transcript.alignment.mismatches] = 0
        self.mismatch_num_distribution[aligned_transcript.alignment.mismatches] += 1
        self.tot_mismatch_num += aligned_transcript.alignment.mismatches

        # target gap content:
        # number of target gap in transcript and number of transcripts having this gaps over all transcripts:
        self.tot_tgap_num += aligned_transcript.alignment.target_fragment.num_insert

        if aligned_transcript.alignment.target_fragment.num_insert not in self.tgap_num_distribution:
            self.tgap_num_distribution[aligned_transcript.alignment.target_fragment.num_insert] = 0
        self.tgap_num_distribution[aligned_transcript.alignment.target_fragment.num_insert] += 1

        self.tot_tgap_len += aligned_transcript.alignment.target_fragment.base_insert
        self.tot_qgap_len += aligned_transcript.alignment.query_fragment.base_insert

        # query gap content:
        # number of query gap in transcript and number of transcripts having this gaps over all transcripts:
        self.tot_qgap_num += aligned_transcript.alignment.query_fragment.num_insert
        if aligned_transcript.alignment.query_fragment.num_insert not in self.qgap_num_distribution:
            self.qgap_num_distribution[aligned_transcript.alignment.query_fragment.num_insert] = 0
        self.qgap_num_distribution[aligned_transcript.alignment.query_fragment.num_insert] += 1

        # get duplication ratio over chromosomes:
        # for i_block in range(aligned_transcript.alignment.blocks_num):
        #     for i_pos in range(aligned_transcript.alignment.target_fragment.starts[i_block], aligned_transcript.alignment.target_fragment.ends[i_block] + 1):
        #         if i_pos not in self.num_covering_pos_distribution[aligned_transcript.strand][aligned_transcript.alignment.target_fragment.name]:
        #             self.num_covered_pos += 1
        #             self.num_covering_pos_distribution[aligned_transcript.strand][aligned_transcript.alignment.target_fragment.name][i_pos] = 0
        #         self.num_covering_pos_distribution[aligned_transcript.strand][aligned_transcript.alignment.target_fragment.name][i_pos] += 1

        # if basic_isoforms_metrics != None:
        #     self.update_metrics_by_best_aligned_transcript_w_annotation(aligned_transcript, basic_isoforms_metrics)

        elapsed_time = datetime.now() - start_time
        return elapsed_time

    # update metrics with alignment and annotation by best mapped aligned transcript
    # (can be only assembled transcript (single union alignment -- non misassembled):
    # def update_metrics_by_best_aligned_transcript_w_annotation(self, aligned_transcript, basic_isoforms_metrics):
    #     if aligned_transcript.reference_len > basic_isoforms_metrics.max_len_w_introns or \
    #                     aligned_transcript.reference_len < basic_isoforms_metrics.min_len_w_introns:
    #         self.num_outliers_ref_len += 1
    #     if aligned_transcript.alignment.query_fragment.size > basic_isoforms_metrics.max_len_wout_introns or \
    #                     aligned_transcript.alignment.query_fragment.size < basic_isoforms_metrics.min_len_wout_introns:
    #         self.num_outliers_len += 1
    #     if aligned_transcript.alignment_len > basic_isoforms_metrics.max_len_wout_introns or \
    #                     aligned_transcript.alignment_len < basic_isoforms_metrics.min_len_wout_introns:
    #         self.num_outliers_alignment_len += 1
    #     if aligned_transcript.alignment.blocks_num > basic_isoforms_metrics.max_exons_num or \
    #                     aligned_transcript.alignment.blocks_num < basic_isoforms_metrics.min_exons_num:
    #         self.num_outliers_blocks_num += 1
    #     if aligned_transcript.alignment.target_fragment.num_insert > basic_isoforms_metrics.max_introns_num or \
    #                     aligned_transcript.alignment.target_fragment.num_insert < basic_isoforms_metrics.min_introns_num:
    #         self.num_outliers_tgap_num += 1


    # update metrics by best mapped alignments
    # (can be only assembled transcript alignments (single union alignment -- not misassembled):
    def update_metrics_by_best_mapped_alignments(self, best_mapped_alignments):
        start_time = datetime.now()

        self.num_alignments += len(best_mapped_alignments)
        if len(best_mapped_alignments) > 0:
            self.ids_aligned.add(best_mapped_alignments[0].query_fragment.name)
        if len(best_mapped_alignments) == 1:
            self.num_unique_aligned += 1
            self.ids_unique_aligned.add(best_mapped_alignments[0].query_fragment.name)
        if len(best_mapped_alignments) > 1:
            self.num_mul_aligned += 1
            self.ids_mul_aligned.add(best_mapped_alignments[0].query_fragment.name)
            if len(best_mapped_alignments) not in self.mult_alignment_distribution:
                self.mult_alignment_distribution[len(best_mapped_alignments)] = 0
            self.mult_alignment_distribution[len(best_mapped_alignments)] += 1

        elapsed_time = datetime.now() - start_time
        return elapsed_time


    # update metrics by best alignments for misassembled transcript:
    def update_metrics_by_misassembled_alignments(self, best_alignments, is_psl_lines):
        self.ids_aligned.add(best_alignments[0].query_fragment.name)
        if is_psl_lines == True:
            self.ids_misassembled_by_blat.add(best_alignments[0].query_fragment.name)
        else:
            self.ids_misassembled_by_blast.add(best_alignments[0].query_fragment.name)

    # update metrics by all alignments
    # def update_metrics_by_all_alignments(self, all_alignments):
    #     TEMPORARY FOR TEST BLAT ALIGNER:
    #     self.update_num_mul_equal_query_aligned(all_alignments)
    #     self.update_num_mul_equal_query_aligned(all_alignments)


    def get_metrics(self, args_blast, reference_dict, transcripts_dict, num_assembled, logger):
        logger.info('  Getting ALIGNMENT metrics...')

        self.num_non_misassembled = len(self.ids_non_misassembled)

        self.num_aligned = len(self.ids_aligned)

        self.ids_unaligned = set(transcripts_dict.keys()) - self.ids_aligned
        self.num_unaligned = num_assembled - self.num_aligned

        # self.ids_misassembled = self.ids_aligned - self.ids_non_misassembled
        self.num_misassembled_by_blat = len(self.ids_misassembled_by_blat)
        self.num_misassembled_by_blast = len(self.ids_misassembled_by_blast)

        if args_blast:
            self.ids_misassembled_together = self.ids_misassembled_by_blat.intersection(self.ids_misassembled_by_blast)
        else:
            self.ids_misassembled_together = self.ids_misassembled_by_blat
        self.num_misassembled_together = len(self.ids_misassembled_together)

        # TEMPORARY FOR TEST BLAT ALIGNER:
        # self.num_mul_equal_query_aligned = len(self.ids_mul_equal_query_aligned)
        # self.ids_mul_equal_query_aligned_mis = self.ids_mul_equal_query_aligned.intersection(self.ids_misassembled)
        # self.num_mul_equal_query_aligned_mis = len(self.ids_mul_equal_query_aligned_mis)

        if num_assembled != 0:
            self.percent_unaligned = self.num_unaligned / num_assembled

        # get averages:
        if self.num_alignments != 0:
            # alignment length:
            self.avg_ref_len = self.tot_ref_len / self.num_alignments
            self.max_ref_len = max(self.ref_len_distribution.keys())
            self.min_ref_len = min(self.ref_len_distribution.keys())

            # alignment length:
            self.avg_alignment_len = self.tot_alignment_len * 1.0 / self.num_alignments
            self.max_alignment_len = max(self.alignment_len_distribution.keys())
            self.min_alignment_len = min(self.alignment_len_distribution.keys())
            # fraction:
            self.avg_fraction = self.tot_fraction * 1.0 / self.num_alignments
            self.max_fraction = max(self.fraction_distribution.keys())
            self.min_fraction = min(self.fraction_distribution.keys())

            # blocks:
            # number:
            self.avg_blocks_num = self.tot_blocks_num * 1.0 / self.num_alignments
            self.max_blocks_num = max(self.blocks_num_distribution.keys())
            self.min_blocks_num = min(self.blocks_num_distribution.keys())
            # len:
            self.avg_block_len = self.tot_alignment_len * 1.0 / self.tot_blocks_num
            self.max_block_len = max(self.blocks_len_distribution.keys())
            self.min_block_len = min(self.blocks_len_distribution.keys())

            self.avg_tgap_num = self.tot_tgap_num * 1.0 / self.num_alignments
            self.min_tgap_num = min(self.tgap_num_distribution.keys())
            self.max_tgap_num = max(self.tgap_num_distribution.keys())

            if self.tot_tgap_num != 0:
                self.min_tgap_len = min(self.tgap_len_distribution.keys())
                self.max_tgap_len = max(self.tgap_len_distribution.keys())
                self.avg_tgap_len = self.tot_tgap_len * 1.0 / self.tot_tgap_num

            self.avg_qgap_num = self.tot_qgap_num * 1.0 / self.num_alignments

            if self.tot_qgap_num != 0:
                self.avg_qgap_len = self.tot_qgap_len * 1.0 / self.tot_qgap_num

            self.avg_mismatch_num = self.tot_mismatch_num * 1.0 / self.num_alignments

            self.avg_len = self.tot_len * 1.0 / self.num_alignments
            self.max_len = max(self.len_distribution.keys())
            self.min_len = min(self.len_distribution.keys())

        # percent of fraction of genome:
        for id_chr in reference_dict:
            self.dict_len_chr[id_chr] = len(reference_dict[id_chr])
            self.sum_len_chr += self.dict_len_chr[id_chr]

        # self.fraction_genome_mapped = self.num_covered_pos * 1.0 / (len(strands) * self.sum_len_chr)

        # for strand in self.num_covering_pos_distribution:
        #     for seq_name in self.num_covering_pos_distribution[strand]:
        #         self.avg_duplication_ratio += sum(self.num_covering_pos_distribution[strand][seq_name].values())
        # if self.num_covered_pos != 0:
        #     self.avg_duplication_ratio /= self.num_covered_pos

        self.na50 = N50.N50(self.alignment_lengths)

        # if args_annotation != None:
        #     self.get_metrics_w_annotation(basic_isoforms_metrics)

        logger.info('  Done.')


    # def get_metrics_w_annotation(self, basic_isoforms_metrics):
    #     if self.num_alignments != 0:
    #         self.percentage_outliers_ref_len = self.num_outliers_ref_len * 1.0 / self.num_alignments
    #         self.percent_outliers_blocks_num = self.num_outliers_blocks_num * 1.0 / self.num_alignments
    #
    #         self.percent_outliers_len = self.num_outliers_len * 1.0 / self.num_alignments
    #         self.percent_outliers_alignment_len = self.num_outliers_alignment_len * 1.0 / self.num_alignments
    #
    #         for blocks_len in self.blocks_len_distribution:
    #             if blocks_len < basic_isoforms_metrics.min_exon_len or blocks_len > basic_isoforms_metrics.max_exon_len:
    #                 self.num_outliers_blocks_len += self.blocks_len_distribution[blocks_len]
    #         self.percent_outliers_blocks_len = self.num_outliers_blocks_len * 1.0 / self.tot_blocks_num
    #
    #     if self.tot_tgap_num != 0:
    #         self.percent_outliers_tgap_num = self.num_outliers_tgap_num * 1.0 / self.tot_tgap_num
    #
    #         for tgap_len in self.tgap_len_distribution:
    #             if tgap_len < basic_isoforms_metrics.min_intron_len or tgap_len > basic_isoforms_metrics.max_intron_len:
    #                 self.num_outliers_tgap_len += self.tgap_len_distribution[tgap_len]
    #             self.percent_outliers_tgap_len = self.num_outliers_tgap_len * 1.0 / self.tot_tgap_num


    # MAYBE ONLY FOR TESTING BLAT ALIGNER:
    # def update_num_exact_mul_equal_query_aligned(self, alignments):
    #     positions = []
    #     for alignment in alignments:
    #         if (alignment.query_fragment.start, alignment.query_fragment.end) not in positions:
    #             positions.append((alignment.query_fragment.start, alignment.query_fragment.end))
    #         else:
    #             self.ids_mul_equal_query_aligned.add(alignment.query_fragment.name)
    #             return


    # MAYBE ONLY FOR TESTING BLAT ALIGNER:
    # def update_num_mul_equal_query_aligned(self, alignments):
    #     for i_alignment in range(len(alignments)):
    #         for j_alignment in range(i_alignment + 1, len(alignments), 1):
    #             start = max(alignments[i_alignment].query_fragment.start, alignments[j_alignment].query_fragment.start)
    #             end = min(alignments[i_alignment].query_fragment.end, alignments[j_alignment].query_fragment.end)
    #             cross = max(0, end - start)
    #             length = max(alignments[i_alignment].query_fragment.end - alignments[i_alignment].query_fragment.start,
    #                          alignments[j_alignment].query_fragment.end - alignments[j_alignment].query_fragment.start)
    #             if cross > 0.5 * length:
    #                 self.ids_mul_equal_query_aligned.add(alignments[0].query_fragment.name)
    #                 return



    def print_fa_transcripts(self, transcripts_dict, ids_transcripts, path_fasta, logger, transcripts_name=''):
        logger.info('    Getting {} transcripts report...'.format(transcripts_name))

        chosen_transcripts = []
        for id_transcript in ids_transcripts:
            chosen_transcripts.append((id_transcript, transcripts_dict[id_transcript]))

        fastaparser.write_fasta(path_fasta, chosen_transcripts)

        logger.info('      saved to {}'.format(path_fasta))
