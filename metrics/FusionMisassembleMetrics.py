__author__ = 'letovesnoi'

import os
import sys
import subprocess

from general import UtilsGeneral
from general import UtilsPipeline

class FusionMisassembleMetrics():
    """Class of metrics for fusion genes or misassembled transcripts"""

    def __init__(self, args_fusion_misassemble_analyze, transcripts_metrics, transcripts_file, sam_file, logger):
        # number of transcripts with non-overlapping (here and after implies slightly overlapping) blocks mapping to isoforms of different genes within a single chromosome.
        # transcripts containing consecutive genes sequences -- neither fusions nor misassemblies.
        self.mapped_one_chr_consecutive_genes_num = 0
        #  transcripts containing non consecutive genes sequences -- fusions or misassemblies.
        self.mapped_one_chr_non_consecutive_genes_num = 0
        # number of transcripts with non-overlapping alignments (from different chromosomes) covering isoforms of different genes -- fusions or misassemblies.
        self.mapped_several_chr_genes_num = 0
        # number of transcripts with non-overlapping blocks covering different isoforms of the same gene (but not mapping fully to any isoform) are candidates for new isoforms.
        self.mapped_several_isoforms_num = 0
        # number of transcripts with non-overlapping blocks/alignments covering both annotated and unannotated regions -- fusion or misassembly.
        self.mapped_ann_unann_num = 0
        # number of transcripts with strongly overlapping blocks/alignments can represent misassemblies (caused by repeats) or map to paralogous  genes
        self.mapped_with_overlap_num = 0

        # number of suspected misassemblies or if --afm option: suspected misassemblies confirmed by paired reads:
        self.misassemble_num = 0
        # number of suspected fusions or if --afm option: suspected fusions confirmed by paired reads:
        self.fusion_num = 0

        self.mis_by_reads_file = None
        self.misassemble_by_reads_dict = None
        #if args_fusion_misassemble_analyze:
        #    self.mis_by_reads_file = UtilsPipeline.get_file_with_misassemblies_by_reads(args, logger, transcripts_metrics, transcripts_file, sam_file)

        #    logger.info('  Setting misassemblies by reads for transcripts in {}:'.format(transcripts_file))
        #    self.misassemble_by_reads_dict = self.parse_out_finding_misassemblies_by_reads(self.mis_by_reads_file)
        #    logger.info('  Done.')


    def get_metrics(self):
        pass


    def is_suspected_fusion(self, best_union_alignments, fusion_err_threshhold):
        strands = ['+', '-']
        q_starts = []
        q_ends = []
        t_starts = {'+': {}, '-': {}}
        t_ends = {'+': {}, '-': {}}
        t_starts_sort_index = {'+': {}, '-': {}}
        t_starts_sort_array = {'+': {}, '-': {}}
        t_ends_sort_index = {'+': {}, '-': {}}
        t_ends_sort_array = {'+': {}, '-': {}}

        for i in range(len(best_union_alignments)):
            q_starts.append(best_union_alignments[i].query_fragment.start)
            q_ends.append(best_union_alignments[i].query_fragment.end)
            strand = best_union_alignments[i].strand
            id_chr = best_union_alignments[i].target_fragment.name
            t_start = best_union_alignments[i].target_fragment.start
            t_end = best_union_alignments[i].target_fragment.end
            if id_chr not in t_starts:
                t_starts[strand][id_chr] = []
                t_ends[strand][id_chr] = []
            t_starts[strand][id_chr].append(t_start)
            t_ends[strand][id_chr].append(t_end)

        q_sort_index, q_sort_array = UtilsGeneral.get_order_indexes_elements(q_starts)
        for i in range(len(q_sort_index) - 1):
            i_alignment0 = q_sort_index[i]
            i_alignment1 = q_sort_index[i + 1]
            start = min(q_ends[i_alignment0], q_ends[i_alignment1])
            end = max(q_starts[i_alignment0], q_starts[i_alignment1])
            if abs(end - start) + 1 > fusion_err_threshhold:
                return 0

        for strand in strands:
            for id_chr in t_starts[strand]:
                t_starts_sort_index[strand][id_chr], t_starts_sort_array[strand][id_chr] = UtilsGeneral.get_order_indexes_elements(t_starts[strand][id_chr])
                t_ends_sort_index[strand][id_chr], t_ends_sort_array[strand][id_chr] = UtilsGeneral.get_order_indexes_elements(t_ends[strand][id_chr])
                if t_starts_sort_index[strand][id_chr] != t_ends_sort_index[strand][id_chr]:
                    return 0

        return 1

    @classmethod
    def parse_out_finding_misassemblies_by_reads(cls, mis_by_reads_file):
        misassemble_by_reads_dict = {}
        query = ''
        with open(mis_by_reads_file, 'r') as fin:
            line = fin.readline().strip().split('\t')
            while line[0] != '':
                if 'Fragments' in line[0]:
                    line = fin.readline().strip().split('\t')
                    name = line
                    line = fin.readline().strip().split('\t')
                if line[0] == 'query':
                    query = line
                    if query[1] not in misassemble_by_reads_dict:
                        misassemble_by_reads_dict[query[1]] = []
                    line = fin.readline().strip().split('\t')
                if 'Fragments' not in line[0] and 'query' not in line[0] and 'start' not in line and line[0] != '':
                    pos = line
                    misassemble_by_reads_dict[query[1]].append((int(pos[0]), int(pos[1])))
                line = fin.readline().strip().split('\t')
        return misassemble_by_reads_dict

    def get_suspected_fusion_misassemble(self, best_union_alignments, fusion_err_threshhold):
        # getting suspected fusions and misassemblies
        is_suspected_misassemble = 0
        is_suspected_fusion = 0

        is_suspected_fusion = self.is_suspected_fusion(best_union_alignments, fusion_err_threshhold)
        self.fusion_num += is_suspected_fusion
        if is_suspected_fusion == 0:
            is_suspected_misassemble = 1
            self.misassemble_num += is_suspected_misassemble

        return is_suspected_fusion, is_suspected_misassemble

    def get_confirmed_fusion_misassemblies(self, best_union_alignments, fusion_err_threshhold):
        # getting confirmed fusions and misassemblies
        q_name = best_union_alignments[0].query_fragment.name
        is_confirmed_fus = 0
        is_confirmed_mis = 0

        is_suspected_fusion, is_suspected_misassemble = self.get_suspected_fusion_misassemble(best_union_alignments, fusion_err_threshhold)
        if is_suspected_misassemble == 1:
            is_confirmed_mis = 1
        elif q_name not in self.misassemble_by_reads_dict and is_suspected_fusion == 1:
            is_confirmed_fus = 1
        elif is_suspected_fusion == 1:
            # say that misassembled if ends of neighboring aligned parts of transcript are into misassembled by reads intervals:
            q_starts_alignments = []
            for i_alignment in range(len(best_union_alignments)):
                q_starts_alignments.append(best_union_alignments[i_alignment].query_fragment.start)
            i_q_sort_starts_alignments, q_sort_starts_alignments = UtilsGeneral.get_order_indexes_elements(q_starts_alignments)
            for i in range(len(i_q_sort_starts_alignments) - 1):
                alignment0 = best_union_alignments[i_q_sort_starts_alignments[i]]
                alignment1 = best_union_alignments[i_q_sort_starts_alignments[i + 1]]
                q_end_alignment = alignment0.query_fragment.end
                q_start_alignment = alignment1.query_fragment.start
                q_start = min(q_end_alignment, q_start_alignment)
                q_end = max(q_end_alignment, q_start_alignment)
                for interval in self.misassemble_by_reads_dict[q_name]:
                    mis_start = interval[0]
                    mis_end = interval[1]
                    if (q_start > mis_start and q_start < mis_end) or (q_end > mis_start and q_end < mis_end):
                        is_confirmed_mis = 1
                        self.fusion_num -= 1
                        self.misassemble_num += 1
                    if is_confirmed_mis == 1:
                        break
                if is_confirmed_mis == 1:
                    break
            if is_confirmed_mis == 0:
                is_confirmed_fus = 1
        #logger.debug('      Done.')
        return is_confirmed_mis, is_confirmed_fus

    def update_metrics(self, args, best_union_lines, best_union_alignments,  single_transcript_lines, transcripts_metrics_report, fusion_err_threshhold):
        if len(best_union_alignments) > 1:
            # updating metrics by best union alignments for misassembled transcript
            is_fusion = 0
            is_misassemble = 0

            # COUNT SUSPECTED FUSIONS AND MISASSEMBLIES:
            #if not args.fusion_misassemble_analyze:
            is_fusion, is_misassemble = \
                self.get_suspected_fusion_misassemble(best_union_alignments, fusion_err_threshhold)

            # COUNT CONFIRMED BY PAIRED READS FUSIONS AND MISASSEMBLIES:
            #elif args.fusion_misassemble_analyze:
            #    #print 'analyze fusions and misassemblies by paired reads...'
            #    is_fusion, is_misassemble = \
            #        self.get_confirmed_fusion_misassemblies(best_union_alignments, fusion_err_threshhold)

            # UPDATE FILES FOR FUSIONS AND MISASSEMBLIES:
            transcripts_metrics_report.update_psl_union_fus_mis_files(is_fusion, is_misassemble, best_union_lines, single_transcript_lines)

    def print_fusion_misassemble_metrics(self, args, path_txt_simple, logger):
        logger.info('      Getting Fusions and misassemblies report...')

        with open(path_txt_simple, 'w') as fout:
            fout.write('METRICS OF FUSIONS AND MISASSEMBLIES:\n')
            #if args.fusion_misassemble_analyze:
            #    fout.write('{:<100}'.format('Misassembly candidates') + str(self.misassemble_num) + '\n\n')

            #    fout.write('{:<100}'.format('Fusion candidates') + str(self.fusion_num) + '\n\n')
            #else:
            fout.write('{:<100}'.format('Misassembly candidates') + str(self.num_misassembled) + '\n\n')

        logger.info('        saved to {}'.format(path_txt_simple))