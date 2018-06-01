__author__ = 'letovesnoi'

from datetime import datetime
import os

from general import UtilsGeneral
from general import best_alignment_set

from objects import Alignment

from quast_libs import fastaparser


class AlignmentsReport():
    """Class which generate alignments report"""

    def __init__(self, psl_file, blast6_file, label, tmp_dir):
        # Blast6s FILES (FOR MISASSEMBLIES SEARCH):
        self.blast6_report = None
        if blast6_file != None:
            self.blast6_report = self.Blast6AlignmentsReport(label, tmp_dir)

        # PSLs FILES:
        self.blat_report = None
        if psl_file != None:
            self.blat_report = self.BlatPSLAlignmentsReport(label, tmp_dir)

    @classmethod
    def get_alignments_report(cls, label, psl_file, blast6_file, transcripts_dict, tmp_dir, min_alignment_threshold,
                              logger, ALIGNMENT_THRESHOLDS):
        alignment_report = cls(psl_file, blast6_file, label, tmp_dir)

        if alignment_report.blast6_report != None:
            alignment_report.blast6_report.get_blast6_alignments_report(label, blast6_file, tmp_dir, min_alignment_threshold, logger, ALIGNMENT_THRESHOLDS)

        if alignment_report.blat_report != None:
            alignment_report.blat_report.get_psl_alignments_report(label, psl_file, transcripts_dict, tmp_dir, min_alignment_threshold, logger, ALIGNMENT_THRESHOLDS)

        return alignment_report

    @classmethod
    def update_alignment_file(cls, single_transcript_lines, fout_file):
        for line in single_transcript_lines:
            fout_file.write(line + '\n')


    class Blast6AlignmentsReport():

        def __init__(self, label, tmp_dir):
            # single alignment file (contains only best alignments for assembled transcripts),
            self.assembled_blast6_file = os.path.join(tmp_dir, '{}.assembled.best.blast6'.format(label))

            # paralogs alignment file (contains only best alignments for assembled transcripts having paralogs)
            self.paralogous_blast6_file = os.path.join(tmp_dir, '{}.paralogs.best.blast6'.format(label))

            # multiple union file (contains only best alignments for misassembled transcripts):
            self.misassembled_blast6_union_file = os.path.join(tmp_dir, '{}.misassembled.best.blast6'.format(label))
            # contains all alignments for misassembled transcripts:
            self.misassembled_blast6_file = os.path.join(tmp_dir, '{}.misassembled.blast6'.format(label))


        # create files with psl alignments report:
        @classmethod
        def get_blast6_alignments_report(cls, label, blast6_file, output_dir, args_min_alignment, logger, ALIGNMENT_THRESHOLDS):
            start_union_time = datetime.now()

            # CREATE TEMPORARY ALIGNMENTS REPORT:
            blast6_report = cls(label, output_dir)

            # open files for report:
            fout_assembled_psl = open(blast6_report.assembled_blast6_file, 'w')
            fout_paralogous_psl = open(blast6_report.paralogous_blast6_file, 'w')
            fout_misassembled_union = open(blast6_report.misassembled_blast6_union_file, 'w')
            fout_misassembled_psl = open(blast6_report.misassembled_blast6_file, 'w')

            logger.print_timestamp('  ')
            logger.info('  Getting BLAST alignments report files...')

            with open(blast6_file, 'r') as fin:
                line1 = fin.readline().strip()
                line2 = fin.readline().strip()
                while line1 != '':
                    single_transcript_lines, single_transcript_alignments, line1, line2 = \
                        get_curr_single_transcript_lines_alignments(args_min_alignment, line1, line2, fin, False)

                    if len(single_transcript_lines) == 0 and len(single_transcript_alignments) == 0:
                        #logger.debug('    Skipping alignments...')
                        continue

                    # GET UNION ALIGNMENTS:
                    best_union_alignments = best_alignment_set.get_best_alignment_set(single_transcript_alignments, ALIGNMENT_THRESHOLDS)
                    best_union_lines = get_best_lines_set(best_union_alignments)

                    # choose over single transcript alignments best union alignments and them lines:
                    # (it can be multiple union alignment cause misassemble or list of single union alignments cause paralogs)
                    # FOR MISASSEMBLIES:
                    if len(best_union_lines) > 1 and len(best_union_alignments) > 1:
                        # UPDATE ALIGNMENTS FILES:
                        AlignmentsReport.update_alignment_file(best_union_lines, fout_misassembled_union)
                        AlignmentsReport.update_alignment_file(single_transcript_lines, fout_misassembled_psl)

                    # FOR ASSEMBLED TRANSCRIPT:
                    elif len(best_union_lines) == 1 and len(best_union_alignments) == 1:
                        is_paralogous = False

                        curr_single_transcript_alignments = single_transcript_alignments[:]

                        best_single_score = best_union_alignments[0].bitscore

                        curr_single_score = best_single_score

                        while best_single_score == curr_single_score and len(best_union_alignments) == 1 and curr_single_transcript_alignments != []:
                            # UPDATE BEST ALIGNMENTS FILE:
                            AlignmentsReport.update_alignment_file(best_union_lines, fout_assembled_psl)

                            prev_best_union_lines = best_union_lines

                            curr_single_transcript_alignments.remove(best_union_alignments[0])

                            if curr_single_transcript_alignments != []:
                                # GET UNION ALIGNMENTS:
                                best_union_alignments = best_alignment_set.get_best_alignment_set(curr_single_transcript_alignments, ALIGNMENT_THRESHOLDS)
                                best_union_lines = get_best_lines_set(best_union_alignments)

                                curr_single_score = best_union_alignments[0].bitscore

                            if best_single_score == curr_single_score and len(best_union_alignments) == 1 and curr_single_transcript_alignments != []:
                                is_paralogous = True

                            if is_paralogous == True:
                                AlignmentsReport.update_alignment_file(prev_best_union_lines, fout_paralogous_psl)


            logger.info('    saved to ' + blast6_report.assembled_blast6_file + ' (contains best alignments for assembled transcripts)\n' +
                        ' ' * 13 + blast6_report.paralogous_blast6_file + ' (contains best alignments for assembled transcripts having paralogs)\n' +
                        ' ' * 13 + blast6_report.misassembled_blast6_union_file + ' (contains best alignments for misassembled transcripts)\n' +
                        ' ' * 13 + blast6_report.misassembled_blast6_file + ' (contains all alignments for misassembled transcripts)')

            # close file for report:
            fout_assembled_psl.close()
            fout_paralogous_psl.close()
            fout_misassembled_union.close()
            fout_misassembled_psl.close()

            end_union_time = datetime.now()
            elapsed_union_time = end_union_time - start_union_time

            logger.debug('  ELAPSED TIME: ' + str(elapsed_union_time))

            return blast6_report


    class BlatPSLAlignmentsReport():

        def __init__(self, label, tmp_dir):
            # remove alignments with negative query gaps lengths:
            self.wout_strange_plsl_file = os.path.join(tmp_dir, '{}.unstrange.psl'.format(label))

            # alignment file contains only best alignments for assembled transcripts,
            self.assembled_psl_file = os.path.join(tmp_dir, '{}.assembled.best.psl'.format(label))

            # alignment file contains only best alignments for uniquely aligned assembled transcripts),
            self.uniquely_psl_file = os.path.join(tmp_dir, '{}.uniquely.psl').format(label)

            # paralogs alignment file (contains only best alignments for assembled transcripts having paralogs)
            self.paralogous_psl_file = os.path.join(tmp_dir, '{}.paralogs.best.psl'.format(label))

            # multiple union file (contains only best alignments for misassembled transcripts):
            self.misassembled_psl_union_file = os.path.join(tmp_dir, '{}.misassembled.best.psl'.format(label))
            # contains all alignments for misassembled transcripts:
            self.misassembled_psl_file = os.path.join(tmp_dir, '{}.misassembled.psl'.format(label))

            # TEMPORARY FILES FOR TESTING:
            self.psl_wout_cross_path = os.path.join(tmp_dir, '{}.woc.psl'.format(label))

            # file with union transcript best alignments which worse splitted by blat:
            self.fake_blat_file = os.path.join(tmp_dir, '{}.fake_blat'.format(label))

            # file with splitted large query gaps in query:
            self.psl_wolqg_file = os.path.join(tmp_dir, '{}.wolqg.psl'.format(label))

            # file with low complexity best alignments:
            self.low_complexity_file = os.path.join(tmp_dir, '{}.low_complexity.fasta'.format(label))


        # create files with psl alignments report:
        @classmethod
        def get_psl_alignments_report(cls, label, psl_file, transcripts_dict, tmp_dir, args_min_alignment, logger, ALIGNMENT_THRESHOLDS):
            start_union_time = datetime.now()

            # CREATE TEMPORARY ALIGNMENTS REPORT:
            blat_report = cls(label, tmp_dir)

            # blat_report.wout_strange_plsl_file = remove_cross_blocks(psl_file, blat_report.wout_strange_plsl_file)

            # open files for report:
            fout_assembled_psl = open(blat_report.assembled_psl_file, 'w')
            fout_uniquely_file = open(blat_report.uniquely_psl_file, 'w')
            fout_paralogous_psl = open(blat_report.paralogous_psl_file, 'w')
            fout_misassembled_union = open(blat_report.misassembled_psl_union_file, 'w')
            fout_misassembled_psl = open(blat_report.misassembled_psl_file, 'w')
            fout_fake_blat = open(blat_report.fake_blat_file, 'w')
            fout_psl_wolqg = open(blat_report.psl_wolqg_file, 'w')
            # fout_low_complexity = open(blat_report.low_complexity_file, 'w')

            logger.print_timestamp('  ')
            logger.info('  Getting GMAP (or BLAT) alignments report files...')

            blat_report.psl_wout_cross_path = \
                remove_strange_psl_alignments(psl_file, blat_report.psl_wout_cross_path)

            with open(blat_report.psl_wout_cross_path, 'r') as fin:
                line1 = fin.readline().strip()
                line2 = fin.readline().strip()
                while line1 != '':
                    single_transcript_lines, single_transcript_alignments, line1, line2 = \
                        get_curr_single_transcript_lines_alignments(args_min_alignment, line1, line2, fin, True)

                    if len(single_transcript_lines) == 0 and len(single_transcript_alignments) == 0:
                        #logger.debug('    Skipping alignments...')
                        continue

                    # SPLIT ALIGNMENTS WITH LARGE QUERY GAPS:
                    single_transcript_lines, single_transcript_alignments = \
                        split_large_query_gap_single_transcript_alignments(single_transcript_alignments, fout_psl_wolqg,
                                                                           ALIGNMENT_THRESHOLDS.QUERY_GAP_THRESHOLD,
                                                                           ALIGNMENT_THRESHOLDS.MIN_SPLIT_ALIGNMENT_THRESHOLD)

                    if len(single_transcript_lines) == 0 and len(single_transcript_alignments) == 0:
                        #logger.debug('    Skipping alignments...')
                        continue

                    # GET UNION ALIGNMENTS:
                    best_union_alignments = best_alignment_set.get_best_alignment_set(single_transcript_alignments, ALIGNMENT_THRESHOLDS)
                    best_union_lines = get_best_lines_set(best_union_alignments)

                    # GET UNION FAKE BLAT ALIGNMENTS:
                    best_union_lines, best_union_alignments, single_transcript_lines, single_transcript_alignments = \
                        get_union_fake_blat_alignments(best_union_lines, best_union_alignments, single_transcript_lines,
                                                       single_transcript_alignments, fout_fake_blat, ALIGNMENT_THRESHOLDS)

                    # REMOVE LOW COMPLEXITY TAILS:
                    best_union_lines, best_union_alignments, single_transcript_lines, single_transcript_alignments = \
                        remove_low_complexity(best_union_lines, best_union_alignments, single_transcript_lines,
                                              single_transcript_alignments, transcripts_dict, blat_report.low_complexity_file,
                                              ALIGNMENT_THRESHOLDS.LOW_COMPLEXITY_LEN_THRESHOLD)

                    # choose over single transcript alignments best union alignments and them lines:
                    # (it can be multiple union alignment cause misassemble or list of single union alignments cause paralogs)
                    # FOR MISASSEMBLIES:
                    if len(best_union_lines) > 1 and len(best_union_alignments) > 1:
                        # UPDATE ALIGNMENTS FILES:
                        AlignmentsReport.update_alignment_file(best_union_lines, fout_misassembled_union)
                        AlignmentsReport.update_alignment_file(single_transcript_lines, fout_misassembled_psl)

                    # FOR ASSEMBLED TRANSCRIPT:
                    elif len(best_union_lines) == 1 and len(best_union_alignments) == 1:
                        is_paralogous = False

                        curr_single_transcript_alignments = single_transcript_alignments[:]
                        curr_single_transcript_lines = single_transcript_lines[:]

                        best_single_score = best_union_alignments[0].matches

                        curr_single_score = best_single_score

                        prev_best_union_lines = best_union_lines

                        while best_single_score == curr_single_score and len(best_union_alignments) == 1 and curr_single_transcript_alignments != []:
                            # UPDATE BEST ALIGNMENTS FILE:
                            AlignmentsReport.update_alignment_file(best_union_lines, fout_assembled_psl)

                            prev_best_union_lines = best_union_lines

                            curr_single_transcript_lines.remove(best_union_lines[0])
                            curr_single_transcript_alignments.remove(best_union_alignments[0])

                            if curr_single_transcript_lines != [] and curr_single_transcript_alignments != []:
                                # GET UNION ALIGNMENTS:
                                best_union_alignments = best_alignment_set.get_best_alignment_set(curr_single_transcript_alignments, ALIGNMENT_THRESHOLDS)
                                best_union_lines = get_best_lines_set(best_union_alignments)

                                # GET UNION FAKE BLAT ALIGNMENTS:
                                best_union_lines, best_union_alignments, curr_single_transcript_lines, curr_single_transcript_alignments = \
                                    get_union_fake_blat_alignments(best_union_lines, best_union_alignments, curr_single_transcript_lines,
                                                                   curr_single_transcript_alignments, fout_fake_blat, ALIGNMENT_THRESHOLDS)

                                curr_single_score = best_union_alignments[0].matches

                                # REMOVE LOW COMPLEXITY TAILS:
                                best_union_lines, best_union_alignments, curr_single_transcript_lines, curr_single_transcript_alignments = \
                                    remove_low_complexity(best_union_lines, best_union_alignments, curr_single_transcript_lines,
                                                          curr_single_transcript_alignments, transcripts_dict, blat_report.low_complexity_file,
                                                          ALIGNMENT_THRESHOLDS.LOW_COMPLEXITY_LEN_THRESHOLD)

                            if best_single_score == curr_single_score and len(best_union_alignments) == 1 and curr_single_transcript_alignments != []:
                                is_paralogous = True

                            if is_paralogous == True:
                                AlignmentsReport.update_alignment_file(prev_best_union_lines, fout_paralogous_psl)

                        if is_paralogous == False:
                            AlignmentsReport.update_alignment_file(prev_best_union_lines, fout_uniquely_file)


            logger.info('    saved to ' + blat_report.assembled_psl_file + ' (contains best alignments for assembled transcripts)\n' +
                        ' ' * 13 + blat_report.uniquely_psl_file + ' (contains best alignments for uniquely aligned assembled transcripts)\n' +
                        ' ' * 13 + blat_report.paralogous_psl_file + ' (contains best alignments for assembled transcripts having paralogs)\n' +
                        ' ' * 13 + blat_report.misassembled_psl_union_file + ' (contains best alignments for misassembled transcripts)\n' +
                        ' ' * 13 + blat_report.misassembled_psl_file + ' (contains all alignments for misassembled transcripts)')

            # close file for report:
            fout_assembled_psl.close()
            fout_paralogous_psl.close()
            fout_misassembled_union.close()
            fout_misassembled_psl.close()
            fout_fake_blat.close()
            fout_psl_wolqg.close()
            # fout_low_complexity.close()

            end_union_time = datetime.now()
            elapsed_union_time = end_union_time - start_union_time

            logger.debug('  ELAPSED TIME: ' + str(elapsed_union_time))

            return blat_report


def get_best_lines_set(alignments):
    lines = []
    for a_i in alignments:
        if a_i.format == 'psl':
            lines.append(a_i.get_psl_line_from_alignment())
        elif a_i.format == 'blast6':
            lines.append(a_i.line)
    return lines


# union alignments crossing at most err_cross in query and distant no more than err_space in query:
# for determine misassemblies find best union with score define by union_penalty for non close (define as fake blat) alignments:
'''def get_greedily_union(single_transcript_lines, single_transcript_alignments, is_psl_lines, err_cross_query_fake, err_space_query_fake,
                       err_cross_target_fake, err_space_target_fake, union_penalty=UNION_PENALTY, err_cross_union=ERR_CROSS_UNION):
    #logger.debug('      Getting best union alignments...')
    q_ends = []
    for i in range(len(single_transcript_alignments)):
        q_ends.append(single_transcript_alignments[i].query_fragment.end)
    q_ends_sort_index, q_ends_sort_array = UtilsGeneral.get_order_indexes_elements(q_ends)

    score = {}
    curr = 0
    last_alignment = {}
    last_line = {}

    q_positions = []
    for i in range(len(single_transcript_alignments)):
        q_positions.append(single_transcript_alignments[i].query_fragment.start + err_cross_union)
        q_positions.append(single_transcript_alignments[i].query_fragment.end)
    q_positions.sort()

    score[q_positions[0]] = 0
    last_alignment[q_positions[0]] = []
    last_line[q_positions[0]] = []
    for i_pos in range(1, len(q_positions)):
        if q_positions[i_pos] in score:
            tmp_score = max(score[q_positions[i_pos - 1]], score[q_positions[i_pos]])
        else:
            tmp_score = score[q_positions[i_pos - 1]]

        if q_positions[i_pos] == q_ends_sort_array[curr]:
            cross = 0
            tmp_union_penalty = union_penalty
            alignment1 = single_transcript_alignments[q_ends_sort_index[curr]]
            line1 = single_transcript_lines[q_ends_sort_index[curr]]
            if len(last_alignment[single_transcript_alignments[q_ends_sort_index[curr]].query_fragment.start + err_cross_union]) > 0:
                alignment0 = last_alignment[single_transcript_alignments[q_ends_sort_index[curr]].query_fragment.start + err_cross_union][-1]
                cross = max(0, alignment0.query_fragment.end - alignment1.query_fragment.start + 1)
                if best_alignment_set.is_union_fake_blat(alignment0, alignment1, err_cross_query_fake, err_space_query_fake,
                                                         err_cross_target_fake, err_space_target_fake) == True and is_psl_lines == True:
                    tmp_union_penalty = 0
            else:
                tmp_union_penalty = 0

            if is_psl_lines == True:
                curr_score = alignment1.matches
            else:
                curr_score = alignment1.bitscore

            union_score = score[alignment1.query_fragment.start + err_cross_union] + curr_score - tmp_union_penalty - cross
            tmp_score = max(tmp_score, union_score)
            curr += 1

        if tmp_score == score[q_positions[i_pos - 1]]:
            score[q_positions[i_pos]] = tmp_score
            last_alignment[q_positions[i_pos]] = last_alignment[q_positions[i_pos - 1]]
            last_line[q_positions[i_pos]] = last_line[q_positions[i_pos - 1]]
        elif tmp_score == union_score:
            score[q_positions[i_pos]] = tmp_score
            last_alignment[q_positions[i_pos]] = last_alignment[alignment1.query_fragment.start + err_cross_union] + [alignment1]
            last_line[q_positions[i_pos]] = last_line[alignment1.query_fragment.start + err_cross_union] + [line1]

    max_score = 0
    max_key = q_positions[q_positions.index(score.keys()[0])]
    for key in score:
        if score[key] > max_score:
            max_score = score[key]
            max_key = key

    union_alignments = last_alignment[max_key]
    union_lines = last_line[max_key]
    #logger.debug('      Done.')
    return union_lines, union_alignments'''


# union close (cross at most err_cross, distance no more than err_space, at one strand and chromosome) alignments:
def get_union_fake_blat_alignments(union_lines, union_alignments, single_transcript_lines, single_transcript_alignments,
                                   fout_fake_blat, ALIGNMENT_THRESHOLDS):
    #logger.debug('      Getting fake blat alignments...')
    new_union_alignments = []
    new_union_lines = []
    q_starts = []
    for alignment in union_alignments:
        q_starts.append(alignment.query_fragment.start)
    sort_index, sort_array = UtilsGeneral.get_order_indexes_elements(q_starts)

    i_alignment0 = sort_index[0]
    alignment0 = union_alignments[i_alignment0]
    alignment_line0 = union_lines[i_alignment0]
    if len(sort_index) == 1:
        new_union_alignments.append(alignment0)
        new_union_lines.append(alignment_line0)
    for i in range(len(sort_index) - 1):
        i_alignment1 = sort_index[i + 1]
        alignment1 = union_alignments[i_alignment1]
        alignment_line1 = union_lines[i_alignment1]
        if best_alignment_set.is_union_fake_blat(alignment0, alignment1, ALIGNMENT_THRESHOLDS):
            fout_fake_blat.write(alignment_line0 + '\n' + alignment_line1 + '\n\n\n')
            alignment0 = get_union_fake_blat_alignment(alignment0, alignment1)
            alignment_line0 = alignment0.get_psl_line_from_alignment()
        else:
            new_union_alignments.append(alignment0)
            new_union_lines.append(alignment_line0)
            i_alignment0 = sort_index[i + 1]
            alignment0 = union_alignments[i_alignment0]
            alignment_line0 = union_lines[i_alignment0]
        if i + 1 == len(sort_index) - 1:
            new_union_alignments.append(alignment0)
            new_union_lines.append(alignment_line0)

    # update single transcript lines / alignments:
    if len(union_lines) != len(new_union_lines):
        for i_alignment in range(len(union_lines)):
            if union_lines[i_alignment] not in new_union_lines:
                single_transcript_lines.remove(union_lines[i_alignment])
                single_transcript_alignments.remove(union_alignments[i_alignment])
        for i_alignment in range(len(new_union_lines)):
            if new_union_lines[i_alignment] not in single_transcript_lines:
                single_transcript_lines.append(new_union_lines[i_alignment])
                single_transcript_alignments.append(new_union_alignments[i_alignment])

    #logger.debug('      Done.')
    return new_union_lines, new_union_alignments, single_transcript_lines, single_transcript_alignments


# alignments are union if they are cross at most err_cross, distance no more than err_space in query, are at one strand and chromosome:
# def is_union_fake_blat(alignment0, alignment1, err_space, err_cross):
#     if alignment0.strand == '+':
#         if alignment0.target_fragment.end < alignment1.target_fragment.start + err_cross and alignment0.strand == alignment1.strand and alignment0.target_fragment.name == alignment1.target_fragment.name and \
#                 ((alignment1.query_fragment.start - alignment0.query_fragment.end - 1 <= err_space and alignment1.query_fragment.start - alignment0.query_fragment.end - 1 >= 0) or
#                      (alignment0.query_fragment.end - alignment1.query_fragment.start + 1 <= err_cross and alignment0.query_fragment.end - alignment1.query_fragment.start + 1 > 0)):
#             return True
#     else:
#         if alignment1.target_fragment.end < alignment0.target_fragment.start + err_cross and alignment0.strand == alignment1.strand and alignment0.target_fragment.name == alignment1.target_fragment.name and \
#                 ((alignment1.query_fragment.start - alignment0.query_fragment.end - 1 <= err_space and alignment1.query_fragment.start - alignment0.query_fragment.end - 1 >= 0) or
#                      (alignment0.query_fragment.end - alignment1.query_fragment.start + 1 <= err_cross and alignment0.query_fragment.end - alignment1.query_fragment.start + 1 > 0)):
#             return True
#     return False


# def is_union_fake_blastn(alignment0, alignment1, err_cross, err_space):
#     if is_union_fake_blat(alignment0, alignment1, err_cross, err_space):
#         if alignment0.strand == '+':
#             if alignment1.target_fragment.start - alignment0.target_fragment.end - 1 <= err_space:
#                 return True
#         else:
#             if alignment0.target_fragment.start - alignment1.target_fragment.end - 1 <= err_space:
#                 return True
#     else:
#         return False


def get_union_fake_blat_alignment(alignment0, alignment1):
    union_alignment = Alignment.PSLFileAlignment()
    strand = alignment0.strand

    # for cross more than one block or equal one block:
    if strand == '+':
        tmp_start = alignment1.query_fragment.starts[0]
        tmp_array = alignment0.query_fragment.starts
    else:
        tmp_start = alignment0.query_fragment.starts[0]
        tmp_array = alignment1.query_fragment.starts
    i_block1 = UtilsGeneral.get_bin_search_position_of_element(tmp_array, tmp_start)

    if strand == '+':
        if i_block1 != alignment0.blocks_num:
            alignment0 = alignment0.get_split_alignment(0, i_block1 - 1)
    else:
        if i_block1 != alignment1.blocks_num:
            alignment1 = alignment1.get_split_alignment(0, i_block1 - 1)

    tmp_qbase_cross = alignment0.query_fragment.end - alignment1.query_fragment.start + 1
    if strand == '+':
        tmp_tbase_cross = alignment0.target_fragment.end - alignment1.target_fragment.start + 1
    else:
        tmp_tbase_cross = alignment1.target_fragment.end - alignment0.target_fragment.start + 1
    cross_bases = max(tmp_qbase_cross, tmp_tbase_cross, 0)

    f_together = False
    if tmp_qbase_cross >= 0 and tmp_tbase_cross >= 0 and tmp_qbase_cross == tmp_tbase_cross:
        f_together = True

    union_alignment.matches = alignment0.matches + alignment1.matches - max(cross_bases, 0)
    union_alignment.mismatches = alignment0.mismatches + alignment1.mismatches
    union_alignment.repmatches = alignment0.repmatches + alignment1.repmatches
    union_alignment.n_num = alignment0.n_num + alignment1.n_num
    union_alignment.strand = strand
    if f_together:
        union_alignment.blocks_num = alignment0.blocks_num + alignment1.blocks_num - 1
    else:
        union_alignment.blocks_num = alignment0.blocks_num + alignment1.blocks_num

    if strand == '+':
        if f_together:
            union_alignment.blocks_sizes = alignment0.blocks_sizes[:-1] + \
                                        [alignment0.blocks_sizes[-1] + alignment1.blocks_sizes[0] - max(cross_bases, 0)] + \
                                        alignment1.blocks_sizes[1:]
        else:
            union_alignment.blocks_sizes = alignment0.blocks_sizes[:-1] + [alignment0.blocks_sizes[-1] - max(cross_bases, 0)] + alignment1.blocks_sizes[:]
    else:
        if f_together:
            union_alignment.blocks_sizes = alignment1.blocks_sizes[:-1] + \
                                        [alignment1.blocks_sizes[-1] + alignment0.blocks_sizes[0] - max(cross_bases, 0)] + \
                                        alignment0.blocks_sizes[1:]
        else:
            union_alignment.blocks_sizes = alignment1.blocks_sizes[:-1] + [alignment1.blocks_sizes[-1] - max(cross_bases, 0)] + alignment0.blocks_sizes[:]

    set_union_query_fragment_attributes(alignment0.query_fragment, alignment1.query_fragment, union_alignment.query_fragment, strand, f_together, cross_bases)

    set_union_target_fragment_attributes(alignment0.target_fragment, alignment1.target_fragment, union_alignment.target_fragment, strand, f_together, cross_bases)

    return union_alignment


def set_union_query_fragment_attributes(query_fragment0, query_fragment1, union_query_fragment, strand, f_together, cross_bases):
    tmp_qbase_ins = query_fragment1.start - query_fragment0.end - 1
    tmp_qbase_cross = query_fragment0.end - query_fragment1.start + 1
    # Number of inserts in query
    if tmp_qbase_ins > 0:
        tmp_num_ins = 1
    else:
        tmp_num_ins = 0
    union_query_fragment.num_insert = query_fragment0.num_insert + query_fragment1.num_insert + tmp_num_ins
    # Number of bases inserted in query
    union_query_fragment.base_insert = query_fragment0.base_insert + query_fragment1.base_insert + max(tmp_qbase_ins, 0)
    # Query sequence name
    union_query_fragment.name = query_fragment0.name
    # Query sequence size
    union_query_fragment.size = query_fragment0.size
    # Alignment start position in query
    union_query_fragment.start = query_fragment0.start
    # Alignment end position in query
    union_query_fragment.end = query_fragment1.end
    # Comma-separated list of starting positions of each block in query
    if strand == '+':
        if f_together:
            union_query_fragment.starts = query_fragment0.starts[:] + query_fragment1.starts[1:]
            union_query_fragment.ends = query_fragment0.ends[:-1] + query_fragment1.ends[:]
        else:
            union_query_fragment.starts = query_fragment0.starts[:] + query_fragment1.starts[:]
            union_query_fragment.ends = query_fragment0.ends[:-1] + [query_fragment0.ends[-1] - max(cross_bases, 0)] + query_fragment1.ends[:]
    else:
        if f_together:
            union_query_fragment.starts = query_fragment1.starts[:] + query_fragment0.starts[1:]
            union_query_fragment.ends = query_fragment1.ends[:-1] + query_fragment0.ends[:]
        else:
            union_query_fragment.starts = query_fragment1.starts[:] + query_fragment0.starts[:]
            union_query_fragment.ends = query_fragment1.ends[:-1] + [query_fragment1.ends[-1] - max(cross_bases, 0)] + query_fragment0.ends[:]


def set_union_target_fragment_attributes(target_fragment0, target_fragment1, union_target_fragment, strand, f_together, cross_bases):
    if strand == '+':
        tmp_tbase_ins = target_fragment1.start - target_fragment0.end - 1
        tmp_tbase_cross = target_fragment0.end - target_fragment1.start + 1
    else:
        tmp_tbase_ins = target_fragment0.start - target_fragment1.end - 1
        tmp_tbase_cross = target_fragment1.end - target_fragment0.start + 1
    # Number of inserts in query
    if tmp_tbase_ins > 0:
        tmp_num_ins = 1
    else:
        tmp_num_ins = 0
    #  Number of inserts in target
    union_target_fragment.num_insert = target_fragment0.num_insert + target_fragment1.num_insert + tmp_num_ins
    # Number of bases inserted in target
    union_target_fragment.base_insert = target_fragment0.base_insert + target_fragment1.base_insert + max(0, tmp_tbase_ins)
    # Target sequence name
    union_target_fragment.name = target_fragment0.name
    # Target sequence size
    union_target_fragment.size = target_fragment0.size
    # Comma-separated list of starting positions of each block in target
    if strand == '+':
        # Alignment start position in target
        union_target_fragment.start = target_fragment0.start
        # Alignment end position in target
        union_target_fragment.end = target_fragment1.end
        if f_together:
            union_target_fragment.starts = target_fragment0.starts[:] + target_fragment1.starts[1:]
            union_target_fragment.ends = target_fragment0.ends[:-1] + target_fragment1.ends[:]
        else:
            union_target_fragment.starts = target_fragment0.starts[:] + target_fragment1.starts[:]
            union_target_fragment.ends = target_fragment0.ends[:-1] + [target_fragment0.ends[-1] - max(cross_bases, 0)] + target_fragment1.ends[:]
    else:
        # Alignment start position in target
        union_target_fragment.start = target_fragment1.start
        # Alignment end position in target
        union_target_fragment.end = target_fragment0.end
        if f_together:
            union_target_fragment.starts = target_fragment1.starts[:] + target_fragment0.starts[1:]
            union_target_fragment.ends = target_fragment1.ends[:-1] + target_fragment0.ends[:]
        else:
            union_target_fragment.starts = target_fragment1.starts[:] + target_fragment0.starts[:]
            union_target_fragment.ends = target_fragment1.ends[:-1] + [target_fragment1.ends[-1] - max(cross_bases, 0)] + target_fragment0.ends[:]


# split alignments containing large gap between blocks to alignments containing at most query_gap_threshold gaps between blocks:
def split_large_query_gap_single_transcript_alignments(single_transcript_alignments, fout_psl,
                                                       query_gap_threshold, min_split_alignment_threshold):
    new_single_transcript_alignments = []
    new_single_transcript_lines = []
    for alignment in single_transcript_alignments:
        tmp_count_split = 0
        curr_iBlock = -1
        for i_block in range(alignment.blocks_num - 1):
            if alignment.query_fragment.starts[i_block + 1] - alignment.query_fragment.ends[i_block] > query_gap_threshold:
                tmp_count_split += 1
                split_alignment = alignment.get_split_alignment(curr_iBlock + 1, i_block)
                # filter small alignments:
                if sum(split_alignment.blocks_sizes) >= min_split_alignment_threshold:
                    new_single_transcript_alignments.append(split_alignment)

                    split_line = split_alignment.get_psl_line_from_alignment()
                    new_single_transcript_lines.append(split_line)
                    fout_psl.write(split_line + '\n')
                curr_iBlock = i_block
        if tmp_count_split != 0:
            split_alignment = alignment.get_split_alignment(curr_iBlock + 1, alignment.blocks_num - 1)
            # filter small alignments:
            if sum(split_alignment.blocks_sizes) >= min_split_alignment_threshold:
                new_single_transcript_alignments.append(split_alignment)

                split_line = split_alignment.get_psl_line_from_alignment()
                new_single_transcript_lines.append(split_line)
                fout_psl.write(split_line + '\n')
        else:
            new_single_transcript_alignments.append(alignment)

            line = alignment.get_psl_line_from_alignment()
            new_single_transcript_lines.append(line)
            fout_psl.write(line + '\n')
    return new_single_transcript_lines, new_single_transcript_alignments


# remove low complexity regions cause polyA/T or something else:
def remove_low_complexity(union_lines, union_alignments, single_transcript_lines, single_transcript_alignments,
                          transcript_dict, out_low_complexity_file, low_complexity_len_threshold):
    transcript_seq = transcript_dict[union_alignments[0].query_fragment.name]

    clear_union_lines, clear_union_alignments, single_transcript_lines, single_transcript_alignments = \
        remove_low_complexity_tail(union_lines, union_alignments, single_transcript_lines, single_transcript_alignments,
                                   transcript_seq, out_low_complexity_file, low_complexity_len_threshold, end_tail=True)

    clear_union_lines, clear_union_alignments, single_transcript_lines, single_transcript_alignments = \
        remove_low_complexity_tail(clear_union_lines, clear_union_alignments, single_transcript_lines, single_transcript_alignments,
                                   transcript_seq, out_low_complexity_file, low_complexity_len_threshold, end_tail=False)

    return clear_union_lines, clear_union_alignments, single_transcript_lines, single_transcript_alignments


# remove low complexity tail:
def remove_low_complexity_tail(union_lines, union_alignments, single_transcript_lines, single_transcript_alignments,
                               transcript_seq, out_low_complexity_file, threshold_block_len, end_tail):
    block_seq = ''
    for i in range(len(union_alignments)):
        if end_tail == True:
            i_alignment = len(union_alignments) - 1 - i
        else:
            i_alignment = i
        psl_alignment = union_alignments[i_alignment]
        for j in range(psl_alignment.blocks_num):
            if end_tail == True:
                i_block = psl_alignment.blocks_num - 1 - j
            else:
                i_block = j

            start = psl_alignment.query_fragment.starts[i_block]
            end = psl_alignment.query_fragment.ends[i_block]

            if psl_alignment.strand == '+':
                block_seq += transcript_seq[start:end + 1]
            else:
                block_seq += UtilsGeneral.rev_comp(transcript_seq)[start:end + 1]

            if len(block_seq) < threshold_block_len:
                continue

            if not is_low_complexity(block_seq):
                if (end_tail == True and i_block == psl_alignment.blocks_num - 1) or (end_tail == False and i_block == 0):
                    latest_alignment = psl_alignment
                    latest_line = union_lines[i_alignment]
                else:
                    if end_tail == True:
                        latest_alignment = psl_alignment.get_split_alignment(0, i_block)
                    else:
                        latest_alignment = psl_alignment.get_split_alignment(i_block, psl_alignment.blocks_num - 1)
                    latest_line = latest_alignment.get_psl_line_from_alignment()

                if end_tail == True:
                    clear_union_alignments = union_alignments[:i_alignment] + [latest_alignment]
                    clear_union_lines = union_lines[:i_alignment] + [latest_line]
                else:
                    clear_union_alignments = [latest_alignment] + union_alignments[i_alignment + 1:]
                    clear_union_lines = [latest_line] + union_lines[i_alignment + 1:]

                # update single transcript lines / alignments:
                for i_alignment in range(len(union_lines)):
                    if union_lines[i_alignment] not in clear_union_lines:
                        single_transcript_lines.remove(union_lines[i_alignment])
                        single_transcript_alignments.remove(union_alignments[i_alignment])
                for i_alignment in range(len(clear_union_lines)):
                    if clear_union_lines[i_alignment] not in single_transcript_lines:
                        single_transcript_lines.append(clear_union_lines[i_alignment])
                        single_transcript_alignments.append(clear_union_alignments[i_alignment])

                return clear_union_lines, clear_union_alignments, single_transcript_lines, single_transcript_alignments,

            else:
                fastaparser.write_fasta(out_low_complexity_file, [('{}_block'.format(psl_alignment.query_fragment.name), block_seq)], mode='a')

            if len(block_seq) > threshold_block_len:
                block_seq = ''

    # update single transcript lines / alignments:
    for i_alignment in range(len(union_lines)):
        single_transcript_lines.remove(union_lines[i_alignment])
        single_transcript_alignments.remove(union_alignments[i_alignment])

    return [], [], single_transcript_lines, single_transcript_alignments


# determine low complexity regions as regions containing mainly two types of kmers:
def is_low_complexity(seq, k=4):
    k_spect = {}
    for i in range(len(seq) - k + 1):
        if seq[i:i + k] not in k_spect:
            k_spect[seq[i:i + k]] = 0
        k_spect[seq[i:i + k]] += 1

    max_key = k_spect.keys()[0]
    max_value = k_spect[k_spect.keys()[0]]
    for key in k_spect:
        if k_spect[key] > max_value:
            max_value = k_spect[key]
            max_key = key
    majority_num = max_value
    k_spect[max_key] = 0
    majority_num += max(k_spect.values())
    if majority_num > 0.6 * (len(seq) - k + 1):
        return True
    else:
        return False


# extract from psl file single transcript lines and get them alignments:
def get_curr_single_transcript_lines_alignments(min_alignment_threshold, line1, line2, fin, is_psl_lines):
    single_transcript_lines = []
    single_transcript_alignments = []

    # get alignment -- the object of Alignment class which represent PSL-line:
    # get alignment from PSL-line means read line in PSL-file and get blocks's and aligned transcript's coordinates:
    if is_psl_lines:
        current_alignment1 = Alignment.PSLFileAlignment.get_alignment_from_psl_line(line1)
    # get alignment -- the object of Alignment class which represent Blast6-line:
    else:
        current_alignment1 = Alignment.BLAST6FileAlignment.get_alignment_from_blast6_line(line1)

    #logger.debug('    Processing {}...'.format(current_alignment1.query_fragment.name))

    if current_alignment1.query_fragment.end - current_alignment1.query_fragment.start >= min_alignment_threshold:
        single_transcript_lines.append(line1)
        single_transcript_alignments.append(current_alignment1)

    if line2 != '':
        if is_psl_lines:
            current_alignment2 = Alignment.PSLFileAlignment.get_alignment_from_psl_line(line2)
        else:
            current_alignment2 = Alignment.BLAST6FileAlignment.get_alignment_from_blast6_line(line2)

        while current_alignment1.query_fragment.name == current_alignment2.query_fragment.name and line2 != '':
            if current_alignment2.query_fragment.end - current_alignment2.query_fragment.start >= min_alignment_threshold:
                single_transcript_lines.append(line2)
                single_transcript_alignments.append(current_alignment2)

            current_alignment1 = current_alignment2
            line2 = fin.readline().strip()
            if line2 != '':
                if is_psl_lines:
                    current_alignment2 = Alignment.PSLFileAlignment.get_alignment_from_psl_line(line2)
                else:
                    current_alignment2 = Alignment.BLAST6FileAlignment.get_alignment_from_blast6_line(line2)

    line1 = line2
    line2 = fin.readline().strip()

    return single_transcript_lines, single_transcript_alignments, line1, line2


# temporary maybe needs some processing this alignments:
def remove_strange_psl_alignments(in_psl_path, out_psl_path):
    in_handle = open(in_psl_path)

    out_handle = open(out_psl_path, 'w')

    for line in in_handle:
        psl_alignment = Alignment.PSLFileAlignment.get_alignment_from_psl_line(line)
        strange = False
        for i_block in range(psl_alignment.blocks_num - 1):
            if psl_alignment.query_fragment.ends[i_block] - psl_alignment.query_fragment.starts[i_block + 1] > 0:
                strange = True
                break
        if strange == False:
            out_handle.write(line)

    in_handle.close()

    out_handle.close()

    return out_psl_path