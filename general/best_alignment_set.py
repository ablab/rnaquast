__author__ = 'letovesnoi'

from general import UtilsGeneral


# union alignments crossing at most err_cross in query and distant no more than err_space in query:
# for determine misassemblies find best union with score define by union_penalty for non close (define as fake blat) alignments:
def get_best_alignment_set(transcript_alignments, ALIGNMENT_THRESHOLDS):
    #logger.debug('      Getting best union alignments...')
    best_score = - float('Inf')
    best_b = None

    q_ends = []
    for i in range(len(transcript_alignments)):
        q_ends.append(transcript_alignments[i].query_fragment.end)
    q_ends_sort_index, q_ends_sort_array = UtilsGeneral.get_order_indexes_elements(q_ends)

    best_list = [[]]
    scores = {str([]): 0}

    for i in q_ends_sort_index:
        a_i = transcript_alignments[i]
        best_i, curr_best_score = \
            get_best_i(best_list, a_i, scores, ALIGNMENT_THRESHOLDS)

        best_list.append(best_i)

        scores[str(best_i)] = curr_best_score

        if curr_best_score >= best_score:
            best_b = best_i
            best_score = curr_best_score

    return best_b


# alignments are union if they are cross at most err_cross, distance no more than err_space in query, are at one strand and chromosome:
def is_union_fake_blat(a_0, a_1, ALIGNMENT_THRESHOLDS):
    query_cross = max(0, a_0.query_fragment.end - a_1.query_fragment.start + 1)

    query_dist = max(0, a_1.query_fragment.start - a_0.query_fragment.end - 1)

    if a_0.strand == '+':
        target_cross = max(0, a_0.target_fragment.end - a_1.target_fragment.start + 1)

        target_dist = max(0, a_1.target_fragment.start - a_0.target_fragment.end - 1)
    else:
        target_cross = max(0, a_1.target_fragment.end - a_0.target_fragment.start + 1)

        target_dist = max(0, a_0.target_fragment.start - a_1.target_fragment.end - 1)

    if (target_cross <= ALIGNMENT_THRESHOLDS.ERR_CROSS_TARGET_FAKE_BLAT) and (target_dist <= ALIGNMENT_THRESHOLDS.ERR_SPACE_TARGET_FAKE_BLAT) \
            and (a_0.strand == a_1.strand) and (a_0.target_fragment.name == a_1.target_fragment.name) \
            and (query_cross <= ALIGNMENT_THRESHOLDS.ERR_CROSS_QUERY_FAKE_BLAT) and (query_dist <= ALIGNMENT_THRESHOLDS.ERR_SPACE_QUERY_FAKE_BLAT):
        return True

    return False


# b -- list of alignments ordered by end crossing at most err_cross in query
# a_i -- one alignment
def get_score_b_a_i(b, a_i, scores, ALIGNMENT_THRESHOLDS):
    if not b:
        return a_i.score

    b_i = b[-1]

    dist = max(0, a_i.query_fragment.start - b_i.query_fragment.end - 1)
    cross = max(0, b_i.query_fragment.end - a_i.query_fragment.start + 1)

    if cross > ALIGNMENT_THRESHOLDS.ERR_CROSS_QUERY_UNION:
        return - float('Inf')

    current_union_penalty = ALIGNMENT_THRESHOLDS.UNION_PENALTY
    if is_union_fake_blat(b_i, a_i, ALIGNMENT_THRESHOLDS) and a_i.format == 'psl':
        current_union_penalty = 0

    score_b_a_i = scores[str(b)] + a_i.score - current_union_penalty - cross

    return score_b_a_i


def get_best_i(best_list, a_i, scores, ALIGNMENT_THRESHOLDS):
    best_score = - float('Inf')

    for b in best_list:
        score_b_a_i = get_score_b_a_i(b, a_i, scores, ALIGNMENT_THRESHOLDS)
        if score_b_a_i >= best_score:
            best_score = score_b_a_i
            best_i = b + [a_i]

    return best_i, best_score


# alignment0 = Alignment.Alignment()
# alignment0.create('+', 'tst', 5, 2, 6)
#
# alignment1 = Alignment.Alignment()
# alignment1.create('+', 'tst', 10, 4, 13)
#
# alignment2 = Alignment.Alignment()
# alignment2.create('+', 'tst',20, 3, 22)
#
# alignment3 = Alignment.Alignment()
# alignment3.create('+', 'tst', 30, 24, 53)
#
# alignment4 = Alignment.Alignment()
# alignment4.create('+', 'tst', 25, 10, 34)
#
# transcript_alignments = [alignment0, alignment1, alignment2, alignment3, alignment4]
#
# best_list, best_scores = get_union(transcript_alignments)
#
# print best_list
# print best_scores
#
# for b in best_list:
#     for b_i in b:
#         print b_i.score
#     print '\n'

# b = get_best_b(best_list, best_scores)
# for b_i in b:
#     print b_i.score