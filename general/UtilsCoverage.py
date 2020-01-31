__author__ = 'lenk'

from general import UtilsGeneral
from general import UtilsAnnotations


# set count of covered bases of exons/blocks by covering them blocks/exons:
def get_coverage_positions(target_ids, target_starts, target_ends, query_ids, query_starts, query_ends):
    target_cov_pos = {}
    query_cov_pos = {}

    tmp_stack_i_exons = []
    tmp_stack_i_blocks = []

    # needs sorted coordinates:
    dict_coordinates = {}
    dict_coordinates['blocks_starts'] = query_starts
    dict_coordinates['blocks_ends'] = query_ends

    dict_coordinates['exons_starts'] = target_starts
    dict_coordinates['exons_ends'] = target_ends

    i_current = {'blocks_starts': 0, 'blocks_ends': 0, 'exons_starts': 0, 'exons_ends': 0}

    # get count of covered bases:
    while dict_coordinates != {}:
        (current, key) = min((dict_coordinates[key0][i_current[key0]], key0) for key0 in dict_coordinates.keys())

        # if start and end are equal, we choose start:
        if key == 'blocks_ends' and 'blocks_starts' in dict_coordinates and \
                        dict_coordinates['blocks_ends'][i_current['blocks_ends']] == dict_coordinates['blocks_starts'][i_current['blocks_starts']]:
                key = 'blocks_starts'
        elif key == 'blocks_ends' and 'exons_starts' in dict_coordinates and \
                        dict_coordinates['exons_starts'][i_current['exons_starts']] == dict_coordinates['blocks_ends'][i_current['blocks_ends']]:
                key = 'exons_starts'
        elif key == 'exons_ends' and 'exons_starts' in dict_coordinates and \
                        dict_coordinates['exons_ends'][i_current['exons_ends']] == dict_coordinates['exons_starts'][i_current['exons_starts']]:
                key = 'exons_starts'
        elif key == 'exons_ends' and 'blocks_starts' in dict_coordinates and \
                        dict_coordinates['exons_ends'][i_current['exons_ends']] == dict_coordinates['blocks_starts'][i_current['blocks_starts']]:
                key = 'blocks_starts'

        if key == 'exons_starts':
            tmp_stack_i_exons.append(i_current[key])

        elif key == 'blocks_starts':
            tmp_stack_i_blocks.append(i_current[key])

        elif key == 'exons_ends':
            for i_block in tmp_stack_i_blocks:
                start_coverage = max(query_starts[i_block], target_starts[i_current[key]])
                end_coverage = min(query_ends[i_block], target_ends[i_current[key]])

                if target_ids[i_current[key]] not in target_cov_pos:
                    target_cov_pos[target_ids[i_current[key]]] = []
                target_cov_pos[target_ids[i_current[key]]].append((start_coverage - target_starts[i_current[key]],
                                                                   end_coverage - target_starts[i_current[key]]))

                if query_ids[i_block] not in query_cov_pos:
                    query_cov_pos[query_ids[i_block]] = []
                query_cov_pos[query_ids[i_block]].append((start_coverage - query_starts[query_ids[i_block]],
                                                          end_coverage - query_starts[query_ids[i_block]]))
            if i_current[key] in tmp_stack_i_exons:
                tmp_stack_i_exons.remove(i_current[key])

        elif key == 'blocks_ends':
            for i_exon in tmp_stack_i_exons:
                start_coverage = max(query_starts[i_current[key]], target_starts[i_exon])
                end_coverage = min(query_ends[i_current[key]], target_ends[i_exon])

                if target_ids[i_exon] not in target_cov_pos:
                    target_cov_pos[target_ids[i_exon]] = []
                target_cov_pos[target_ids[i_exon]].append((start_coverage - target_starts[i_exon],
                                                           end_coverage - target_starts[i_exon]))

                if query_ids[i_current[key]] not in query_cov_pos:
                    query_cov_pos[query_ids[i_current[key]]] = []
                query_cov_pos[query_ids[i_current[key]]].append((start_coverage - query_starts[query_ids[i_current[key]]],
                                                                 end_coverage - query_starts[query_ids[i_current[key]]]))
            if i_current[key] in tmp_stack_i_blocks:
                tmp_stack_i_blocks.remove(i_current[key])
        if i_current[key] == len(dict_coordinates[key]) - 1:
            del dict_coordinates[key]
        else:
            i_current[key] += 1

    return target_cov_pos, query_cov_pos


# def get_starts_ends_coverage(start_coverage, end_coverage, starts_coverage, ends_coverage):
#     positions = []
#     for i in range(len(starts_coverage)):
#         positions.append(starts_coverage[i])
#         positions.append(ends_coverage[i])
#
#     i_start = UtilsGeneral.get_bin_search_position_of_element(positions, start_coverage)
#     i_end = UtilsGeneral.get_bin_search_position_of_element(positions, end_coverage)
#
#     new_positions = []
#     if i_start == len(positions):
#         new_positions = positions[:] + [start_coverage]
#     else:
#         if start_coverage == positions[i_start] and i_start % 2 == 0:
#             new_positions = positions[:i_start + 1]
#         if start_coverage == positions[i_start] and i_start % 2 == 1:
#             new_positions = positions[:i_start]
#         if start_coverage != positions[i_start] and i_start % 2 == 0:
#             new_positions = positions[:i_start] + [start_coverage]
#         if start_coverage != positions[i_start] and i_start % 2 == 1:
#             new_positions = positions[:i_start]
#
#     if i_end == len(positions):
#         new_positions += [end_coverage]
#     else:
#         if end_coverage == positions[i_end] and i_end % 2 == 0:
#             new_positions += positions[i_end + 1:]
#         if end_coverage == positions[i_end] and i_end % 2 == 1:
#             new_positions += positions[i_end:]
#         if end_coverage != positions[i_end] and i_end % 2 == 0:
#             new_positions += [end_coverage] + positions[i_end:]
#         if end_coverage != positions[i_end] and i_end % 2 == 1:
#             new_positions += positions[i_end:]
#
#     new_starts_coverage = []
#     new_ends_coverage = []
#     for i in range(len(new_positions) / 2):
#         new_starts_coverage.append(new_positions[2 * i])
#         new_ends_coverage.append(new_positions[2 * i + 1])
#
#     return new_starts_coverage, new_ends_coverage


def get_internal_exons(strand, sqlite3_db_genes, target_name, target_starts, target_ends):
    internal_exons = set()

    for i_block in range(len(target_starts)):
        region = sqlite3_db_genes.region(seqid=target_name, featuretype=UtilsAnnotations.type_exons,
                                 start=target_starts[i_block],
                                 end=target_ends[i_block],
                                 strand=strand, completely_within=False)
        internal_exons.update(list(region))

    return internal_exons


def get_internal_isoforms(sqlite3_db_genes, type_isoforms, internal_exons):
    internal_isoforms = set()
    for exon in internal_exons:
        if exon.featuretype in UtilsAnnotations.default_type_exons:
            internal_isoforms.update(list(sqlite3_db_genes.parents(exon.id, featuretype=type_isoforms)))
        # for prokaryotes:
        else:
            internal_isoforms.update([exon])

    return internal_isoforms


# get internal exons faster (without gff utils binning strategy -- usage sorted coordinates of exons):
def get_internal_exons_faster(sqlite3_db_genes, sorted_exons_attr, alignment_t_starts, alignment_t_ends, strand, id_chr):
    ids_internal_exons = set()
    internal_exons = set()

    # ids_internal_exons_tmp = set()

    for i_block in range(len(alignment_t_starts)):
        bin_start_i_in_ends, bin_end_i_in_ends = \
            get_bin_indexes(alignment_t_starts[i_block], sorted_exons_attr.sort_target_ends[str(strand)][id_chr],
                            sorted_exons_attr.index_sort_ends[str(strand)][id_chr], sorted_exons_attr.index_step[id_chr])
        if bin_start_i_in_ends is not None and bin_end_i_in_ends is not None:
            bin_ends = sorted_exons_attr.sort_target_ends[str(strand)][id_chr][bin_start_i_in_ends:bin_end_i_in_ends]
            begin = UtilsGeneral.get_bin_search_position_of_element(bin_ends, alignment_t_starts[i_block]) + bin_start_i_in_ends
            ids_ends_set = set(sorted_exons_attr.ids_by_end[strand][id_chr][begin:])
        else:
            ids_ends_set = set()


        bin_start_i_in_starts, bin_end_i_in_starts = \
            get_bin_indexes(alignment_t_ends[i_block], sorted_exons_attr.sort_target_starts[str(strand)][id_chr],
                            sorted_exons_attr.index_sort_starts[str(strand)][id_chr], sorted_exons_attr.index_step[id_chr])
        if bin_start_i_in_starts is not None and bin_end_i_in_starts is not None:
            bin_starts = sorted_exons_attr.sort_target_starts[str(strand)][id_chr][bin_start_i_in_starts:bin_end_i_in_starts]
            end = UtilsGeneral.get_bin_search_position_of_element(bin_starts, alignment_t_ends[i_block]) + bin_start_i_in_starts
            while end != len(sorted_exons_attr.sort_target_starts[strand][id_chr]) and end == sorted_exons_attr.sort_target_starts[strand][id_chr][end]:
                end += 1
            ids_starts_set = set(sorted_exons_attr.ids_by_start[strand][id_chr][:end])
        else:
            ids_starts_set = set(sorted_exons_attr.ids_by_start[strand][id_chr])

        ids_internal_exons = ids_internal_exons.union(ids_ends_set.intersection(ids_starts_set))

        for id_exon in ids_internal_exons:
            internal_exons.add(sqlite3_db_genes[id_exon])

        # TODO: delete
        # begin_tmp = UtilsGeneral.get_bin_search_position_of_element(sorted_exons_attr.sort_target_ends[str(strand)][id_chr], alignment_t_starts[i_block])
        # ids_ends_set_tmp = set(sorted_exons_attr.ids_by_end[strand][id_chr][begin_tmp:])
        #
        #
        # end_tmp = UtilsGeneral.get_bin_search_position_of_element(sorted_exons_attr.sort_target_starts[str(strand)][id_chr], alignment_t_ends[i_block])
        # while end_tmp != len(sorted_exons_attr.sort_target_starts[strand][id_chr]) and end_tmp == sorted_exons_attr.sort_target_starts[strand][id_chr][end_tmp]:
        #     end_tmp += 1
        # ids_starts_set_tmp = set(sorted_exons_attr.ids_by_start[strand][id_chr][:end_tmp])
        #
        # ids_internal_exons_tmp = ids_internal_exons_tmp.union(ids_ends_set_tmp.intersection(ids_starts_set_tmp))
        #
        # if ids_internal_exons != ids_internal_exons_tmp:
        #     print '!!!!'
        #     import sys
        #     sys.exit()

    return list(internal_exons)


def get_bin_indexes(alignment_pos, arr, index_arr, step_i):
    int_part_start = alignment_pos // step_i

    if int_part_start >= len(index_arr):
        return None, None
    start_index = index_arr[int_part_start]
    end_index = max(len(arr), start_index + step_i)

    return start_index, end_index


# get id best mapped isoform or transcript:
# choose between internal isoforms for mapped transcript or between transcripts mapped to isoforms:
def get_ids_best_mapped(transcripts_covered_bases, isoforms_covered_fraction):
    # choose isoforms ids that correspond isoforms which covered whole transcript in the best way:
    ids_isoforms_max, max_transcript_coverage = get_keys_corr_max_value(transcripts_covered_bases)

    # form isoforms coverage dictionary obtained from isoforms ids corresponded max whole transcript coverage:
    covered_isoform_fraction_max = {id_isoform: isoforms_covered_fraction[id_isoform] for id_isoform in ids_isoforms_max}

    # choose max isoform coverage over previously obtained dictionary:
    ids_max, max = get_keys_corr_max_value(covered_isoform_fraction_max)

    return ids_max


def get_keys_corr_max_value(dictionary):
    max_value = None
    if len(dictionary.keys()) != 0:
        max_value = max(dictionary.values())

    keys_corr_max = []
    for key in dictionary:
        if dictionary[key] == max_value:
            keys_corr_max.append(key)

    return keys_corr_max, max_value