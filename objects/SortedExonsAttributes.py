__author__ = 'letovesnoi'

from general import UtilsAnnotations
from general import rqconfig

class SortedExonsAttributes():

    def __init__(self, sqlite3_db_genes, strands, ids_chrs, reference_dict, logger):
        # from datetime import datetime
        # start_time = datetime.now()

        logger.print_timestamp()
        logger.info('Sorting exons attributes...')

        if len(ids_chrs) > 100:
            logger.info()
            logger.warning('Number of chromosomes / scaffolds more than 100.')


        if len(list(sqlite3_db_genes.features_of_type(UtilsAnnotations.type_exons))) != 0:
            type = UtilsAnnotations.type_exons
        # for prokaryotes or bad annotated species:
        elif len(list(sqlite3_db_genes.features_of_type(UtilsAnnotations.type_isoforms))) != 0:
            type = UtilsAnnotations.type_isoforms
        elif len(list(sqlite3_db_genes.features_of_type(UtilsAnnotations.type_genes))) != 0:
            type = UtilsAnnotations.type_genes
        else:
            logger.error('Annotated exons not founded.', exit_with_code=2, to_stderr=True)

        self.ids_by_start = {}
        self.sort_target_starts = {}
        self.index_sort_starts = {}

        self.ids_by_end = {}
        self.sort_target_ends = {}
        self.index_sort_ends = {}

        self.index_step = {}

        for strand in strands:
            self.ids_by_start[str(strand)] = {}
            self.sort_target_starts[str(strand)] = {}
            self.index_sort_starts[str(strand)] = {}

            self.ids_by_end[str(strand)] = {}
            self.sort_target_ends[str(strand)] = {}
            self.index_sort_ends[str(strand)] = {}

            for id_chr in ids_chrs:
                self.index_sort_starts[str(strand)][id_chr] = []
                self.index_sort_ends[str(strand)][id_chr] = []

                self.index_step[id_chr] = len(reference_dict[id_chr]) / rqconfig.SORT_INDEX_LEN

                exons_by_start = list(sqlite3_db_genes.features_of_type(type, order_by='start', strand=strand, limit=(id_chr, 0, len(reference_dict[id_chr]) - 1)))
                self.sort_target_starts[str(strand)][id_chr] = [exon.start - 1 for exon in exons_by_start]
                self.ids_by_start[str(strand)][id_chr] = [exon.id for exon in exons_by_start]

                exons_by_end = list(sqlite3_db_genes.features_of_type(type, order_by='end', strand=strand, limit=(id_chr, 0, len(reference_dict[id_chr]) - 1)))
                self.sort_target_ends[str(strand)][id_chr] = [exon.end - 1 for exon in exons_by_end]
                self.ids_by_end[str(strand)][id_chr] = [exon.id for exon in exons_by_end]

                set_index_array(self.sort_target_starts[str(strand)][id_chr], self.index_sort_starts[str(strand)][id_chr], self.index_step[id_chr])

                set_index_array(self.sort_target_ends[str(strand)][id_chr], self.index_sort_ends[str(strand)][id_chr], self.index_step[id_chr])

                logger.info('  Sorted in {}.'.format(id_chr))

        logger.info('Done.')
        # elapsed_time = datetime.now() - start_time
        # print elapsed_time
        # import sys
        # sys.exit()


def set_index_array(arr, index_arr, step_i):
    curr_value = 0
    for i in range(len(arr)):
        while arr[i] >= curr_value:
            index_arr.append(i)
            curr_value += step_i
