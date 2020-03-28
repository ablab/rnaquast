__author__ = 'letovesnoi'

from general import UtilsTools
from general import UtilsCoverage
from general import UtilsAnnotations

from objects import Alignment


class ReadsCoverage():
    """Class, which represent coverage of annotated isoforms by aligned exons"""

    def __init__(self, sorted_sam_path, reference_path, single_reads, left_reads, right_reads,
                 reference_dict, sqlite3_db_genes, type_isoforms, sorted_exons_attr, strand_specific, tot_isoforms_len, genome_len,
                 output_dir, threads, WELL_FULLY_COVERAGE_THRESHOLDS, logger, log_dir):
        # COVERAGE BY READS (upper bound):
        # GENES:
        self.ids_well_expressed_genes = set()
        self.num_well_expressed_genes = 0

        self.ids_fully_expressed_genes = set()
        self.num_fully_expressed_genes = 0


        # ISOFORMS:
        self.ids_well_expressed_isoforms = set()
        self.num_well_expressed_isoforms = 0

        self.ids_fully_expressed_isoforms = set()
        self.num_fully_expressed_isoforms = 0


        # EXONS:
        self.ids_well_expressed_exons = set()
        self.num_well_expressed_exons = 0

        self.ids_fully_expressed_exons = set()
        self.num_fully_expressed_exons = 0


        # DATABASE COVERAGE:
        self.num_reads_covered_pos = {}
        self.num_expressed_pos_at_least_one_by_reads = 0
        self.fraction_annotation_mapped_by_reads = 0.0

        self.expressed_fraction_isoform = {}
        self.expressed_fraction_exon = {}

        self.get_database_coverage_by_reads(sorted_sam_path, reference_path, single_reads, left_reads,
                                            right_reads, reference_dict, sqlite3_db_genes, type_isoforms, sorted_exons_attr,
                                            strand_specific, tot_isoforms_len, genome_len, output_dir, threads,
                                            WELL_FULLY_COVERAGE_THRESHOLDS, logger, log_dir)


    def get_database_coverage_by_reads(self, sam_path, reference_path, single_reads, left_reads,
                                       right_reads, reference_dict, sqlite3_db_genes, type_isoforms, sorted_exons_attr,
                                       strand_specific, tot_isoforms_len, genome_len, output_dir, threads,
                                       WELL_FULLY_COVERAGE_THRESHOLDS, logger, log_dir):
            if sam_path is None:
                sam_path = \
                    UtilsTools.get_sam_by_STAR(threads, reference_path, None, single_reads, left_reads, right_reads,
                                                output_dir, None, None, genome_len, logger, log_dir)

            if sam_path is None:
                return

            logger.print_timestamp()
            logger.info('Getting database coverage by reads...')

            with open(sam_path, 'r') as in_handle:
                for line in in_handle:
                    if line[0] == '@':
                        continue

                    curr_sam_alignment = Alignment.SAMFileAlignment.get_alignment_from_sam_line(line, logger, reference_dict)

                    if strand_specific:
                        strand = curr_sam_alignment.strand
                    else:
                        strand = None

                    internal_exons = \
                        UtilsCoverage.get_internal_exons_faster(sqlite3_db_genes, sorted_exons_attr, curr_sam_alignment.target_fragment.starts,
                                                                    curr_sam_alignment.target_fragment.ends, str(strand), curr_sam_alignment.target_fragment.name)

                    internal_isoforms = list(UtilsCoverage.get_internal_isoforms(sqlite3_db_genes, type_isoforms, internal_exons))

                    for isoform in internal_isoforms:
                        children_exons = list(sqlite3_db_genes.children(isoform.id, featuretype=UtilsAnnotations.default_type_exons, order_by='start'))
                        # for prokaryotes:
                        if len(children_exons) == 0:
                            children_exons = [isoform]

                        exon_starts = [exon.start for exon in children_exons]
                        exon_ends = [exon.end for exon in children_exons]
                        exon_ids = [exon.id for exon in children_exons]

                        target_cov_pos, query_cov_pos = \
                            UtilsCoverage.get_coverage_positions(exon_ids, exon_starts, exon_ends, range(len(curr_sam_alignment.target_fragment.starts)),
                                                                 curr_sam_alignment.target_fragment.starts, curr_sam_alignment.target_fragment.ends)

                        if isoform.id not in self.num_reads_covered_pos:
                            self.num_reads_covered_pos[isoform.id] = {}
                        for id_exon in target_cov_pos:
                            if id_exon not in self.num_reads_covered_pos[isoform.id]:
                                self.num_reads_covered_pos[isoform.id][id_exon] = [0] * len(sqlite3_db_genes[id_exon])
                            for i_cov_pos in range(len(target_cov_pos[id_exon])):
                                for i_pos in range(target_cov_pos[id_exon][i_cov_pos][0], target_cov_pos[id_exon][i_cov_pos][1] + 1):
                                    self.num_reads_covered_pos[isoform.id][id_exon][i_pos] += 1

            for id_isoform in self.num_reads_covered_pos:
                parent_genes = list(sqlite3_db_genes.parents(id_isoform, featuretype=UtilsAnnotations.default_type_genes))
                if parent_genes == []:
                    parent_gene_id = id_isoform
                else:
                    parent_gene_id = parent_genes[0].id

                self.expressed_fraction_isoform[id_isoform] = 0.0
                len_isoform = 0
                for id_exon in self.num_reads_covered_pos[id_isoform]:
                    len_exon = len(self.num_reads_covered_pos[id_isoform][id_exon])
                    num_uncovered_pos = self.num_reads_covered_pos[id_isoform][id_exon].count(0)
                    num_expressed_bases = len_exon - num_uncovered_pos

                    if id_exon not in self.expressed_fraction_exon:
                        self.expressed_fraction_exon[id_exon] = num_expressed_bases * 1.0 / len_exon

                        if self.expressed_fraction_exon[id_exon] >= WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold:
                            self.ids_well_expressed_exons.add(id_exon)

                        if self.expressed_fraction_exon[id_exon] >= WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold:
                            self.ids_fully_expressed_exons.add(id_exon)

                    len_isoform += len_exon

                    self.num_expressed_pos_at_least_one_by_reads += num_expressed_bases

                    self.expressed_fraction_isoform[id_isoform] += num_expressed_bases

                self.expressed_fraction_isoform[id_isoform] /= len_isoform

                if self.expressed_fraction_isoform[id_isoform] >= WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold:
                    self.ids_well_expressed_genes.add(parent_gene_id)

                    self.ids_well_expressed_isoforms.add(id_isoform)

                if self.expressed_fraction_isoform[id_isoform] >= WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold:
                    self.ids_fully_expressed_genes.add(parent_gene_id)

                    self.ids_fully_expressed_isoforms.add(id_isoform)

            self.num_well_expressed_genes = len(self.ids_well_expressed_genes)
            self.num_well_expressed_isoforms = len(self.ids_well_expressed_isoforms)
            self.num_well_expressed_exons = len(self.ids_well_expressed_exons)

            self.num_fully_expressed_genes = len(self.ids_fully_expressed_genes)
            self.num_fully_expressed_isoforms = len(self.ids_fully_expressed_isoforms)
            self.num_fully_expressed_exons = len(self.ids_fully_expressed_exons)

            if tot_isoforms_len != 0:
                self.fraction_annotation_mapped_by_reads = self.num_expressed_pos_at_least_one_by_reads * 1.0 / tot_isoforms_len

            logger.info('Done.')

    def print_fully_expressed_isoforms(self, path_fully_expressed_list, logger):
        logger.info('    Getting Fully covered isoforms list by reads...')

        with open(path_fully_expressed_list, 'w') as fout:
            for id_isoform in self.ids_fully_expressed_isoforms:
                fout.write(id_isoform + '\n')

        logger.info('      saved to {}'.format(path_fully_expressed_list))


    def print_well_expressed_isoforms(self, path_well_expressed_list, logger):
        logger.info('    Getting Well covered isoforms list by reads...')

        with open(path_well_expressed_list, 'w') as fout:
            for id_isoform in self.ids_well_expressed_isoforms:
                fout.write(id_isoform + '\n')

        logger.info('      saved to {}'.format(path_well_expressed_list))