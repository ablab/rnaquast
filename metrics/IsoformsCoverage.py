__author__ = 'lenk'

from general import UtilsAnnotations


class IsoformsCoverage():
    """Class, which represent coverage of annotated isoforms by aligned exons"""

    def __init__(self):
        # COVERAGE BY TRANSCRIPTS:
        # CONSIDER COVERED BASES (BY ALL MAPPED TRANSCRIPTS) WITHOUT OVERLAPS:
        # GENES:
        self.ids_well_covered_genes = set()
        self.num_well_covered_genes = 0

        self.ids_fully_covered_genes = set()
        self.num_fully_covered_genes = 0

        # ISOFORMS:
        self.avg_covered_fraction = 0.0
        self.covered_fraction = {}
        self.covered_fraction_distribution = {}

        self.ids_well_covered_isoforms = set()
        self.num_well_covered_isoforms = 0

        self.ids_fully_covered_isoforms = set()
        self.num_fully_covered_isoforms = 0

        # EXONS:
        self.avg_covered_fraction_exons = 0.0
        self.covered_fraction_exons = {}
        self.covered_fraction_exons_distribution = {}

        self.ids_well_covered_exons = set()
        self.num_well_covered_exons = 0

        self.ids_fully_covered_exons = set()
        self.num_fully_covered_exons = 0


        # number and percentage of well-covered exons for each isoform:
        self.num_isoform_well_covered_exons = {}
        self.percentage_isoform_well_covered_exons = {}
        # number and percentage of fully-covered exons for each isoform:
        self.num_isoform_fully_covered_exons = {}
        self.percentage_isoform_fully_covered_exons = {}

        # average over annotated isoforms percentage of well-covered exons and distribution:
        self.avg_percentage_isoform_well_covered_exons = 0.0
        self.percentage_isoform_well_covered_exons_distribution = {}
        # average over annotated isoforms percentage of fully-covered exons and distribution:
        self.avg_percentage_isoform_fully_covered_exons = 0.0
        self.percentage_isoform_fully_covered_exons_distribution = {}


        # CONSIDER COVERED BASES BY EACH MAPPED TRANSCRIPT SEPARATELY:
        # GENES:
        self.num_assembled_genes = 0
        self.ids_assembled_genes = set()

        self.ids_well_assembled_genes = set()
        self.num_well_assembled_genes = 0

        self.ids_fully_assembled_genes = set()
        self.num_fully_assembled_genes = 0

        # ISOFORMS:
        # number of annotated isoforms mapped to at least one transcript:
        self.num_assembled_isoforms = 0
        self.ids_assembled_isoforms = set()

        self.ids_well_assembled_isoforms = set()
        self.num_well_assembled_isoforms = 0

        self.ids_fully_assembled_isoforms = set()
        self.num_fully_assembled_isoforms = 0

        # average over annotated isoforms percentage of assembled bases:
        self.avg_assembled_fraction = 0.0
        # percentage and distribution of assembled bases by aligned transcript:
        self.assembled_fraction = {}
        self.assembled_fraction_distribution = {}

        # EXONS:
        self.num_assembled_exons = 0
        self.ids_assembled_exons = set()

        self.avg_assembled_fraction_exons = 0.0
        self.assembled_bases_exons = {}
        self.assembled_fraction_exons = {}
        self.assembled_fraction_exons_distribution = {}

        self.num_well_assembled_exons = 0
        self.ids_well_assembled_exons = set()

        self.num_fully_assembled_exons = 0
        self.ids_fully_assembled_exons = set()


        # dictionary of isoforms and ids of mapped to them transcripts:
        # self.ids_transcripts_mapped_to_isoform = {}
        # dictionary of mapped isoforms ids and number of transcripts mapped to annotated isoforms:
        self.num_transcripts_mapped_to_isoform = {}

        # number of aligned transcripts covered this position:
        self.num_transcripts_covered_pos = {}

        # number of transcripts mapped to at least one isoform:
        self.annotated_transcripts_num = 0

        # FOR ALL TRANSCRIPT ALIGNMENTS:
        # dictionary of number of transcripts mapped to isoform and number of isoforms with this number of mapped transcripts:
        self.num_transcripts_mapped_to_isoform_distribution = {}
        # average duplication ratio:
        self.avg_num_transcripts_mapped_to_isoform = 0.0

        # DATABASE COVERAGE:
        # percentages of covered bases of annotated isoforms without repeats of coverages by all aligned transcripts:
        self.fraction_annotation_mapped = 0.0

        # DUPLICATION RATIO OVER ISOFORMS:
        self.num_covered_pos_at_least_one = 0
        self.avg_duplication_ratio = 0.0

        self.relative_database_coverage = None


    class RelativeDatabaseCoverage():

        def __init__(self, reads_coverage, isoforms_coverage):
            self.database_coverage = 0.0

            # GENES:
            # CONSIDER COVERED BASES BY EACH MAPPED TRANSCRIPT SEPARATELY:
            self.percentage_well_assembled_genes = 0.0
            self.percentage_fully_assembled_genes = 0.0

            # CONSIDER COVERED BASES (BY ALL ALIGNMENTS) WITHOUT OVERLAPS:
            self.percentage_well_covered_genes = 0.0
            self.percentage_fully_covered_genes = 0.0


            # ISOFORMS:
            # CONSIDER COVERED BASES BY EACH ALIGNMENT SEPARATELY:
            self.percentage_well_assembled_isoforms = 0.0
            self.percentage_fully_assembled_isoforms = 0.0

            # CONSIDER COVERED BASES (BY ALL ALIGNMENTS) WITHOUT OVERLAPS:
            self.percentage_well_covered_isoforms = 0.0
            self.percentage_fully_covered_isoforms = 0.0


            # EXONS:
            # CONSIDER COVERED BASES BY EACH ALIGNMENT SEPARATELY:
            self.percentage_well_assembled_exons = 0.0
            self.percentage_fully_assembled_exons = 0.0

            # CONSIDER COVERED BASES (BY ALL ALIGNMENTS) WITHOUT OVERLAPS:
            self.percentage_well_covered_exons = 0.0
            self.percentage_fully_covered_exons = 0.0

            self.get_relative_database_coverage(reads_coverage, isoforms_coverage)


        def get_relative_database_coverage(self, reads_coverage, isoforms_coverage):
            if reads_coverage.fraction_annotation_mapped_by_reads != 0:
                self.database_coverage = isoforms_coverage.fraction_annotation_mapped / reads_coverage.fraction_annotation_mapped_by_reads

            # GENES:
            if reads_coverage.num_well_expressed_genes != 0:
                self.percentage_well_assembled_genes = isoforms_coverage.num_well_assembled_genes * 1.0 / reads_coverage.num_well_expressed_genes
                self.percentage_well_covered_genes = isoforms_coverage.num_well_covered_genes * 1.0 / reads_coverage.num_well_expressed_genes

            if reads_coverage.num_fully_expressed_genes != 0:
                self.percentage_fully_assembled_genes = isoforms_coverage.num_fully_assembled_genes * 1.0 / reads_coverage.num_fully_expressed_genes
                self.percentage_fully_covered_genes = isoforms_coverage.num_fully_covered_genes * 1.0 / reads_coverage.num_fully_expressed_genes

            # ISOFORMS:
            if reads_coverage.num_well_expressed_isoforms != 0:
                self.percentage_well_assembled_isoforms = isoforms_coverage.num_well_assembled_isoforms * 1.0 / reads_coverage.num_well_expressed_isoforms
                self.percentage_well_covered_isoforms = isoforms_coverage.num_well_covered_isoforms * 1.0 / reads_coverage.num_well_expressed_isoforms

            if reads_coverage.num_fully_expressed_isoforms != 0:
                self.percentage_fully_assembled_isoforms = isoforms_coverage.num_fully_assembled_isoforms * 1.0 / reads_coverage.num_fully_expressed_isoforms
                self.percentage_fully_covered_isoforms = isoforms_coverage.num_fully_covered_isoforms * 1.0 / reads_coverage.num_fully_expressed_isoforms

            # EXONS:
            if reads_coverage.num_well_expressed_exons != 0:
                self.percentage_well_assembled_exons = isoforms_coverage.num_well_assembled_exons * 1.0 / reads_coverage.num_well_expressed_exons
                self.percentage_well_covered_exons = isoforms_coverage.num_well_covered_exons * 1.0 / reads_coverage.num_well_expressed_exons

            if reads_coverage.num_fully_expressed_exons != 0:
                self.percentage_fully_assembled_exons = isoforms_coverage.num_fully_assembled_exons * 1.0 / reads_coverage.num_fully_expressed_exons
                self.percentage_fully_covered_exons = isoforms_coverage.num_fully_covered_exons * 1.0 / reads_coverage.num_fully_expressed_exons


    def update_isoforms_coverage_by_specific_isoform(self, sqlite3_db_genes, internal_isoforms_coverage, id_isoform):
        self.annotated_transcripts_num += 1

        if id_isoform not in self.num_transcripts_mapped_to_isoform:
            parent_genes = list(sqlite3_db_genes.parents(id_isoform, featuretype=UtilsAnnotations.default_type_genes))
            if parent_genes == []:
                parent_gene_id = id_isoform
            else:
                parent_gene_id = parent_genes[0].id

            self.num_transcripts_mapped_to_isoform[id_isoform] = 0

            self.ids_assembled_isoforms.add(id_isoform)

            self.ids_assembled_genes.add(parent_gene_id)

            self.assembled_fraction[id_isoform] = internal_isoforms_coverage.assembled_fraction[id_isoform]

            self.num_transcripts_covered_pos[id_isoform] = {}

        self.num_transcripts_mapped_to_isoform[id_isoform] += 1

        if internal_isoforms_coverage.assembled_fraction[id_isoform] > self.assembled_fraction[id_isoform]:
            self.assembled_fraction[id_isoform] = internal_isoforms_coverage.assembled_fraction[id_isoform]

        for id_exon in internal_isoforms_coverage.num_transcripts_covered_pos[id_isoform]:
            self.ids_assembled_exons.add(id_exon)

            len_exon = len(internal_isoforms_coverage.num_transcripts_covered_pos[id_isoform][id_exon])

            if id_exon not in self.num_transcripts_covered_pos[id_isoform]:
                self.num_transcripts_covered_pos[id_isoform][id_exon] = [0] * len_exon

            if id_exon not in self.assembled_bases_exons or \
                            internal_isoforms_coverage.assembled_bases_exons[id_isoform][id_exon] > self.assembled_bases_exons[id_exon]:
                self.assembled_bases_exons[id_exon] = internal_isoforms_coverage.assembled_bases_exons[id_isoform][id_exon]
                self.assembled_fraction_exons[id_exon] = internal_isoforms_coverage.assembled_fraction_exons[id_isoform][id_exon]

            self.num_transcripts_covered_pos[id_isoform][id_exon] = \
                map(lambda a, b: a + b, self.num_transcripts_covered_pos[id_isoform][id_exon],
                    internal_isoforms_coverage.num_transcripts_covered_pos[id_isoform][id_exon])


    def get_isoforms_coverage(self, sqlite3_db_genes, tot_isoforms_len, reads_coverage, WELL_FULLY_COVERAGE_THRESHOLDS):
        self.num_assembled_genes = len(self.ids_assembled_genes)

        self.num_assembled_isoforms = len(self.ids_assembled_isoforms)

        self.num_assembled_exons = len(self.ids_assembled_exons)

        self.covered_fraction_exons = {}
        for id_isoform in self.num_transcripts_mapped_to_isoform:
            isoform = sqlite3_db_genes[id_isoform]

            parent_genes = list(sqlite3_db_genes.parents(id_isoform, featuretype=UtilsAnnotations.default_type_genes))
            if parent_genes == []:
                parent_gene_id = id_isoform
            else:
                parent_gene_id = parent_genes[0].id

            children_exons = list(sqlite3_db_genes.children(isoform.id, featuretype=UtilsAnnotations.default_type_exons, order_by='start'))
            # for prokaryotes:
            if len(children_exons) == 0:
                children_exons = [isoform]

            len_isoform = 0
            for exon in children_exons:
                len_isoform += len(exon)

            self.num_isoform_well_covered_exons[id_isoform] = 0
            self.percentage_isoform_well_covered_exons[id_isoform] = 0.0

            self.num_isoform_fully_covered_exons[id_isoform] = 0
            self.percentage_isoform_fully_covered_exons[id_isoform] = 0.0

            if self.num_transcripts_mapped_to_isoform[id_isoform] not in self.num_transcripts_mapped_to_isoform_distribution:
                self.num_transcripts_mapped_to_isoform_distribution[self.num_transcripts_mapped_to_isoform[id_isoform]] = 0
            self.num_transcripts_mapped_to_isoform_distribution[self.num_transcripts_mapped_to_isoform[id_isoform]] += 1

            # CONSIDER COVERED BASES (BY ALL MAPPED TRANSCRIPTS) WITHOUT OVERLAPS:
            self.covered_fraction[id_isoform] = 0.0
            for id_exon in self.num_transcripts_covered_pos[id_isoform]:
                len_exon = len(self.num_transcripts_covered_pos[id_isoform][id_exon])
                num_uncovered_bases = self.num_transcripts_covered_pos[id_isoform][id_exon].count(0)
                num_covered_bases = len_exon - num_uncovered_bases
                covered_fraction_exon = num_covered_bases * 1.0 / len_exon

                if covered_fraction_exon >= WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold:
                    self.num_isoform_well_covered_exons[id_isoform] += 1

                if covered_fraction_exon >= WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold:
                    self.num_isoform_fully_covered_exons[id_isoform] += 1

                if id_exon not in self.covered_fraction_exons or covered_fraction_exon > self.covered_fraction_exons[id_exon]:
                    self.covered_fraction_exons[id_exon] = covered_fraction_exon
                self.covered_fraction[id_isoform] += num_covered_bases

                self.num_covered_pos_at_least_one += num_covered_bases
                self.avg_duplication_ratio += sum(self.num_transcripts_covered_pos[id_isoform][id_exon])

            self.covered_fraction[id_isoform] /= len_isoform
            if self.covered_fraction[id_isoform] not in self.covered_fraction_distribution:
                self.covered_fraction_distribution[self.covered_fraction[id_isoform]] = 0
            self.covered_fraction_distribution[self.covered_fraction[id_isoform]] += 1

            if len(children_exons) != 0:
                self.percentage_isoform_well_covered_exons[id_isoform] = \
                    self.num_isoform_well_covered_exons[id_isoform] * 1.0 / len(children_exons)
                self.percentage_isoform_fully_covered_exons[id_isoform] = \
                    self.num_isoform_fully_covered_exons[id_isoform] * 1.0 / len(children_exons)

            # isoform well-covered by transcripts if this isoform have more then well-threshold * sumLengthOfIsoform covered bases:
            if self.covered_fraction[id_isoform] >= WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold:
                self.ids_well_covered_isoforms.add(id_isoform)
                self.ids_well_covered_genes.add(parent_gene_id)

            # isoform fully-covered by transcripts if this isoform have more then fully-threshold * sumLengthOfIsoform covered bases:
            if self.covered_fraction[id_isoform] >= WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold:
                self.ids_fully_covered_isoforms.add(id_isoform)
                self.ids_fully_covered_genes.add(parent_gene_id)


            # CONSIDER COVERED BASES BY EACH MAPPED TRANSCRIPT SEPARATELY:
            # isoform well-covered by transcripts if this isoform have more then well-threshold * sumLengthOfIsoform covered bases:
            if self.assembled_fraction[id_isoform] >= WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold:
                self.ids_well_assembled_isoforms.add(id_isoform)
                self.ids_well_assembled_genes.add(parent_gene_id)

            # isoform fully-covered by transcripts if this isoform have more then fully-threshold * sumLengthOfIsoform covered bases:
            if self.assembled_fraction[id_isoform] >= WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold:
                self.ids_fully_assembled_isoforms.add(id_isoform)
                self.ids_fully_assembled_genes.add(parent_gene_id)

            # self.avg_assembled_bases += self.assembled_bases[id_isoform]
            self.avg_assembled_fraction += self.assembled_fraction[id_isoform]

            if self.assembled_fraction[id_isoform] not in self.assembled_fraction_distribution:
                self.assembled_fraction_distribution[self.assembled_fraction[id_isoform]] = 0
            self.assembled_fraction_distribution[self.assembled_fraction[id_isoform]] += 1

            # average over annotated isoforms number, percentage of well-covered exons and distribution:
            self.avg_percentage_isoform_well_covered_exons += self.percentage_isoform_well_covered_exons[id_isoform]
            if self.percentage_isoform_well_covered_exons[id_isoform] not in self.percentage_isoform_well_covered_exons_distribution:
                self.percentage_isoform_well_covered_exons_distribution[self.percentage_isoform_well_covered_exons[id_isoform]] = 0
            self.percentage_isoform_well_covered_exons_distribution[self.percentage_isoform_well_covered_exons[id_isoform]] += 1

            # average over annotated isoforms number, percentage of fully-covered exons and distribution:
            self.avg_percentage_isoform_fully_covered_exons += self.percentage_isoform_fully_covered_exons[id_isoform]
            if self.percentage_isoform_fully_covered_exons[id_isoform] not in self.percentage_isoform_fully_covered_exons_distribution:
                self.percentage_isoform_fully_covered_exons_distribution[self.percentage_isoform_fully_covered_exons[id_isoform]] = 0
            self.percentage_isoform_fully_covered_exons_distribution[self.percentage_isoform_fully_covered_exons[id_isoform]] += 1


        # CONSIDER COVERED BASES BY EACH MAPPED TRANSCRIPT SEPARATELY:
        self.num_well_assembled_isoforms = len(self.ids_well_assembled_isoforms)
        self.num_fully_assembled_isoforms = len(self.ids_fully_assembled_isoforms)
        self.num_well_assembled_genes = len(self.ids_well_assembled_genes)
        self.num_fully_assembled_genes = len(self.ids_fully_assembled_genes)

        # CONSIDER COVERED BASES (BY ALL MAPPED TRANSCRIPTS) WITHOUT OVERLAPS:
        self.num_well_covered_isoforms = len(self.ids_well_covered_isoforms)
        self.num_fully_covered_isoforms = len(self.ids_fully_covered_isoforms)
        self.num_well_covered_genes = len(self.ids_well_covered_genes)
        self.num_fully_covered_genes = len(self.ids_fully_covered_genes)

        if tot_isoforms_len != 0:
            self.fraction_annotation_mapped = self.num_covered_pos_at_least_one * 1.0 / tot_isoforms_len

        if self.num_covered_pos_at_least_one != 0:
            self.avg_duplication_ratio /= self.num_covered_pos_at_least_one

        for key in self.num_transcripts_mapped_to_isoform_distribution:
            self.avg_num_transcripts_mapped_to_isoform += key * self.num_transcripts_mapped_to_isoform_distribution[key]
        if sum(self.num_transcripts_mapped_to_isoform_distribution.values()) != 0:
            self.avg_num_transcripts_mapped_to_isoform /= sum(self.num_transcripts_mapped_to_isoform_distribution.values())


        # ISOFORMS:
        if self.num_assembled_isoforms != 0:
            # CONSIDER COVERED BASES BY EACH MAPPED TRANSCRIPT SEPARATELY:
            # self.avg_assembled_bases = self.avg_assembled_bases * 1.0 / self.num_assembled
            self.avg_assembled_fraction /= self.num_assembled_isoforms

            # CONSIDER COVERED BASES (BY ALL MAPPED TRANSCRIPTS) WITHOUT OVERLAPS:
            self.avg_covered_fraction = sum(self.covered_fraction.values()) / self.num_assembled_isoforms

            # average over annotated isoforms percentage of well-covered exons:
            self.avg_percentage_isoform_well_covered_exons /= self.num_assembled_isoforms
            # average over annotated isoforms percentage of fully-covered exons:
            self.avg_percentage_isoform_fully_covered_exons /= self.num_assembled_isoforms


        # EXONS:
        for id_exon in self.covered_fraction_exons:
            if self.covered_fraction_exons[id_exon] >= WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold:
                self.ids_well_covered_exons.add(id_exon)

            if self.covered_fraction_exons[id_exon] >= WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold:
                self.ids_fully_covered_exons.add(id_exon)

            if self.assembled_fraction_exons[id_exon] >= WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold:
                self.ids_well_assembled_exons.add(id_exon)

            if self.assembled_fraction_exons[id_exon] >= WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold:
                self.ids_fully_assembled_exons.add(id_exon)

            if self.covered_fraction_exons[id_exon] not in self.covered_fraction_exons_distribution:
                self.covered_fraction_exons_distribution[self.covered_fraction_exons[id_exon]] = 0
            self.covered_fraction_exons_distribution[self.covered_fraction_exons[id_exon]] += 1

            if self.assembled_fraction_exons[id_exon] not in self.assembled_fraction_exons_distribution:
                self.assembled_fraction_exons_distribution[self.assembled_fraction_exons[id_exon]] = 0
            self.assembled_fraction_exons_distribution[self.assembled_fraction_exons[id_exon]] += 1

        self.num_well_covered_exons = len(self.ids_well_covered_exons)

        self.num_fully_covered_exons = len(self.ids_fully_covered_exons)

        self.num_well_assembled_exons = len(self.ids_well_assembled_exons)

        self.num_fully_assembled_exons = len(self.ids_fully_assembled_exons)

        if self.num_assembled_exons != 0:
            self.avg_covered_fraction_exons = sum(self.covered_fraction_exons.values()) / self.num_assembled_exons

            self.avg_assembled_fraction_exons = sum(self.assembled_fraction_exons.values()) / self.num_assembled_exons


        # RELATIVE DATABASE COVERAGE:
        if reads_coverage is not None:
            self.relative_database_coverage = IsoformsCoverage.RelativeDatabaseCoverage(reads_coverage, self)


    def print_fully_assembled_isoforms(self, path_fully_assembled_list, logger):
        logger.info('    Getting Fully assembled isoforms list...')

        with open(path_fully_assembled_list, 'w') as fout:
            for id_isoform in self.ids_fully_assembled_isoforms:
                fout.write(id_isoform + '\n')

        logger.info('      saved to {}'.format(path_fully_assembled_list))