__author__ = 'lenk'


class InternalIsoformsCoverage(object):
    """Class, which represent coverage of annotated genes/internal for aligned transcript genes's(whom) by aligned exons(who)"""

    def __init__(self, internal_isoforms):
        # number of annotated isoforms in aligned transcript or in annotations:
        self.isoforms_num = len(internal_isoforms)

        # CONSIDER COVERED BASES BY ONE MAPPED TRANSCRIPT:
        # number and percentage of covered bases for each isoform:
        self.assembled_fraction = {}

        self.assembled_bases_exons = {}
        self.assembled_fraction_exons = {}

        self.num_transcripts_covered_pos = {}


    def update_internal_isoforms_coverage(self, sqlite3_db_genes, id_isoform, exon_cov_pos):
        if id_isoform not in self.assembled_fraction:
            self.assembled_fraction[id_isoform] = 0.0
            self.assembled_bases_exons[id_isoform] = {}

            self.num_transcripts_covered_pos[id_isoform] = {}

        for id_exon in exon_cov_pos:
            if id_exon not in self.assembled_bases_exons[id_isoform]:
                self.assembled_bases_exons[id_isoform][id_exon] = 0

                self.num_transcripts_covered_pos[id_isoform][id_exon] = [0] * len(sqlite3_db_genes[id_exon])

            for i_pos in range(len(exon_cov_pos[id_exon])):
                start_coverage = exon_cov_pos[id_exon][i_pos][0]
                end_coverage = exon_cov_pos[id_exon][i_pos][1]

                self.assembled_fraction[id_isoform] += end_coverage - start_coverage + 1


                self.assembled_bases_exons[id_isoform][id_exon] += end_coverage - start_coverage + 1

                for i_pos in range(start_coverage, end_coverage + 1):
                    self.num_transcripts_covered_pos[id_isoform][id_exon][i_pos] += 1


    def get_internal_isoforms_coverage(self, internal_isoforms, children_exons_dict):
        for isoform in internal_isoforms:
            id_isoform = isoform.id

            self.assembled_fraction_exons[id_isoform] = {}
            for id_exon in self.assembled_bases_exons[id_isoform]:
                exon_len = len(self.num_transcripts_covered_pos[id_isoform][id_exon])

                self.assembled_fraction_exons[id_isoform][id_exon] = \
                    self.assembled_bases_exons[id_isoform][id_exon] * 1.0 / exon_len

            len_isoform = 0
            for exon in children_exons_dict[id_isoform]:
                len_isoform += len(exon)
            self.assembled_fraction[id_isoform] /= len_isoform