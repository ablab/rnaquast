__author__ = 'letovesnoi'

from general import UtilsAnnotations


class GeneDatabaseMetrics():
    """Class of basic gene database metrics"""

    def __init__(self, sqlite3_db_genes, logger):
        # number of genes in annotations:
        self.genes_num = 0
        # number of protein coding genes in annotations:
        self.protein_coding_genes_num = 0

        # number of isoforms in annotations:
        self.isoforms_num = 0
        # number of protein coding isoforms in annotations:
        self.protein_coding_isoforms_num = 0

        self.tot_exons_num = 0

        self.tot_introns_num = 0

        # isoform length without introns and isoforms number having this length over all isoforms:
        self.avg_isoform_len = 0.0
        self.tot_isoforms_len = 0
        self.isoforms_len_distribution = {}

        # exon length and exons number having this length over all exons:
        self.avg_exon_len = 0.0
        self.exons_len_distribution = {}

        # intron length:
        self.avg_intron_len = 0.0

        # exons number and isoforms number having this number over all isoforms:
        self.avg_exons_num = 0.0
        self.max_exons_num = 0
        self.exons_num_distribution = {}

        # isoform length with introns and isoforms number having this size over all isoforms:
        # self.len_w_introns_distribution = {}
        # self.min_len_w_introns = 0
        # self.max_len_w_introns = 0
        # self.avg_len_w_introns = 0.0
        # self.tot_len_w_introns = 0

        # self.min_len = 0
        # self.max_len = 0

        # self.min_exons_num = 0
        # self.min_exon_len = 0
        # self.max_exon_len = 0

        # introns number and isoforms number having this number over all isoforms:
        # self.introns_num_distribution = {}
        # self.min_introns_num = 0
        # self.max_introns_num = 0
        # self.avg_introns_num = 0.0

        # intron length and introns number having this length over all introns over all isoforms:
        # self.introns_len_distribution = {}
        # self.min_intron_len = 0
        self.max_intron_len = 0
        # self.tot_introns_len = 0

        logger.print_timestamp()
        logger.info('Getting GENE DATABASE metrics...')

        genes = list(sqlite3_db_genes.features_of_type(UtilsAnnotations.default_type_genes))

        isoforms = list(sqlite3_db_genes.features_of_type(UtilsAnnotations.default_type_isoforms))

        # exons = list(sqlite3_db_genes.features_of_type(type_exons))

        self.genes_num = len(genes)
        self.isoforms_num = len(isoforms)

        for gene in genes:
            if gene.source == 'protein_coding' or \
                    ('gene_biotype' in gene.attributes and gene.attributes['gene_biotype'][0] == 'protein_coding') \
                    or ('biotype' in gene.attributes and gene.attributes['biotype'][0] == 'protein_coding'):
                self.protein_coding_genes_num += 1

        for transcript in isoforms:
            if transcript.source == 'protein_coding' or \
                    ('transcript_biotype' in transcript.attributes and transcript.attributes['transcript_biotype'][0] == 'protein_coding') \
                    or ('biotype' in transcript.attributes and transcript.attributes['biotype'][0] == 'protein_coding'):
                self.protein_coding_isoforms_num += 1

            exons = list(sqlite3_db_genes.children(transcript.id, featuretype=UtilsAnnotations.default_type_exons, order_by='start'))

            self.tot_exons_num += len(list(exons))
            isoform_len = 0
            for exon in exons:
                self.avg_exon_len += len(exon)
                isoform_len += len(exon)
                if len(exon) not in self.exons_len_distribution:
                    self.exons_len_distribution[len(exon)] = 0
                self.exons_len_distribution[len(exon)] += 1

            self.tot_isoforms_len += isoform_len
            if isoform_len not in self.isoforms_len_distribution:
                self.isoforms_len_distribution[isoform_len] = 0
            self.isoforms_len_distribution[isoform_len] += 1

            if len(list(exons)) not in self.exons_num_distribution:
                self.exons_num_distribution[len(list(exons))] = 0
            self.exons_num_distribution[len(list(exons))] += 1

            introns = list(sqlite3_db_genes.interfeatures(exons))
            self.tot_introns_num += len(list(introns))
            for intron in introns:
                intron_len = len(intron)

                self.avg_intron_len += intron_len

                if intron_len > self.max_intron_len:
                    self.max_intron_len = intron_len

        if self.isoforms_num != 0:
            self.avg_isoform_len = float(self.tot_isoforms_len) / self.isoforms_num
            self.avg_exons_num = float(self.tot_exons_num) / self.isoforms_num

        self.max_exons_num = max(self.exons_num_distribution.keys())

        if self.tot_exons_num != 0:
            self.avg_exon_len /= self.tot_exons_num

        if self.tot_introns_num != 0:
            self.avg_intron_len = float(self.avg_intron_len) / self.tot_introns_num

        logger.info('Done.')