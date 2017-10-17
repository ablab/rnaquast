__author__ = 'lenk'

from general import UtilsAnnotations
from general import UtilsCoverage


class AlignedTranscript(object):
    """class of aligned gene (aligned transcript), which represent line in PSL-file with alignments"""

    def __init__(self, psl_alignment, sorted_exons_attr, strand_specific, sqlite3_db_genes):
        # getting aligned transcript:
        self.alignment = psl_alignment

        # size of aligned transcript in genome reference:
        self.reference_len = self.alignment.target_fragment.end - self.alignment.target_fragment.start + 1

        self.alignment_len = sum(psl_alignment.blocks_sizes)
        self.fraction = self.alignment_len * 1.0 / self.alignment.query_fragment.size

        self.internal_isoforms = None
        self.ids_internal_isoforms = set()

        self.children_exons_dict = {}
        self.ids_children_exons_dict = {}

        if sqlite3_db_genes is not None:
            # internal_exons = self.get_internal_exons(db_genes, strand_specific)

            if strand_specific:
                strand = self.alignment.strand
            else:
                strand = str(None)
            internal_exons = \
                UtilsCoverage.get_internal_exons_faster(sqlite3_db_genes, sorted_exons_attr, self.alignment.target_fragment.starts,
                                                        self.alignment.target_fragment.ends, strand, self.alignment.target_fragment.name)
            self.internal_isoforms = list(UtilsCoverage.get_internal_isoforms(sqlite3_db_genes, internal_exons))
            for internal_isoform in self.internal_isoforms:
                self.ids_internal_isoforms.add(internal_isoform.id)

                self.children_exons_dict[internal_isoform.id] = list(sqlite3_db_genes.children(internal_isoform.id, featuretype=UtilsAnnotations.default_type_exons, order_by='start'))
                # for prokaryotes:
                if len(self.children_exons_dict[internal_isoform.id]) == 0:
                    self.children_exons_dict[internal_isoform.id] = [internal_isoform]

                self.ids_children_exons_dict[internal_isoform.id] = set([exon.id for exon in self.children_exons_dict[internal_isoform.id]])