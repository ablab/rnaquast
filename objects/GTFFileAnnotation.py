__author__ = 'lenk'


class GTFFileAnnotation(object):
    """class of GTF line, which represent annotation"""
    def __init__(self):
        # The name of the sequence. Must be a chromosome or scaffold
        self.seqname = ''
        # The program that generated this feature
        self.source = ''
        # The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon",
        # "stop_codon", and "exon"
        self.feature = ''
        # The starting position of the feature in the sequence. The first base is numbered 1
        self.start = 0
        # The ending position of the feature (inclusive)
        self.end = 0
        # A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set,
        # the score value will determine the level of gray in which this feature is displayed (higher numbers = darker
        # gray). If there is no score value, enter "."
        self.score = 0
        # Valid entries include '+', '-', or '.' (for don't know/don't care)
        self.strand = '.'
        # If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of
        # the first base. If the feature is not a coding exon, the value should be '.'
        self.frame = '.'
        # This list of attributes contain two mandatory attributes:
        # gene_id_value - a globally unique identifier for the genomic source of the sequence and
        # transcript_id_value - a globally unique identifier for the predicted transcript
        # may also contain exon_number and etc.
        self.attributes_line = ''
        self.attributes_list = {}

    # get attributes of class annotations from line in GTF-file:
    @classmethod
    def get_annotation_from_gtf_gff_file(cls, text_line):
        gtf_annotation = cls()
        # for GTF-file:
        if text_line[-1] == ';':
            text_line = text_line[:-1]

        parameters_list = text_line[:].split('\t')
        gtf_annotation.seqname = parameters_list[0]
        gtf_annotation.source = parameters_list[1]
        gtf_annotation.feature = parameters_list[2]
        gtf_annotation.start = int(parameters_list[3]) - 1
        gtf_annotation.end = int(parameters_list[4]) - 1
        gtf_annotation.score = parameters_list[5]
        if gtf_annotation.score != '.':
            gtf_annotation.score = int(gtf_annotation.score)
        gtf_annotation.strand = parameters_list[6]
        gtf_annotation.frame = parameters_list[7]
        if gtf_annotation.frame != '.':
            gtf_annotation.frame = int(gtf_annotation.frame)
        gtf_annotation.attributes_line = parameters_list[8]

        type_value_list = gtf_annotation.attributes_line.split('; ')
        for type_value in type_value_list:
            # for GTF file:
            if ' ' in type_value:
                temp = type_value.split(' ')
                gtf_annotation.attributes_list[temp[0]] = temp[1].replace('"', '')

            # for GFF file:
            elif '=' in type_value:
                temp = type_value.split('=')
                gtf_annotation.attributes_list[temp[0]] = str(temp[1])
        return gtf_annotation


    def get_line_from_gtf_annotation(self):
        return self.seqname + '\t' + self.source + '\t' + self.feature + '\t' + str(self.start + 1) + '\t' + str(self.end + 1) + \
               '\t' + str(self.score) + '\t' + self.strand + '\t' + str(self.frame) + '\t' + self.attributes_line