import matplotlib
matplotlib.use('Agg')

__author__ = 'lenk'

from abc import ABCMeta, abstractmethod

from pylab import *


class Alignment(object):
    """Class, which represent alignment"""

    def __init__(self):
        # psl or blast6 formats:
        self.format = ''

        self.score = '.'

        # Number of bases that don't match
        self.mismatches = 0
        self.strand = '.'
        # attributes of query
        self.query_fragment = self.QueryAttributes()
        # attributes of target
        self.target_fragment = self.TargetAttributes()

    # procedure for testing and creating manually:
    def create(self, strand, name, size, start, end):
        self.strand = strand
        self.query_fragment.create(name, size, start, end)
        self.target_fragment.create(name, size, start, end)
        self.score = size


    class Attributes:
        """Abstract class of Ouery's or Target's attributes, which represent alignment line"""

        __metaclass__ = ABCMeta

        def __init__(self):
            self.name = ''
            self.size = 0
            self.start = 0
            self.end = 0

    class QueryAttributes(Attributes):
        """Class of Query's attributes"""

        def create(self, name, size, start, end):
            self.name = name
            self.size = size
            self.start = start
            self.end = end


    class TargetAttributes(Attributes):
        """Class of Target's attributes"""

        def create(self, name, size, start, end):
            self.name = name
            self.size = size
            self.start = start
            self.end = end


class PSLFileAlignment(Alignment):
    """Class of PSL line, which represent alignment"""

    def __init__(self):
        Alignment.__init__(self)

        self.format = 'psl'

        # Number of bases that match that aren't repeats
        self.matches = 0
        # Number of bases that don't match
        # self.mismatches = 0
        # Number of bases that match but are part of repeats
        self.repmatches = 0
        # Number of 'N' bases
        self.n_num = 0
        # '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
        # The coordinates for a negative strand in a PSL line are handled in a special way:
        # in the qStart and qEnd fields, the coordinates indicate the position where the query
        # matches from the point of view of the forward strand, even when the match is on the
        # reverse strand. However, in the qStarts list, the coordinates are reversed.
        # self.strand = '+'
        # Number of blocks in the alignment (a block contains no gaps)
        self.blocks_num = 0
        # Comma-separated list of len of each block
        self.blocks_sizes = []
        # attributes of query
        self.query_fragment = self.QueryAttributes()
        # attributes of target
        self.target_fragment = self.TargetAttributes()


    class PSLAttributes(Alignment.Attributes):
        """Abstract class of Ouery's or Target's attributes, which represent alignment line in PSL-file"""

        __metaclass__ = ABCMeta

        def __init__(self):
            Alignment.Attributes.__init__(self)

            self.num_insert = 0
            self.base_insert = 0
            # self.name = ''
            # self.size = 0
            # self.start = 0
            # self.end = 0
            self.starts = []


    class QueryAttributes(PSLAttributes):
        """Class of Query's (aligned on genome transcripts) attributes, which represent line in PSL-file"""

        def set_from_psl_line(self, parameters_list, blocks_sizes):
            # Number of inserts in query
            self.num_insert = int(parameters_list[4])
            # Number of bases inserted in query
            self.base_insert = int(parameters_list[5])
            # Query sequence name
            self.name = parameters_list[9]
            # Query sequence size
            self.size = int(parameters_list[10])
            # Alignment start position in query
            self.start = int(parameters_list[11])
            # Alignment end position in query
            self.end = int(parameters_list[12]) - 1
            # Comma-separated list of starting positions of each block in query
            self.starts = parameters_list[19].split(',')[:-1]
            self.ends = []
            for i in range(len(self.starts)):
                self.starts[i] = int(self.starts[i])
                self.ends.append(self.starts[i] + blocks_sizes[i] - 1)


    class TargetAttributes(PSLAttributes):
        """Class of Target's (genome) attributes, which represent line in PSL-file"""

        def set_from_psl_line(self, parameters_list, blocks_sizes):
            #  Number of inserts in target
            self.num_insert = int(parameters_list[6])
            # Number of bases inserted in target
            self.base_insert = int(parameters_list[7])
            # Target sequence name
            self.name = parameters_list[13]
            # Target sequence size
            self.size = int(parameters_list[14])
            # Alignment start position in target
            self.start = int(parameters_list[15])
            # Alignment end position in target
            self.end = int(parameters_list[16]) - 1
            # Comma-separated list of starting positions of each block in target
            self.starts = parameters_list[20].split(',')[:-1]
            self.ends = []
            for i in range(len(self.starts)):
                self.starts[i] = int(self.starts[i])
                self.ends.append(self.starts[i] + blocks_sizes[i] - 1)


    # get attributes of class Alignment from line in PSL-file:
    @classmethod
    def get_alignment_from_psl_line(cls, text_line):
        psl_alignment = cls()
        psl_alignment.set_alignment_from_psl_line(text_line)
        return psl_alignment


    # set attributes of class Alignment from line in PSL-file:
    def set_alignment_from_psl_line(self, text_line):
        parameters_list = text_line.split('\t')
        self.matches = int(parameters_list[0])
        self.mismatches = int(parameters_list[1])
        self.repmatches = int(parameters_list[2])
        self.n_num = int(parameters_list[3])
        self.strand = parameters_list[8]
        self.blocks_num = int(parameters_list[17])
        self.blocks_sizes = parameters_list[18][:-1].split(',')
        for i in range(len(self.blocks_sizes)):
            self.blocks_sizes[i] = int(self.blocks_sizes[i])
        self.query_fragment.set_from_psl_line(parameters_list, self.blocks_sizes)
        self.target_fragment.set_from_psl_line(parameters_list, self.blocks_sizes)

        self.score = self.matches


    def get_psl_line_from_alignment(self):
        psl_line = str(self.matches) + '\t' + str(self.mismatches) + '\t' + str(self.repmatches) + '\t' + str(self.n_num) + '\t' + \
                   str(self.query_fragment.num_insert) + '\t' + str(self.query_fragment.base_insert) + '\t' + \
                   str(self.target_fragment.num_insert) + '\t' + str(self.target_fragment.base_insert) + '\t' + self.strand + '\t' + \
                   self.query_fragment.name + '\t' + str(self.query_fragment.size) + '\t' + str(self.query_fragment.start) + '\t' + str(self.query_fragment.end + 1) + '\t' + \
                   self.target_fragment.name + '\t' + str(self.target_fragment.size) + '\t' + str(self.target_fragment.start) + '\t' + str(self.target_fragment.end + 1) + '\t' + \
                   str(self.blocks_num) + '\t'

        for i_block in range(self.blocks_num):
            psl_line += str(self.blocks_sizes[i_block]) + ','
        psl_line += '\t'

        for i_block in range(self.blocks_num):
            psl_line += str(self.query_fragment.starts[i_block]) + ','
        psl_line += '\t'

        for i_block in range(self.blocks_num):
            psl_line += str(self.target_fragment.starts[i_block]) + ','
        # psl_line += '\n'

        return psl_line


    # get alignment contains blocks from iBlock0 to iBlock1:
    def get_split_alignment(self, i_block0, i_block1):
        split_alignment = PSLFileAlignment()
        # !!! TEMPORARY SOLUTION, needs ssw alignment:
        split_alignment.matches = 0
        for i_block in range(i_block1 - i_block0 + 1):
            split_alignment.matches += self.query_fragment.ends[i_block0 + i_block] - self.query_fragment.starts[i_block0 + i_block] + 1
        split_alignment.mismatches = 0
        split_alignment.repmatches = 0
        split_alignment.n_num = 0
        split_alignment.strand = self.strand
        split_alignment.blocks_num = i_block1 - i_block0 + 1
        split_alignment.blocks_sizes = self.blocks_sizes[i_block0:i_block1 + 1]

        split_alignment.score = split_alignment.matches

        # QUERY:
        # Number of inserts in query
        split_alignment.query_fragment.num_insert = 0
        # Number of bases inserted in query
        split_alignment.query_fragment.base_insert = 0
        for i_block in range(i_block1 - i_block0 + 1 - 1):
            query_base_insert = self.query_fragment.starts[i_block0 + i_block + 1] - self.query_fragment.ends[i_block0 + i_block] - 1
            split_alignment.query_fragment.base_insert += query_base_insert
            if query_base_insert > 0:
                split_alignment.query_fragment.num_insert += 1
        # Query sequence name
        split_alignment.query_fragment.name = self.query_fragment.name
        # Query sequence size
        split_alignment.query_fragment.size = self.query_fragment.size

        if self.strand == '+':
            # Alignment start position in query
            split_alignment.query_fragment.start = self.query_fragment.starts[i_block0]
            # Alignment end position in query
            split_alignment.query_fragment.end = self.query_fragment.ends[i_block1]
        else:
            # Alignment start position in query
            split_alignment.query_fragment.end = self.query_fragment.size - self.query_fragment.starts[i_block0] - 1
            # Alignment end position in query
            split_alignment.query_fragment.start = self.query_fragment.size - self.query_fragment.ends[i_block1] - 1

        # Comma-separated list of starting positions of each block in query
        split_alignment.query_fragment.starts = self.query_fragment.starts[i_block0:i_block1 + 1]
        split_alignment.query_fragment.ends = self.query_fragment.ends[i_block0:i_block1 + 1]

        # TARGET:
        #  Number of inserts in target
        split_alignment.target_fragment.num_insert = 0
        # Number of bases inserted in target
        split_alignment.target_fragment.base_insert = 0
        for i_block in range(i_block1 - i_block0 + 1 - 1):
            target_base_insert = self.target_fragment.starts[i_block0 + i_block + 1] - self.target_fragment.ends[i_block0 + i_block] - 1
            split_alignment.target_fragment.base_insert += target_base_insert
            if target_base_insert > 0:
                split_alignment.target_fragment.num_insert += 1
        # Target sequence name
        split_alignment.target_fragment.name = self.target_fragment.name
        # Target sequence size
        split_alignment.target_fragment.size = self.target_fragment.size

        # Alignment start position in target
        split_alignment.target_fragment.start = self.target_fragment.starts[i_block0]
        # Alignment end position in target
        split_alignment.target_fragment.end = self.target_fragment.ends[i_block1]

        # Comma-separated list of starting positions of each block in target
        split_alignment.target_fragment.starts = self.target_fragment.starts[i_block0:i_block1 + 1]
        split_alignment.target_fragment.ends = self.target_fragment.ends[i_block0:i_block1 + 1]

        ## DEBUGGING:
        #if split_alignment.targetFragment.start > split_alignment.targetFragment.end or \
        #    split_alignment.queryFragment.start > split_alignment.queryFragment.end:
        #    print(getAlignmentLine(split_alignment))
        #    import sys
        #    sys.exit(1)

        return split_alignment


class BLAST6FileAlignment(Alignment):
    """Class which represent alignment line by blast in outfmt 6"""

    def __init__(self):
        Alignment.__init__(self)

        self.format = 'blast6'

        self.line = ''

        self.pident = 0.0
        self.alen = 0
        # self.mismatches = 0
        self.gapopen = 0
        self.evalue = 0.0
        self.bitscore = 0.0
        # self.strand = ''

        # attributes of query
        self.query_fragment = self.QueryAttributes()
        # attributes of target
        self.target_fragment = self.TargetAttributes()


    class QueryAttributes(Alignment.QueryAttributes):
        """Class of Query's attributes, which represent line in BLAST6-file"""
        def set_from_blast6_line(self, parameters_list):
            # Query sequence name
            self.name = parameters_list[0]
            # Query sequence size
            # self.size = int(parameters_list[3])
            # Alignment start position in query
            self.start = int(parameters_list[6]) - 1
            # Alignment end position in query
            self.end = int(parameters_list[7]) - 1


    class TargetAttributes(Alignment.TargetAttributes):
        """Class of Target's attributes, which represent line in BLAST6-file"""
        def set_from_blast6_line(self, parameters_list):
            # Query sequence name
            self.name = parameters_list[1]
            # Query sequence size
            # self.size = int(parameters_list[4])
            # Alignment start position in query
            self.start = int(parameters_list[8]) - 1
            # Alignment end position in query
            self.end = int(parameters_list[9]) - 1
            if parameters_list[12] == 'minus':
                tmp = self.start
                self.start = self.end
                self.end = tmp


    # set attributes of class Alignment from line in blast6-file:
    def set_alignment_from_blast6_line(self, text_line):
        self.line = text_line.strip()

        parameters_list = self.line.split('\t')

        self.pident = float(parameters_list[2])
        self.alen = int(parameters_list[3])
        self.mismatches = int(parameters_list[4])
        self.gapopen = int(parameters_list[5])
        self.evalue = float(parameters_list[10])
        self.bitscore = float(parameters_list[11])
        tmp_strand = parameters_list[12]
        if tmp_strand == 'plus':
            self.strand == '+'
        else:
            self.strand == '-'

        self.query_fragment.set_from_blast6_line(parameters_list)
        self.target_fragment.set_from_blast6_line(parameters_list)

        self.score = self.bitscore


    # get attributes of class Alignment from line in blast6-file:
    @classmethod
    def get_alignment_from_blast6_line(cls, text_line):
        blast6_alignment = cls()
        blast6_alignment.set_alignment_from_blast6_line(text_line)
        return blast6_alignment


    # def addAlignmentFromKey(self, alignments, key):
    #     self.qseqid.append(alignments.qseqids[key])
    #     self.sseqid.append(alignments.sseqids[key])
    #     self.pident.append(alignments.pident[key])
    #     self.qlen.append(alignments.qlens[key])
    #     self.slen.append(alignments.slens[key])
    #     self.alen.append(alignments.alens[key])
    #     self.mismatches.append(alignments.mismatches[key])
    #     self.qstart.append(alignments.qstarts[key])
    #     self.qend.append(alignments.qends[key])
    #     self.sstart.append(alignments.sstarts[key])
    #     self.send.append(alignments.sends[key])
    #     self.evalue.append(alignments.evalues[key])
    #     self.bitscore.append(alignments.bitscores[key])


    # @classmethod
    # def getBestAlignments(cls, alignments):
    #     bestAlignments = BLAST6FileAlignment()
    #     iAlignment = 0
    #     while iAlignment in range(len(alignments.qseqids) - 1):
    #         min = abs(alignments.qlens[iAlignment]  - alignments.slens[iAlignment]) - \
    #               alignments.alens[iAlignment] + alignments.mismatches[iAlignment]
    #         key = iAlignment
    #         while alignments.qseqids[iAlignment] == alignments.qseqids[iAlignment + 1]:
    #             value = abs(alignments.qlens[iAlignment + 1] - alignments.slens[iAlignment + 1]) - \
    #                     alignments.alens[iAlignment + 1] + alignments.mismatches[iAlignment + 1]
    #             if value < min:
    #                 min = value
    #                 key = iAlignment + 1
    #             iAlignment += 1
    #             if iAlignment + 1 == len(alignments.qseqids):
    #                 break
    #         bestAlignments.addAlignmentFromKey(alignments, key)
    #         iAlignment += 1
    #     return bestAlignments


    # def get_similarity(self, outdir):
    #     similarityDistribution = {}
    #     querySimilarityDistribution = {}
    #     subjectSimilarityDistribution = {}
    #     for iAlignment in range(len(self.qseqid)):
    #         subjectSimilarity = (self.alen[iAlignment] - self.mismatches[iAlignment]) * 1.0 / self.slen[iAlignment]
    #         querySimilarity = (self.alen[iAlignment] - self.mismatches[iAlignment]) * 1.0 / self.qlen[iAlignment]
    #         similarity = min(querySimilarity, subjectSimilarity)
    #
    #         if subjectSimilarity not in subjectSimilarityDistribution:
    #             subjectSimilarityDistribution[subjectSimilarity] = 0
    #         subjectSimilarityDistribution[subjectSimilarity] += 1
    #
    #         if querySimilarity not in querySimilarityDistribution:
    #             querySimilarityDistribution[querySimilarity] = 0
    #         querySimilarityDistribution[querySimilarity] += 1
    #
    #         if similarity not in similarityDistribution:
    #             similarityDistribution[similarity] = 0
    #         similarityDistribution[similarity] += 1
    #
    #     figure()
    #     step = 0.1
    #     tempShowSimilarityDistribution = GeneralUtils.showDistribution(similarityDistribution, step)
    #
    #     bar(tempShowSimilarityDistribution.keys(), tempShowSimilarityDistribution.values(), width=step)
    #     title('Similarity distribution')
    #     xlabel('similarity')
    #     ylabel('number of transcripts')
    #     savefig('{}/Similarity.png'.format(outdir))
    #
    #     figure()
    #     title('Similarity distribution')
    #     step = 0.1
    #     tempShowSubjectSimilarityDistribution = GeneralUtils.showDistribution(subjectSimilarityDistribution, step)
    #     tempShowQuerySimilarityDistribution = GeneralUtils.showDistribution(querySimilarityDistribution, step)
    #     tempShowSimilarityDistribution = GeneralUtils.showDistribution(similarityDistribution, step)
    #
    #     plot(tempShowSubjectSimilarityDistribution.keys(), tempShowSubjectSimilarityDistribution.values(), '.', label='subject similarity')
    #     plot(tempShowQuerySimilarityDistribution.keys(), tempShowQuerySimilarityDistribution.values(), '^', label='query similarity')
    #     plot(tempShowSimilarityDistribution.keys(), tempShowSimilarityDistribution.values(), '^', label='similarity')
    #
    #     xBegin = min(tempShowSimilarityDistribution.keys())
    #     xEnd = max(tempShowSimilarityDistribution.keys())
    #     yBegin = min(tempShowSimilarityDistribution.values())
    #     yEnd = max(tempShowSimilarityDistribution.values())
    #     axis([xBegin - step, xEnd + step, yBegin - 1, yEnd + 1])
    #
    #     legend(fontsize='small')
    #
    #     xlabel('similarity')
    #     ylabel('number of transcripts')
    #
    #     savefig('{}/AllSimilarities.png'.format(outdir))


class SAMFileAlignment(Alignment):
    '''Docs: https://samtools.github.io/hts-specs/SAMv1.pdf'''

    def __init__(self):
        Alignment.__init__(self)

        self.format = 'sam'

        self.flag = 0
        self.mapq = 0
        self.cigar = ""
        self.rnext = "*"
        self.pnext = 0
        self.tlen = 0
        self.qual = "*"
        self.comments = "AS:i:0"

        # attributes of query
        self.query_fragment = self.QueryAttributes()
        # attributes of target
        self.target_fragment = self.TargetAttributes()

        self.cigar_commands = []


    class QueryAttributes(Alignment.QueryAttributes):

        def __init__(self):
            Alignment.QueryAttributes.__init__(self)

            self.seq = ""
            self.aligned_seq = ''

            self.starts = []
            self.ends = []

        def set_from_sam_line(self, parameters_list, cigar_commands):
            self.name = parameters_list[0]
            self.seq = parameters_list[9]

            curr_pos = 0
            for command in cigar_commands:
                if command[1] == 'M' or command[1] == '=' or command[1] == 'X':
                    self.starts.append(curr_pos)
                    self.ends.append(curr_pos + command[0] - 1)
                    curr_pos = self.ends[-1] + 1
                elif command[1] == 'S' or command[1] == 'I':
                    curr_pos += command[0]


    class TargetAttributes(Alignment.TargetAttributes):

        def __init__(self):
            Alignment.TargetAttributes.__init__(self)

            self.aligned_seq = ''

            self.starts = []
            self.ends = []

        def set_from_sam_line(self, parameters_list, cigar_commands):
            self.name = parameters_list[2]
            self.start = int(parameters_list[3]) - 1

            curr_pos = self.start
            for command in cigar_commands:
                if command[1] == 'M' or command[1] == '=' or command[1] == 'X':
                    self.starts.append(curr_pos)
                    self.ends.append(curr_pos + command[0] - 1)
                    curr_pos = self.ends[-1] + 1
                elif command[1] == 'D' or command[1] == 'N':
                    curr_pos += command[0]


    def set_alignment_from_sam_line(self, text_line, logger, reference_dict):
        self.line = text_line

        parameters_list = text_line.split('\t')

        self.flag = int(parameters_list[1])
        if self.flag == 16:
            self.strand = '-'
        elif self.flag == 0:
            self.strand = '+'


        self.mapq = int(parameters_list[4])

        self.cigar = parameters_list[5]
        self.cigar_commands = self.get_commands_from_cigar(logger)

        self.rnext = parameters_list[6]
        self.pnext = int(parameters_list[7])
        self.tlen = int(parameters_list[8])
        self.qual = parameters_list[10]

        self.query_fragment.set_from_sam_line(parameters_list, self.cigar_commands)
        self.target_fragment.set_from_sam_line(parameters_list, self.cigar_commands)

        if reference_dict is not None:
            self.set_aligned_seqs(reference_dict, logger)


    @classmethod
    def get_alignment_from_sam_line(cls, text_line, logger, reference_dict):
        sam_alignment = cls()
        sam_alignment.set_alignment_from_sam_line(text_line, logger, reference_dict)
        return sam_alignment


    def get_commands_from_cigar(self, logger):
        commands = []
        i = 0
        curr_num = ''
        curr_command = ''
        while i < len(self.cigar):
            if self.cigar[i] in ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']:
                curr_command = self.cigar[i]
                commands.append((int(curr_num), curr_command))
                curr_num = ''
            elif self.cigar[i] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                curr_num += self.cigar[i]
            else:
                logger.warning('Strange CIGAR string contains {}'.format(self.cigar[i]))
            i += 1
        return commands


    def set_aligned_seqs(self, reference_dict, logger):
        # print self.cigar, self.cigar_commands
        # print 'query starts: ', self.query_fragment.starts
        # print 'query ends: ', self.query_fragment.ends
        # print 'target_starts: ', self.target_fragment.starts
        # print 'target_ends: ', self.target_fragment.ends
        for i_pos in range(len(self.query_fragment.starts)):
            self.query_fragment.aligned_seq += self.query_fragment.seq[self.query_fragment.starts[i_pos]:self.query_fragment.ends[i_pos] + 1]
            self.target_fragment.aligned_seq += reference_dict[self.target_fragment.name][self.target_fragment.starts[i_pos]:self.target_fragment.ends[i_pos] + 1]
        # print 'query_seq: ', self.query_fragment.seq
        # print 'query_aligned_seq: ', query_aligned_seq
        # print 'target_aligned_seq: ', target_aligned_seq
        if len(self.query_fragment.aligned_seq) != len(self.target_fragment.aligned_seq):
            logger.warning('Unconsistent alignment lengths: q_alen={}, t_alen={}, cigar={}'.format(len(self.query_fragment.aligned_seq), len(self.target_fragment.aligned_seq), self.cigar))
            print self.target_fragment.starts, self.target_fragment.ends, len(reference_dict[self.target_fragment.name])
            print self.line
            print '\n'