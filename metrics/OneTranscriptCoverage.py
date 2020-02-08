__author__ = 'lenk'

from general import UtilsCoverage


class OneTranscriptCoverage(object):
    """Class, which represent coverage of internal aligned blocks in aligned transcript by annotated isoforms and exons"""

    def __init__(self, ids_isoforms, blocks_num):
        # number of blocks in aligned transcripts:
        self.blocks_num = 0

        # length of aligned blocks of aligned transcript:
        self.alignment_len = 0
        # length of aligned transcript:
        self.transcript_len = 0

        # DISTRIBUTIONS FOR TRANSCRIPTS ALIGNMENTS:
        # coverage by annotated isoforms:
        # number of covered bases of transcript by internal annotated isoforms:
        self.covered_bases = {}
        # percentage of covered bases of transcript of the aligned part by internal annotated isoforms:
        self.covered_fraction_whole_transcript = {}
        # percentage of covered bases of transcript of the transcript length by internal annotated isoforms:
        self.covered_fraction_aligned_part = {}

        # coverage of blocks:
        # number of covered bases for each blocks:
        self.covered_bases_blocks = {}
        # percentage of covered bases for each block:
        self.covered_fraction_blocks = {}
        # average over blocks in aligned transcript percentage of covered bases of block by internal annotated isoform:
        self.avg_covered_fraction_block = {}

        # number and percentage of well-covered by internal annotated isoforms blocks in aligned transcript:
        self.num_well_covered_blocks = {}
        self.percentage_well_covered_blocks = {}
        # number and percentage of fully-covered by internal annotated isoforms blocks in aligned transcript:
        self.num_fully_covered_blocks = {}
        self.percentage_fully_covered_blocks = {}

        # id isoform mapped to transcript:
        self.id_mapped_isoform = None

        self.blocks_num = blocks_num
        for id_isoform in ids_isoforms:
            self.covered_bases[id_isoform] = 0

            self.avg_covered_fraction_block[id_isoform] = 0.0

            self.num_well_covered_blocks[id_isoform] = 0
            self.num_fully_covered_blocks[id_isoform] = 0

            self.covered_bases_blocks[id_isoform] = []
            self.covered_fraction_blocks[id_isoform] = []
            for i_block in range(self.blocks_num):
                self.covered_bases_blocks[id_isoform].append(0)
                self.covered_fraction_blocks[id_isoform].append(0.0)


    # update coverage of blocks:
    def update_transcript_coverage(self, id_isoform, block_cov_pos):
        for i_block in block_cov_pos:
            for i_pos in range(len(block_cov_pos[i_block])):
                start_coverage = block_cov_pos[i_block][i_pos][0]
                end_coverage = block_cov_pos[i_block][i_pos][1]

                self.covered_bases_blocks[id_isoform][i_block] += end_coverage - start_coverage + 1


    # set coverage means get all attributes of coverage of alignedTranscript:
    def get_transcript_coverage(self, aligned_transcript, internal_isoforms_coverage, WELL_FULLY_COVERAGE_THRESHOLDS):
        self.transcript_len = aligned_transcript.alignment.query_fragment.size
        self.alignment_len = aligned_transcript.alignment_len

        for i_isoform in range(len(aligned_transcript.internal_isoforms)):
            id_isoform = aligned_transcript.internal_isoforms[i_isoform].id

            # get metrics of coverage excluded covered bases of blocks:
            for i_block in range(self.blocks_num):
                self.covered_bases[id_isoform] += self.covered_bases_blocks[id_isoform][i_block]

                self.covered_fraction_blocks[id_isoform][i_block] = self.covered_bases_blocks[id_isoform][i_block] * 1.0 / aligned_transcript.alignment.blocks_sizes[i_block]

                self.avg_covered_fraction_block[id_isoform] += self.covered_fraction_blocks[id_isoform][i_block]

                if self.covered_fraction_blocks[id_isoform][i_block] >= WELL_FULLY_COVERAGE_THRESHOLDS.well_block_threshold:
                    self.num_well_covered_blocks[id_isoform] += 1
                if self.covered_fraction_blocks[id_isoform][i_block] >= WELL_FULLY_COVERAGE_THRESHOLDS.fully_block_threshold:
                    self.num_fully_covered_blocks[id_isoform] += 1

            self.covered_fraction_whole_transcript[id_isoform] = self.covered_bases[id_isoform] * 1.0 / self.transcript_len
            self.covered_fraction_aligned_part[id_isoform] = self.covered_bases[id_isoform] * 1.0 / self.alignment_len

            self.avg_covered_fraction_block[id_isoform] /= self.blocks_num

            self.percentage_well_covered_blocks[id_isoform] = self.num_well_covered_blocks[id_isoform] * 1.0 / self.blocks_num
            self.percentage_fully_covered_blocks[id_isoform] = self.num_fully_covered_blocks[id_isoform] * 1.0 / self.blocks_num

        # get mapped to transcript isoform (TEMPORARY: can have better solution):
        if len(aligned_transcript.internal_isoforms) != 0:
            self.id_mapped_isoform = sorted(UtilsCoverage.get_ids_best_mapped(self.covered_bases, internal_isoforms_coverage.assembled_fraction))[0]