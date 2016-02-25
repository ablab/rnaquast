__author__ = 'letovesnoi'

import os

class TXTMetricsReport():
    """Class which generate txt reports"""

    def __init__(self, txt_report_dir, transcripts_metrics, db_genes_metrics, reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS):
        logger.print_timestamp('  ')
        logger.info('  Getting TXT report...')

        self.txt_reports_dir = txt_report_dir

        if db_genes_metrics is not None:
            self.path_txt_database_metrics = os.path.join(self.txt_reports_dir, 'database_metrics.txt')
            self.get_db_genes_metrics_report(db_genes_metrics, logger, PRECISION)

        if reads_coverage is not None:
            self.path_txt_reads_coverage_metrics = self.path_txt_database_metrics
            self.get_reads_coverage_report(reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

            self.path_txt_relative_database_coverage = os.path.join(self.txt_reports_dir, 'relative_database_coverage.txt')
            self.get_relative_database_coverage_report(transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

        if transcripts_metrics[0].basic_metrics is not None:
            self.path_txt_basic = os.path.join(self.txt_reports_dir, 'basic_metrics.txt')
            self.get_basic_metrics_report(transcripts_metrics, logger, TRANSCRIPT_LENS, PRECISION)

        if transcripts_metrics[0].simple_metrics is not None:
            self.path_txt_alignment = os.path.join(self.txt_reports_dir, 'alignment_metrics.txt')
            self.get_alignment_metrics_report(transcripts_metrics, logger, PRECISION)

            self.path_txt_misassemblies = os.path.join(self.txt_reports_dir, 'misassemblies.txt')
            self.get_misassemblies_report(transcripts_metrics, logger)

        if transcripts_metrics[0].assembly_completeness_metrics is not None:
                self.path_txt_sensitivity = os.path.join(self.txt_reports_dir, 'sensitivity.txt')
                self.get_sensitivity_report(transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

        if transcripts_metrics[0].assembly_correctness_metrics is not None:
            self.path_txt_specificity = os.path.join(self.txt_reports_dir, 'specificity.txt')
            self.get_specificity_report(transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

        logger.info('  Done.')


    # Gene database metrics report
    # Number of genes / protein coding genes
    # Number of isoforms / protein coding isoforms
    # Total number of exons / introns
    # Total / average isoform length, bp
    # Total / average exon length, bp
    # Total / average intron length, bp
    # Average / max number of exons per isoform
    # Average number of introns per isoform
    def get_db_genes_metrics_report(self, db_genes_metrics, logger, PRECISION):
        logger.info('    Getting GENE DATABASE METRICS report...')

        with open(self.path_txt_database_metrics, 'w') as fout:
            fout.write(' == GENE DATABASE METRICS ==\n')
            fout.write('{:<100}'.format('Genes') + str(db_genes_metrics.genes_num) + '\n')
            # gff source field unfortunately contain database name instead of transcripts type protein_coding:
            if db_genes_metrics.protein_coding_genes_num != 0:
                fout.write('{:<100}'.format('Protein coding genes') + str(db_genes_metrics.protein_coding_genes_num) + '\n\n')

            fout.write('{:<100}'.format('Isoforms') + str(db_genes_metrics.isoforms_num) + '\n')
            # gff source field unfortunately contain database name instead of transcripts type protein_coding:
            if db_genes_metrics.protein_coding_isoforms_num != 0:
                fout.write('{:<100}'.format('Protein coding isoforms') + str(db_genes_metrics.protein_coding_isoforms_num) + '\n\n')

            fout.write('{:<100}'.format('Exons') + str(db_genes_metrics.tot_exons_num) + '\n\n')

            fout.write('{:<100}'.format('Introns') + str(db_genes_metrics.tot_introns_num) + '\n\n')

            fout.write('{:<100}'.format('Avg. length of all isoforms') + str(round(db_genes_metrics.avg_isoform_len, PRECISION)) + '\n\n')
            fout.write('{:<100}'.format('Total length of all isoforms') + str(db_genes_metrics.tot_isoforms_len) + '\n')

            fout.write('{:<100}'.format('Avg. exon length') + str(round(db_genes_metrics.avg_exon_len, PRECISION)) + '\n\n')

            fout.write('{:<100}'.format('Avg. intron length') + str(round(db_genes_metrics.avg_intron_len, PRECISION)) + '\n\n')

            fout.write('{:<100}'.format('Avg. number of exons per isoform') + str(round(db_genes_metrics.avg_exons_num, PRECISION)) + '\n')
            fout.write('{:<100}'.format('Max. number of exons per isoform') + str(db_genes_metrics.max_exons_num) + '\n\n')

            # fout.write('{:<100}'.format('Avg. number of introns per isoform') + str(round(basic_isoforms_metrics.avg_introns_num, precision)) + '\n\n')

            # fout.write('ISOFORMS:\n')
            # fout.write('{:<100}'.format('Total length with introns') + str(basic_isoforms_metrics.tot_len_w_introns) + '\n')
            # fout.write('{:<100}'.format('Min. length with introns') + str(basic_isoforms_metrics.min_len_w_introns) + '\n')
            # fout.write('{:<100}'.format('Max. length with introns') + str(basic_isoforms_metrics.max_len_w_introns) + '\n')
            # fout.write('{:<100}'.format('Avg. length with introns') + str(round(basic_isoforms_metrics.avg_len_w_introns, precision)) + '\n\n')

            # fout.write('{:<100}'.format('Min. isoform length') + str(basic_isoforms_metrics.min_len_wout_introns) + '\n')
            # fout.write('{:<100}'.format('Max. isoform length') + str(basic_isoforms_metrics.max_len_wout_introns) + '\n')

            # fout.write('EXONS:\n')
            # fout.write('{:<100}'.format('Min. exons per isoform') + str(basic_isoforms_metrics.min_exons_num) + '\n')

            # fout.write('{:<100}'.format('Min. exon length') + str(basic_isoforms_metrics.min_exon_len) + '\n')
            # fout.write('{:<100}'.format('Max. exon length') + str(basic_isoforms_metrics.max_exon_len) + '\n')

            # fout.write('INTRONS:\n')
            # fout.write('{:<100}'.format('Min. introns per isoform') + str(basic_isoforms_metrics.min_introns_num) + '\n')
            # fout.write('{:<100}'.format('Max. introns per isoform') + str(basic_isoforms_metrics.max_introns_num) + '\n')

            # fout.write('{:<100}'.format('Min. intron length') + str(basic_isoforms_metrics.min_intron_len) + '\n')
            # fout.write('{:<100}'.format('Max. intron length') + str(basic_isoforms_metrics.max_intron_len) + '\n')
            # fout.write('{:<100}'.format('Total introns length [bp]') + str(basic_isoforms_metrics.tot_introns_len) + '\n')

            logger.info('      saved to {}'.format(self.path_txt_database_metrics))


    def get_reads_coverage_report(self, reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION):
        logger.info('    Getting COVERAGE BY READS report...')

        with open(self.path_txt_reads_coverage_metrics, 'a') as fout:
            fout.write(' == COVERAGE BY READS ==\n')
            fout.write('{:<100}'.format('Database coverage') + str(round(reads_coverage.fraction_annotation_mapped_by_reads, PRECISION)) + '\n\n')

            fout.write('{:<100}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered genes') + str(reads_coverage.num_well_expressed_genes) + '\n')
            fout.write('{:<100}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered genes') + str(reads_coverage.num_fully_expressed_genes) + '\n\n')

            fout.write('{:<100}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered isoforms') + str(reads_coverage.num_well_expressed_isoforms) + '\n')
            fout.write('{:<100}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered isoforms') + str(reads_coverage.num_fully_expressed_isoforms) + '\n\n')

            fout.write('{:<100}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-covered exons') + str(reads_coverage.num_well_expressed_exons) + '\n')
            fout.write('{:<100}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold * 100)) + '%-covered exons') + str(reads_coverage.num_fully_expressed_exons) + '\n\n')

        logger.info('      saved to {}'.format(self.path_txt_reads_coverage_metrics))


    def get_relative_database_coverage_report(self, transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION):
        logger.info('    Getting RELATIVE DATABASE COVERAGE report...')

        name_str = '{:<50}'.format('METRICS/TRANSCRIPTS')

        database_coverage_str = '{:<50}'.format('Relative database coverage')

        well_assembled_genes_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled genes')
        fully_assembled_genes_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled genes')

        well_covered_genes_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered genes')
        fully_covered_genes_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered genes')


        well_assembled_isoforms_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled isoforms')
        fully_assembled_isoforms_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled isoforms')

        well_covered_isoforms_str = '{:<50}'.format(str('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100))) + '%-covered isoforms')
        fully_covered_isoforms_str = '{:<50}'.format(str('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100))) + '%-covered isoforms')


        well_assembled_exons_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-assembled exons')
        fully_assembled_exons_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold * 100)) + '%-assembled exons')

        well_covered_exons_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-covered exons')
        fully_covered_exons_str = '{:<50}'.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold * 100)) + '%-covered exons')

        for i_transcripts in range(len(transcripts_metrics)):
            relative_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.relative_database_coverage

            name_str += '{:<25}'.format(transcripts_metrics[i_transcripts].label)

            database_coverage_str += '{:<25}'.format(round(relative_metrics.database_coverage, PRECISION))

            well_assembled_genes_str += '{:<25}'.format(round(relative_metrics.percentage_well_assembled_genes, PRECISION))
            fully_assembled_genes_str += '{:<25}'.format(round(relative_metrics.percentage_fully_assembled_genes, PRECISION))

            well_covered_genes_str += '{:<25}'.format(round(relative_metrics.percentage_well_covered_genes, PRECISION))
            fully_covered_genes_str += '{:<25}'.format(round(relative_metrics.percentage_fully_covered_genes, PRECISION))


            well_assembled_isoforms_str += '{:<25}'.format(round(relative_metrics.percentage_well_assembled_isoforms, PRECISION))
            fully_assembled_isoforms_str += '{:<25}'.format(round(relative_metrics.percentage_fully_assembled_isoforms, PRECISION))

            well_covered_isoforms_str += '{:<25}'.format(round(relative_metrics.percentage_well_covered_isoforms, PRECISION))
            fully_covered_isoforms_str += '{:<25}'.format(round(relative_metrics.percentage_fully_covered_isoforms, PRECISION))


            well_assembled_exons_str += '{:<25}'.format(round(relative_metrics.percentage_well_assembled_exons, PRECISION))
            fully_assembled_exons_str += '{:<25}'.format(round(relative_metrics.percentage_fully_assembled_exons, PRECISION))

            well_covered_exons_str += '{:<25}'.format(round(relative_metrics.percentage_well_covered_exons, PRECISION))
            fully_covered_exons_str += '{:<25}'.format(round(relative_metrics.percentage_fully_covered_exons, PRECISION))

        with open(self.path_txt_relative_database_coverage, 'w') as out_handle:
            out_handle.write(' == RELATIVE DATABASE COVERAGE ==\n')
            out_handle.write(name_str + '\n')

            out_handle.write(database_coverage_str + '\n\n')

            out_handle.write(well_assembled_genes_str + '\n')
            out_handle.write(fully_assembled_genes_str + '\n')
            out_handle.write(well_covered_genes_str + '\n')
            out_handle.write(fully_covered_genes_str + '\n\n')

            out_handle.write(well_assembled_isoforms_str + '\n')
            out_handle.write(fully_assembled_isoforms_str + '\n')
            out_handle.write(well_covered_isoforms_str + '\n')
            out_handle.write(fully_covered_isoforms_str + '\n\n')

            out_handle.write(well_assembled_exons_str + '\n')
            out_handle.write(fully_assembled_exons_str + '\n')
            out_handle.write(well_covered_exons_str + '\n')
            out_handle.write(fully_covered_exons_str + '\n\n')

        logger.info('      saved to {}'.format(self.path_txt_relative_database_coverage))


    # Basic transcript metrics (calculated without reference genome and gene database)
    # Number of assembled transcripts
    # #transcripts > 500 bp
    # #transcripts > 1000 bp
    # Average length of assembled transcripts
    # Longest transcript
    # Total length
    # Transcript N50: a maximal number X, such that the total length of all transcripts longer than X bp is at least 50% of the total length of all transcripts
    def get_basic_metrics_report(self, transcripts_metrics, logger, TRANSCRIPT_LENS, PRECISION):
        logger.info('    Getting BASIC TRANSCRIPTS METRICS report...')

        name_str = '{:<50}'.format('METRICS/TRANSCRIPTS')

        num_transcripts_str = '{:<50}'.format('Transcripts')
        num_transcripts_500_str = '{:<50}'.format('Transcripts > {} bp'.format(str(TRANSCRIPT_LENS[0])))
        num_transcripts_1000_str = '{:<50}'.format('Transcripts > {} bp'.format(str(TRANSCRIPT_LENS[1])))
        avg_len_transcripts_str = '{:<50}'.format('Average length of assembled transcripts')
        longest_transcript_str = '{:<50}'.format('Longest transcript')
        tot_len_str = '{:<50}'.format('Total length')
        n50_str = '{:<50}'.format('Transcript N50')

        for i_transcripts in range(len(transcripts_metrics)):
            basic_metrics = transcripts_metrics[i_transcripts].basic_metrics

            name_str += '{:<25}'.format(transcripts_metrics[i_transcripts].label)

            num_transcripts_str += '{:<25}'.format(basic_metrics.number)
            num_transcripts_500_str += '{:<25}'.format(basic_metrics.num_transcripts_500)
            num_transcripts_1000_str += '{:<25}'.format(basic_metrics.num_transcripts_1000)
            avg_len_transcripts_str += '{:<25}'.format(round(basic_metrics.avg_len, PRECISION))
            longest_transcript_str += '{:<25}'.format(basic_metrics.max_len)
            tot_len_str += '{:<25}'.format(basic_metrics.tot_len)
            n50_str += '{:<25}'.format(basic_metrics.n50)

            # fout.write('{:<100}'.format('Min. length') + str(self.min_len) + '\n')
            # if args.annotation is not None:
                # number with length inside isoform length range:
                # fout.write('{:<100}'.format('Outliers') + str(self.num_outliers_len) + '\n')
            # fout.write('\n')

        fout = open(self.path_txt_basic, 'w')

        fout.write(' == BASIC TRANSCRIPTS METRICS (calculated without reference genome and gene database) ==\n')
        fout.write(name_str + '\n')
        fout.write(num_transcripts_str + '\n\n')
        fout.write(num_transcripts_500_str + '\n')
        fout.write(num_transcripts_1000_str + '\n\n')
        fout.write(avg_len_transcripts_str + '\n')
        fout.write(longest_transcript_str + '\n')
        fout.write(tot_len_str + '\n\n')
        fout.write(n50_str + '\n')

        fout.close()

        logger.info('      saved to {}'.format(self.path_txt_basic))


    # Alignment metrics (calculated with reference genome but without gene database)
    # Transcripts = Unaligned + Aligned = Unaligned + (Uniquely aligned + Multiply aligned + Misassembly candidates reported by BLAT)
    # Aligned: the number of transcripts having at least 1 significant alignment.
    # Uniquely aligned: the number of transcripts having a single significant alignment.
    # Multiply aligned: the number of transcripts having 2 or more significant alignments.
    # Misassembly candidates reported by BLAT: transcripts are aligned to the reference genome with BLAT and discordant partial alignments are selected as misassembly candidates.
    # Unaligned: the number of transcripts without any significant alignments.

    # Alignment metrics for uniquely aligned transcripts
    # Genome fraction: the total number of aligned bases divided by the reference genome length.
    # Duplication ratio w/o database: the total number of aligned bases in transcripts divided by the total number of bases covered by these alignments in reference genome. Note that duplication ratio can be increased for organisms having simultaneous expression of various isoforms of the same gene.
    # Average aligned fraction. Aligned fraction for a single transcript is defined as total number of aligned bases divided by the total transcript length.
    # Average alignment length.
    # Average blocks per alignment. A block is defined as a continuous alignment fragment without indels.
    # Average block length (see above).
    # Average mismatches per transcript: average number of single nucleotide differences with reference genome per transcript.
    # NA50: N50 for alignments
    def get_alignment_metrics_report(self, transcripts_metrics, logger, PRECISION):
        logger.info('    Getting ALIGNMENT METRICS report...')

        name_str = '{:<50}'.format('METRICS/TRANSCRIPTS')

        num_transcripts_str = '{:<50}'.format('Transcripts')
        num_aligned_str = '{:<50}'.format('Aligned')
        num_uniquely_aligned_str = '{:<50}'.format('Uniquely aligned')
        num_multiply_aligned_str = '{:<50}'.format('Multiply aligned')
        num_misassembled_by_blat_str = '{:<50}'.format('Misassembly candidates reported by BLAT')
        num_unaligned_str = '{:<50}'.format('Unaligned')

        # Alignment metrics for uniquely aligned transcripts
        # genome_fraction_str = '{:<50}'.format('Genome fraction')
        # duplication_ratio_wo_database_str = '{:<50}'.format('Duplication ratio w/o database')
        avg_aligned_fraction_str = '{:<50}'.format('Average aligned fraction')
        avg_alignment_len_str = '{:<50}'.format('Average alignment length')
        avg_blocks_per_alignment = '{:<50}'.format('Average blocks per alignment')
        avg_block_len_str = '{:<50}'.format('Average block length')
        avg_mismatches_per_transcript_str = '{:<50}'.format('Average mismatches per transcript')
        na50_str = '{:<50}'.format('NA50')
        # avg_gaps_per_alignment_str = '{:<50}'.format('Average gaps per alignment')
        # avg_gap_length_str = '{:<50}'.format('Average gap length')

        for i_transcripts in range(len(transcripts_metrics)):
            basic_metrics = transcripts_metrics[i_transcripts].basic_metrics
            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics

            name_str += '{:<25}'.format(transcripts_metrics[i_transcripts].label)

            num_transcripts_str += '{:<25}'.format(basic_metrics.number)
            num_aligned_str += '{:<25}'.format(simple_metrics.num_aligned)
            num_uniquely_aligned_str += '{:<25}'.format(simple_metrics.num_unique_aligned)
            num_multiply_aligned_str += '{:<25}'.format(simple_metrics.num_mul_aligned)
            num_misassembled_by_blat_str += '{:<25}'.format(simple_metrics.num_misassembled_by_blat)
            num_unaligned_str += '{:<25}'.format(simple_metrics.num_unaligned)

            # number of alignments in PSL-file:
            # fout.write('{:<100}'.format('Alignments') + str(self.num_alignments) + '\n\n')

            # number of not misassembled aligned transcripts (best union is single union):
            # fout.write('{:<100}'.format('Not misassembled aligned') + str(self.num_non_misassembled) + '\n\n')

            # TEMPORARY FOR TEST BLAT ALIGNER:
            # fout.write('{:<100}'.format('Multiple equal query aligned') + str(simple_metrics.num_mul_equal_query_aligned) + '\n')
            # fout.write('{:<100}'.format('Multiple equal query aligned over misassemblies') + str(simple_metrics.num_mul_equal_query_aligned_mis) + '\n\n')

            # genome_fraction_str += '{:<25}'.format(round(simple_metrics.fraction_genome_mapped, precision))
            # duplication_ratio_wo_database_str += '{:<25}'.format(round(simple_metrics.avg_duplication_ratio, precision))
            avg_aligned_fraction_str += '{:<25}'.format(round(simple_metrics.avg_fraction, PRECISION))
            avg_alignment_len_str += '{:<25}'.format(round(simple_metrics.avg_alignment_len, PRECISION))
            avg_blocks_per_alignment += '{:<25}'.format(round(simple_metrics.avg_blocks_num, PRECISION))
            avg_block_len_str += '{:<25}'.format(round(simple_metrics.avg_block_len, PRECISION))
            avg_mismatches_per_transcript_str += '{:<25}'.format(round(simple_metrics.avg_mismatch_num, PRECISION))
            na50_str += '{:<25}'.format(simple_metrics.na50)
            # avg_gaps_per_alignment_str +=  '{:<25}'.format(round(simple_metrics.avg_qgap_num, precision))
            # avg_gap_length_str += '{:<25}'.format(round(simple_metrics.avg_qgap_len, precision))

            # LENGTH:
            # fout.write('{:<100}'.format('Total length') + str(self.tot_len) + '\n')
            # fout.write('{:<100}'.format('Min. length') + str(self.min_len) + '\n')
            # fout.write('{:<100}'.format('Max. length') + str(self.max_len) + '\n')
            # fout.write('{:<100}'.format('Avg. length') + str(round(self.avg_len, precision)) + '\n')
            # if args.annotation is not None:
                # number of aligned transcripts len out isoforms without introns length:
                # fout.write('{:<100}'.format('Outliers') + str(self.num_outliers_len) + '\n')
            # fout.write('\n')

            # REFERENCE LENGTH (from begin to end of alignment):
            # fout.write('{:<100}'.format('Total reference length') + str(self.tot_ref_len) + '\n')
            # fout.write('{:<100}'.format('Min. reference length') + str(self.min_ref_len) + '\n')
            # fout.write('{:<100}'.format('Max. reference length') + str(self.max_ref_len) + '\n')
            # fout.write('{:<100}'.format('Avg. reference length') + str(round(self.avg_ref_len, precision)) + '\n')
            # if args.annotation is not None:
                # number of aligned transcripts with length more then max isoform with introns length or less then min isoform with introns length:
                # fout.write('{:<100}'.format('Outliers') + str(self.num_outliers_ref_len) + '\n')
            # fout.write('\n')

            # ALIGNMENT LENGTH (summary length of aligned blocks):
            # fout.write('{:<100}'.format('Total alignment length') + str(self.tot_alignment_len) + '\n')
            # fout.write('{:<100}'.format('Min. alignment length') + str(self.min_alignment_len) + '\n')
            # fout.write('{:<100}'.format('Max. alignment length') + str(self.max_alignment_len) + '\n')

            # if args.annotation is not None:
                # number of aligned transcripts with alignment length less then min isoform without introns length or more then max isoform without introns length:
                # fout.write('{:<100}'.format('Outliers') + str(self.num_outliers_alignment_len) + '\n')
            # fout.write('\n')

            # ALIGNMENT FRACTION:
            # fout.write('{:<100}'.format('Min. alignment fraction') + str(round(self.min_fraction, precision)) + '\n')
            # fout.write('{:<100}'.format('Max. alignment fraction') + str(round(self.max_fraction, precision)) + '\n')

            # BLOCKS NUMBER:
            # fout.write('{:<100}'.format('Total blocks number') + str(self.tot_blocks_num) + '\n')
            # fout.write('{:<100}'.format('Min. blocks number') + str(self.min_blocks_num) + '\n')
            # fout.write('{:<100}'.format('Max. blocks number') + str(self.max_blocks_num) + '\n')
            # if args.annotation is not None:
                # number and percentage of aligned transcripts with number of blocks more then max number of exons per isoform or less then min:
                # fout.write('{:<100}'.format('Outliers') + str(self.num_outliers_blocks_num) + '\n')
            # fout.write('\n')

            # BLOCKS LENGTH:
            # fout.write('{:<100}'.format('Min. block length') + str(self.min_block_len) + '\n')
            # fout.write('{:<100}'.format('Max. block length') + str(self.max_block_len) + '\n')
            # if args.annotation is not None:
                # number and percentage of blocks with length more then max length of exon and less then min length of exon:
                # fout.write('{:<100}'.format('Outliers') + str(self.num_outliers_blocks_len) + '\n')
            # fout.write('\n')

            # MISMATCHES:
            # fout.write('{:<100}'.format('Total mismatch number') + str(self.tot_mismatch_num) + '\n')

            # TARGET GAPS NUMBER:
            # fout.write('{:<100}'.format('Total target insert number') + str(self.tot_tgap_num) + '\n')
            # fout.write('{:<100}'.format('Min. target insert number') + str(self.min_tgap_num) + '\n')
            # fout.write('{:<100}'.format('Max. target insert number') + str(self.max_tgap_num) + '\n')
            # fout.write('{:<100}'.format('Avg. target insert number') + str(round(self.avg_tgap_num, precision)) + '\n')
            # if args.annotation is not None:
                # number and percentage of transcripts with number of target gaps more then max number of introns per isoform and less then min:
                # fout.write('{:<100}'.format('Outliers') + str(self.num_outliers_tgap_num) + '\n')
            # fout.write('\n')

            # TARGET GAPS LENGTH:
            # fout.write('{:<100}'.format('Total target insert length') + str(self.tot_tgap_len) + '\n')
            # fout.write('{:<100}'.format('Min. target insert length') + str(self.min_tgap_len) + '\n')
            # fout.write('{:<100}'.format('Max. target insert length') + str(self.max_tgap_len) + '\n')
            # fout.write('{:<100}'.format('Avg. target insert length') + str(round(self.avg_tgap_len, precision)) + '\n')
            # if args.annotation is not None:
                # number and percentage of target gaps with length more then max length of intron and less then min length of intron:
                # fout.write('{:<100}'.format('Outliers') + str(self.num_outliers_tgap_len) + '\n')
            # fout.write('\n')

            # QUERY GAPS NUMBER:
            # fout.write('{:<100}'.format('Total query insert number') + str(self.tot_qgap_num) + '\n')
            # QUERY GAPS LENGTH:
            # fout.write('{:<100}'.format('Total query insert length') + str(self.tot_qgap_len) + '\n')

        fout = open(self.path_txt_alignment, 'w')

        fout.write(' == ALIGNMENT METRICS (calculated with reference genome but without gene database) ==\n')
        fout.write(name_str + '\n')
        fout.write(num_transcripts_str + '\n\n')
        fout.write(num_aligned_str + '\n\n')
        fout.write(num_uniquely_aligned_str + '\n')
        fout.write(num_multiply_aligned_str + '\n')
        fout.write(num_misassembled_by_blat_str + '\n\n')
        fout.write(num_unaligned_str + '\n\n')

        fout.write(' == ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS ==\n')
        # fout.write(genome_fraction_str + '\n')
        # fout.write(duplication_ratio_wo_database_str + '\n\n')
        fout.write(avg_aligned_fraction_str + '\n')
        fout.write(avg_alignment_len_str + '\n\n')
        fout.write(avg_blocks_per_alignment + '\n')
        fout.write(avg_block_len_str + '\n\n')
        fout.write(avg_mismatches_per_transcript_str + '\n\n')
        fout.write(na50_str + '\n\n')
        # fout.write(avg_gaps_per_alignment_str + '\n')
        # fout.write(avg_gap_length_str + '\n')

        fout.close()

        logger.info('      saved to {}'.format(self.path_txt_alignment))


    # Alignment metrics for misassembled (chimeric) transcripts
    # Misassembly candidates reported by BLAT: transcripts are aligned to the reference genome with BLAT and discordant partial alignments are selected as misassembly candidates.
    # Misassembly candidates reported by BLASTN: transcripts are aligned to the isoform sequences (extracted from the genome using gene database) with BLASTN and discordant partial alignments are selected as misassembly candidate.
    # Misassemblies: misassembly candidates confirmed by both methods described above. Using both methods simultaneously allows to avoid considering misalignments that can be caused, for  example, by paralogous genes.
    def get_misassemblies_report(self, transcripts_metrics, logger):
        logger.info('    Getting ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS report...')

        name_str = '{:<50}'.format('METRICS/TRANSCRIPTS')

        num_transcripts_str = '{:<50}'.format('Transcripts')

        num_misassembled_by_blat_str = '{:<50}'.format('Misassembly candidates reported by BLAT')
        num_misassembled_by_blast_str = '{:<50}'.format('Misassembly candidates reported by BLASTN')
        num_misassembled_together_str = '{:<50}'.format('Misassemblies')

        for i_transcripts in range(len(transcripts_metrics)):
            basic_metrics = transcripts_metrics[i_transcripts].basic_metrics

            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics

            name_str += '{:<25}'.format(transcripts_metrics[i_transcripts].label)
            num_transcripts_str += '{:<25}'.format(basic_metrics.number)

            num_misassembled_by_blat_str += '{:<25}'.format(simple_metrics.num_misassembled_by_blat)
            num_misassembled_by_blast_str += '{:<25}'.format(simple_metrics.num_misassembled_by_blast)
            num_misassembled_together_str += '{:<25}'.format(simple_metrics.num_misassembled_together)

        fout = open(self.path_txt_misassemblies, 'w')

        fout.write(' == ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS (calculated with reference genome but without gene database) ==\n')
        fout.write(name_str + '\n')
        fout.write(num_transcripts_str + '\n\n')

        fout.write(num_misassembled_by_blat_str + '\n')
        fout.write(num_misassembled_by_blast_str + '\n')
        fout.write(num_misassembled_together_str + '\n\n')

        fout.close()

        logger.info('      saved to {}'.format(self.path_txt_misassemblies))


    # Assembly completeness (sensitivity). For the following metrics rnaQUAST attempts to select best-matching isoforms for every transcript. Note that a single transcript can contribute to multiple isoforms in the case of, for example, paralogous genes or genomic repeats.
    # Database coverage: the total number of covered bases in all isoforms divided by the total length of all isoforms.
    # Duplication ratio: total number of aligned bases in assembled transcripts divided by the total number of isoform covered bases. Since rnaQUAST selects the best-matching isoform for every transcript, shared exons covered by a single transcript are not counted more than once.
    # Average number of transcripts mapped to one isoform.
    # 30%-assembled: number of isoforms from the database that have at least 30% captured by a single assembled transcript.
    # 90%-assembled: same as above, but for 90%.
    # 30%-covered: number of isoforms from the database that have at least 30% of bases covered by all alignments.
    # 90%-covered: same as above, but for 90%.
    # 30%-assembled exons: number of exons from the database that have at least 30%  captured by a single assembled transcript.
    # 90%-assembled exons: same as above, but for 90%.
    # Mean isoform coverage: average x value for all isoforms (with > 0 bases covered)\, where x is the number of isoform bases covered by all assembled transcripts divided by the isoform length.
    # Mean isoform assembly: average x value for all isoforms with > 0 bases covered, where x is the number of isoform bases captured by a single assembled transcript divided by the isoform length.
    # Mean exon coverage: average x value for all exons with > 0 bases covered, where x is the number of exon bases covered by all assembled transcripts divided by the exon length.
    # Average percentage of isoform 30%-covered exons. For each isoform rnaQUAST calculates the number of 30%-covered exons divided by the total number of exons. Afterwards it computes average value for all covered isoforms.
    # Average percentage of isoform 90%-covered exons: same as above but for 90%-covered exons.

    # CEGMA metrics (http://korflab.ucdavis.edu/Datasets/cegma/). The following metrics are calculated only when --cegma option is used (see options for details).
    # Complete: predicted proteins in the set of 248 CEGs that when aligned to the HMM for the KOG for that protein-family, give an alignment length that is 70% of the protein length.
    # Partial:  the rest of predicted proteins in the set of 248 CEGs that have an alignment.
    def get_sensitivity_report(self, transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION):
        logger.info('    Getting ASSEMBLY COMPLETENESS (SENSITIVITY) report...')

        name_str = '{:<50}'.format('METRICS/TRANSCRIPTS')

        database_coverage_str = '{:<50}'.format('Database coverage')

        duplication_str = '{:<50}'.format('Duplication ratio')
        avg_num_transcripts_in_isoform_str = '{:<50}'.format('Avg. number of transcripts mapped to one isoform')

        assembled_well_genes_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled genes')
        assembled_fully_genes_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled genes')
        covered_well_genes_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered genes')
        covered_fully_genes_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered genes')

        assembled_well_isoforms_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled isoforms')
        assembled_fully_isoforms_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled isoforms')
        covered_well_isoforms_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered isoforms')
        covered_fully_isoforms_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered isoforms')

        assembled_exons_well_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-assembled exons')
        assembled_exons_fully_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled exons')
        mean_isoform_assembly_str = '{:<50}'.format('Mean isoform assembly')
        mean_isoform_cov_str = '{:<50}'.format('Mean isoform coverage')
        mean_exon_cov_str = '{:<50}'.format('Mean exon coverage')
        isoform_well_cov_exons_str = '{:<50}'.format('Avg. percentage of isoform ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-covered exons')
        isoform_fully_cov_exons_str = '{:<50}'.format('Avg. percentage of isoform ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold * 100)) + '%-covered exons')

        # cegma_complete_str = '{:<50}'.format('Complete')
        # cegma_partial_str = '{:<50}'.format('Partial')

        busco_complete_str = '{:<50}'.format('Complete')
        busco_partial_str = '{:<50}'.format('Partial')

        GeneMarkS_T_genes_str = '{:<50}'.format('Genes')

        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].assembly_completeness_metrics is not None:
                name_str += '{:<25}'.format(transcripts_metrics[i_transcripts].label)

                isoforms_coverage = transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage
                if isoforms_coverage is not None:
                    database_coverage_str += '{:<25}'.format(round(isoforms_coverage.fraction_annotation_mapped, PRECISION))

                    duplication_str += '{:<25}'.format(round(isoforms_coverage.avg_duplication_ratio, PRECISION))
                    avg_num_transcripts_in_isoform_str += '{:<25}'.format(round(isoforms_coverage.avg_num_transcripts_mapped_to_isoform, PRECISION))

                    assembled_well_genes_str += '{:<25}'.format(isoforms_coverage.num_well_assembled_genes)
                    assembled_fully_genes_str += '{:<25}'.format(isoforms_coverage.num_fully_assembled_genes)
                    covered_well_genes_str += '{:<25}'.format(isoforms_coverage.num_well_covered_genes)
                    covered_fully_genes_str += '{:<25}'.format(isoforms_coverage.num_fully_covered_genes)

                    assembled_well_isoforms_str += '{:<25}'.format(isoforms_coverage.num_well_assembled_isoforms)
                    assembled_fully_isoforms_str += '{:<25}'.format(isoforms_coverage.num_fully_assembled_isoforms)
                    covered_well_isoforms_str += '{:<25}'.format(isoforms_coverage.num_well_covered_isoforms)
                    covered_fully_isoforms_str += '{:<25}'.format(isoforms_coverage.num_fully_covered_isoforms)

                    assembled_exons_well_str += '{:<25}'.format(isoforms_coverage.num_well_assembled_exons)
                    assembled_exons_fully_str += '{:<25}'.format(isoforms_coverage.num_fully_assembled_exons)
                    mean_isoform_cov_str += '{:<25}'.format(round(isoforms_coverage.avg_covered_fraction, PRECISION))
                    mean_isoform_assembly_str += '{:<25}'.format(round(isoforms_coverage.avg_assembled_fraction, PRECISION))
                    mean_exon_cov_str += '{:<25}'.format(round(isoforms_coverage.avg_covered_fraction_exons, PRECISION))
                    isoform_well_cov_exons_str += '{:<25}'.format(round(isoforms_coverage.avg_percentage_isoform_well_covered_exons, PRECISION))
                    isoform_fully_cov_exons_str += '{:<25}'.format(round(isoforms_coverage.avg_percentage_isoform_fully_covered_exons, PRECISION))

            # cegma_metrics = transcripts_metrics[i_transcripts].cegma_metrics
            # if cegma_metrics is not None:
            #     cegma_complete_str += '{:<25}'.format(cegma_metrics.complete_completeness)
            #     cegma_partial_str += '{:<25}'.format(cegma_metrics.partial_completeness)

            busco_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.busco_metrics
            if busco_metrics is not None:
                busco_complete_str += '{:<25}'.format(round(busco_metrics.complete_completeness, PRECISION))
                busco_partial_str += '{:<25}'.format(round(busco_metrics.partial_completeness, PRECISION))

            geneMarkS_T_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.geneMarkS_T_metrics
            if geneMarkS_T_metrics is not None:
                GeneMarkS_T_genes_str += '{:<25}'.format(geneMarkS_T_metrics.genes)

        fout = open(self.path_txt_sensitivity, 'w')

        fout.write(' == ASSEMBLY COMPLETENESS (SENSITIVITY) ==\n')

        fout.write(name_str + '\n')

        if transcripts_metrics[0].assembly_completeness_metrics.isoforms_coverage is not None:
            fout.write(database_coverage_str + '\n')

            fout.write(duplication_str + '\n')
            fout.write(avg_num_transcripts_in_isoform_str + '\n\n')

            fout.write(assembled_well_genes_str + '\n')
            fout.write(assembled_fully_genes_str + '\n')
            fout.write(covered_well_genes_str + '\n')
            fout.write(covered_fully_genes_str + '\n\n')

            fout.write(assembled_well_isoforms_str + '\n')
            fout.write(assembled_fully_isoforms_str + '\n')
            fout.write(covered_well_isoforms_str + '\n')
            fout.write(covered_fully_isoforms_str + '\n\n')

            fout.write(assembled_exons_well_str + '\n')
            fout.write(assembled_exons_fully_str + '\n\n')

            fout.write(mean_isoform_assembly_str + '\n\n')
            fout.write(mean_isoform_cov_str + '\n')

            fout.write(mean_exon_cov_str + '\n\n')

            fout.write(isoform_well_cov_exons_str + '\n')
            fout.write(isoform_fully_cov_exons_str + '\n\n')

        # if transcripts_metrics[0].assembly_completeness_metrics.cegma_metrics is not None:
        #     fout.write(' == CEGMA METRICS ==\n')
        #     fout.write(cegma_complete_str + '\n')
        #     fout.write(cegma_partial_str + '\n\n')

        if transcripts_metrics[0].assembly_completeness_metrics.busco_metrics is not None:
            fout.write(' == BUSCO METRICS ==\n')
            fout.write(busco_complete_str + '\n')
            fout.write(busco_partial_str + '\n\n')

        if transcripts_metrics[0].assembly_completeness_metrics.geneMarkS_T_metrics is not None:
            fout.write(' == GeneMarkS-T METRICS ==\n')
            fout.write(GeneMarkS_T_genes_str + '\n\n')

        fout.close()

        logger.info('      saved to {}'.format(self.path_txt_sensitivity))


        # fout.write('{:<100}'.format('Assembled isoforms') + str(isoforms_coverage.num_assembled) + '\n\n')

        # fout.write('{:<100}'.format('Assembled exons') + str(exons_coverage.num_assembled) + '\n\n')

        # fout.write('COVERAGE BY EACH MAPPED TRANSCRIPT SEPARATELY:\n')
        # fout.write('DISTRIBUTIONS FOR ALIGNED TRANSCRIPTS:\n')

        # fout.write('{:<100}'.format('Avg. fraction of exon assembled in isoform') + str(round(isoforms_coverage.avg_assembled_fraction_exon, precision)) + '\n\n')

        # fout.write('{:<100}'.format('Avg. percentage of isoform well-assembled exons') + str(round(isoforms_coverage.avg_percentage_well_assembled_exons, precision)) + '\n')
        # fout.write('{:<100}'.format('Avg. percentage of isoform fully-assembled exons') + str(round(isoforms_coverage.avg_percentage_fully_assembled_exons, precision)) + '\n\n')

        # fout.write('DISTRIBUTIONS FOR MAPPED ISOFORMS:\n')
        # fout.write('{:<100}'.format('Avg. fraction of isoform assembled') + str(round(isoforms_coverage.avg_over_isoforms_assembled_fraction, precision)) + '\n\n')

        # fout.write('{:<100}'.format('Avg. fraction of exon assembled in isoform') + str(round(isoforms_coverage.avg_over_isoforms_avg_assembled_fraction_exon, precision)) + '\n\n')

        # fout.write('{:<100}'.format('Avg. percentage of isoform well-assembled exons') + str(round(isoforms_coverage.avg_over_isoforms_percentage_well_assembled_exons, precision)) + '\n')
        # fout.write('{:<100}'.format('Avg. percentage of isoform fully-assembled exons') + str(round(isoforms_coverage.avg_over_isoforms_percentage_fully_assembled_exons, precision)) + '\n\n\n')


        # fout.write('DISTRIBUTIONS FOR MAPPED EXONS:\n')
        # fout.write('{:<100}'.format('Avg. fraction of exon assembled') + str(round(exons_coverage.avg_over_exons_assembled_fraction, precision)) + '\n\n\n')

        # fout.write('FOR ALL ALIGNED TRANSCRIPTS:\n')
        # fout.write('{:<100}'.format('Well-assembled in average isoforms') + str(isoforms_coverage.num_well_assembled_in_average) + '\n')
        # fout.write('{:<100}'.format('Fully-assembled in average isoforms') + str(isoforms_coverage.num_fully_assembled_in_average) + '\n\n')

        # fout.write('{:<100}'.format('Well-assembled in average exons') + str(exons_coverage.num_well_assembled_in_average) + '\n')
        # fout.write('{:<100}'.format('Fully-assembled in average exons') + str(exons_coverage.num_fully_assembled_in_average) + '\n\n')

        # fout.write('{:<100}'.format('Avg. fraction of exon assembled') + str(round(exons_coverage.avg_assembled_fraction, precision)) + '\n\n\n')

        # fout.write('COVERAGE BY ALL MAPPED TRANSCRIPTS AT ONCE:\n')
        # fout.write('DISTRIBUTIONS FOR ALIGNED TRANSCRIPTS:\n')

        # fout.write('{:<100}'.format('Avg. fraction of exon covered in isoform') + str(round(isoforms_coverage.avg_covered_fraction_exon, precision)) + '\n\n')

        # fout.write('FOR ALL ALIGNED TRANSCRIPTS:\n')
        # fout.write('{:<100}'.format('Well-covered exons') + str(exons_coverage.num_well_covered) + '\n')
        # fout.write('{:<100}'.format('Fully-covered exons') + str(exons_coverage.num_fully_covered) + '\n\n')


    # Assembly specificity. To compute the following metrics we use only transcripts that have at least one significant alignment and are not misassembled.
    # Unannotated: total number of transcripts that do not cover any isoform from the database.
    # 30%-matched: total number of transcripts that have at least 30% covering an isoform from the database.
    # 90%-matched: same as above, but for 90%.
    # Mean transcript match: average x value for all transcripts (with > 0 bases matched), where x is the number of transcript bases covering an isoform from the database divided by the transcript length.
    # Mean block match: average x value for all blocks (with > 0 bases matched), where x is the number of block bases covering an isoform from the database divided by the block length.
    # 30%-matched blocks: total number of blocks that have at least 30% covering an isoform from the database.
    # 90%-matched blocks: same as above, but for 90%.
    def get_specificity_report(self, transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, precision):
        logger.info('    Getting ASSEMBLY SPECIFICITY report...')

        name_str = '{:<50}'.format('METRICS/TRANSCRIPTS')

        unannotated_str = '{:<50}'.format('Unannotated')
        matched_well_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_transcript_threshold * 100)) + '%-matched')
        matched_fully_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_transcript_threshold * 100)) + '%-matched')
        mean_transcript_match_str = '{:<50}'.format('Mean fraction of transcript matched')
        mean_block_match_str = '{:<50}'.format('Mean fraction of block matched')
        matched_blocks_well_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_block_threshold * 100)) + '%-matched blocks')
        matched_blocks_fully_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_block_threshold * 100)) + '%-matched blocks')
        matched_len_str = '{:<50}'.format('Matched length')
        unmatched_len_str = '{:<50}'.format('Unmatched length')

        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_coverage = transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage

            name_str += '{:<25}'.format(transcripts_metrics[i_transcripts].label)

            if transcripts_coverage is not None:
                unannotated_str += '{:<25}'.format(transcripts_coverage.num_unannotated_transcripts)
                matched_well_str += '{:<25}'.format(transcripts_coverage.num_well_covered_transcripts)
                matched_fully_str += '{:<25}'.format(transcripts_coverage.num_fully_covered_transcripts)
                mean_transcript_match_str += '{:<25}'.format(round(transcripts_coverage.avg_covered_fraction_whole_transcript, precision))
                mean_block_match_str += '{:<25}'.format(round(transcripts_coverage.avg_covered_fraction_block, precision))
                matched_blocks_well_str += '{:<25}'.format(round(transcripts_coverage.percentage_well_covered_blocks, precision))
                matched_blocks_fully_str += '{:<25}'.format(round(transcripts_coverage.percentage_fully_covered_blocks, precision))
                matched_len_str += '{:<25}'.format(transcripts_coverage.matched_len)
                unmatched_len_str += '{:<25}'.format(transcripts_coverage.unmatched_len)

        fout = open(self.path_txt_specificity, 'w')


        if transcripts_metrics[0].assembly_correctness_metrics.transcripts_coverage:
            fout.write(' == ASSEMBLY SPECIFICITY ==\n')

            fout.write(name_str + '\n')

            fout.write(unannotated_str + '\n\n')

            fout.write(matched_well_str + '\n')
            fout.write(matched_fully_str + '\n\n')

            fout.write(mean_transcript_match_str + '\n\n')

            fout.write(mean_block_match_str + '\n\n')

            fout.write(matched_blocks_well_str + '\n')
            fout.write(matched_blocks_fully_str + '\n\n')

            fout.write(matched_len_str + '\n')
            fout.write(unmatched_len_str + '\n\n')

        fout.close()

        logger.info('      saved to {}'.format(self.path_txt_specificity))

        # FOR ALL TRANSCRIPTS ALIGNMENTS:
        # fout.write('FOR ALL TRANSCRIPTS ALIGNMENTS:\n')

        # fout.write('{:<100}'.format('Annotated transcripts') + str(transcripts_coverage.num_annotated_transcripts) + '\n')

        # fout.write('DISTRIBUTIONS FOR TRANSCRIPTS ALIGNMENTS:\n')
        # fout.write('{:<100}'.format('Avg. fraction of aligned part annotated') + str(round(transcripts_coverage.avg_covered_fraction_aligned_part, precision)) + '\n\n')

        # fout.write('{:<100}'.format('Avg. fraction of block annotated in transcript:') + str(round(transcripts_coverage.avg_covered_fraction_block_in_t, precision)) + '\n\n')

        # fout.write('{:<100}'.format('Avg. percentage of well-annotated blocks in transcript') + str(round(transcripts_coverage.avg_percentage_well_covered_blocks_in_t, precision)) + '\n')
        # fout.write('{:<100}'.format('Avg. percentage of fully-annotated blocks in transcript') + str(round(transcripts_coverage.avg_percentage_fully_covered_blocks_in_t, precision)) + '\n\n\n')

        # fout.write('FOR ALL TRANSCRIPTS MAPPED TO SPECIFIC ISOFORM:\n')
        # fout.write('{:<100}'.format('Avg. fraction of whole transcript annotated') + str(round(transcripts_coverage.avg_over_isoforms_covered_fraction_whole_transcript, precision)) + '\n')
        # fout.write('{:<100}'.format('Avg. fraction of aligned part annotated') + str(round(transcripts_coverage.avg_over_isoforms_covered_fraction_aligned_part, precision)) + '\n\n')

        # fout.write('{:<100}'.format('Avg. fraction of total transcript annotated length over total whole transcript lengths') + str(round(transcripts_coverage.avg_covered_fraction_whole_transcript_over_tot_len, precision)) + '\n')
        # fout.write('{:<100}'.format('Fraction of total transcript annotated length over total aligned parts') + str(round(transcripts_coverage.avg_covered_fraction_aligned_part_over_tot_len, precision)) + '\n\n\n')