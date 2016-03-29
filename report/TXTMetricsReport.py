__author__ = 'letovesnoi'

import os

from general import rqconfig


class TXTMetricsReport():
    """Class which generate txt reports"""

    def __init__(self, args_blast, txt_report_dir, labels, transcripts_metrics, db_genes_metrics, reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS):
        logger.print_timestamp('  ')
        logger.info('  Getting TXT report...')

        self.txt_reports_dir = txt_report_dir

        self.widths = TXTMetricsReport.get_column_widths(labels)

        if db_genes_metrics is not None:
            self.path_txt_database_metrics = os.path.join(self.txt_reports_dir, 'database_metrics.txt')
            self.get_db_genes_metrics_report(db_genes_metrics, logger, PRECISION)

        if reads_coverage is not None:
            self.path_txt_reads_coverage_metrics = self.path_txt_database_metrics
            self.get_reads_coverage_report(reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

            self.path_txt_relative_database_coverage = os.path.join(self.txt_reports_dir, 'relative_database_coverage.txt')
            self.get_relative_database_coverage_report(transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

        self.path_txt_basic = os.path.join(self.txt_reports_dir, 'basic_metrics.txt')
        self.get_basic_metrics_report(transcripts_metrics, logger, TRANSCRIPT_LENS, PRECISION)

        self.path_txt_alignment = os.path.join(self.txt_reports_dir, 'alignment_metrics.txt')
        self.get_alignment_metrics_report(transcripts_metrics, logger, PRECISION)

        self.path_txt_misassemblies = os.path.join(self.txt_reports_dir, 'misassemblies.txt')
        self.get_misassemblies_report(args_blast, transcripts_metrics, logger)

        self.path_txt_sensitivity = os.path.join(self.txt_reports_dir, 'sensitivity.txt')
        self.get_sensitivity_report(transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

        self.path_txt_specificity = os.path.join(self.txt_reports_dir, 'specificity.txt')
        self.get_specificity_report(transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

        logger.info('  Done.')


    @classmethod
    def get_column_widths(cls, labels):
        widths = []

        widths.append(rqconfig.space_label)
        for t_label in labels:
            curr_width = len(t_label) + 2
            widths.append(max(curr_width, rqconfig.space_value))

        return widths


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

        with open(self.path_txt_database_metrics, 'w') as out_handle:
            out_handle.write(' == GENE DATABASE METRICS ==\n')

            column_width_str = '{:<' + str(self.widths[0]) + '}'

            out_handle.write(column_width_str.format('Genes') + str(db_genes_metrics.genes_num) + '\n')
            # gff source field unfortunately contain database name instead of transcripts type protein_coding:
            if db_genes_metrics.protein_coding_genes_num != 0:
                out_handle.write(column_width_str.format('Protein coding genes') + str(db_genes_metrics.protein_coding_genes_num) + '\n\n')

            out_handle.write(column_width_str.format('Isoforms') + str(db_genes_metrics.isoforms_num) + '\n')
            # gff source field unfortunately contain database name instead of transcripts type protein_coding:
            if db_genes_metrics.protein_coding_isoforms_num != 0:
                out_handle.write(column_width_str.format('Protein coding isoforms') + str(db_genes_metrics.protein_coding_isoforms_num) + '\n\n')

            out_handle.write(column_width_str.format('Exons') + str(db_genes_metrics.tot_exons_num) + '\n\n')

            out_handle.write(column_width_str.format('Introns') + str(db_genes_metrics.tot_introns_num) + '\n\n')

            out_handle.write(column_width_str.format('Avg. length of all isoforms') + str(round(db_genes_metrics.avg_isoform_len, PRECISION)) + '\n\n')
            out_handle.write(column_width_str.format('Total length of all isoforms') + str(db_genes_metrics.tot_isoforms_len) + '\n')

            out_handle.write(column_width_str.format('Avg. exon length') + str(round(db_genes_metrics.avg_exon_len, PRECISION)) + '\n\n')

            out_handle.write(column_width_str.format('Avg. intron length') + str(round(db_genes_metrics.avg_intron_len, PRECISION)) + '\n\n')

            out_handle.write(column_width_str.format('Avg. number of exons per isoform') + str(round(db_genes_metrics.avg_exons_num, PRECISION)) + '\n')
            out_handle.write(column_width_str.format('Max. number of exons per isoform') + str(db_genes_metrics.max_exons_num) + '\n\n')

            # out_handle.write('{:<100}'.format('Avg. number of introns per isoform') + str(round(basic_isoforms_metrics.avg_introns_num, precision)) + '\n\n')

            # out_handle.write('ISOFORMS:\n')
            # out_handle.write('{:<100}'.format('Total length with introns') + str(basic_isoforms_metrics.tot_len_w_introns) + '\n')
            # out_handle.write('{:<100}'.format('Min. length with introns') + str(basic_isoforms_metrics.min_len_w_introns) + '\n')
            # out_handle.write('{:<100}'.format('Max. length with introns') + str(basic_isoforms_metrics.max_len_w_introns) + '\n')
            # out_handle.write('{:<100}'.format('Avg. length with introns') + str(round(basic_isoforms_metrics.avg_len_w_introns, precision)) + '\n\n')

            # out_handle.write('{:<100}'.format('Min. isoform length') + str(basic_isoforms_metrics.min_len_wout_introns) + '\n')
            # out_handle.write('{:<100}'.format('Max. isoform length') + str(basic_isoforms_metrics.max_len_wout_introns) + '\n')

            # out_handle.write('EXONS:\n')
            # out_handle.write('{:<100}'.format('Min. exons per isoform') + str(basic_isoforms_metrics.min_exons_num) + '\n')

            # out_handle.write('{:<100}'.format('Min. exon length') + str(basic_isoforms_metrics.min_exon_len) + '\n')
            # out_handle.write('{:<100}'.format('Max. exon length') + str(basic_isoforms_metrics.max_exon_len) + '\n')

            # out_handle.write('INTRONS:\n')
            # out_handle.write('{:<100}'.format('Min. introns per isoform') + str(basic_isoforms_metrics.min_introns_num) + '\n')
            # out_handle.write('{:<100}'.format('Max. introns per isoform') + str(basic_isoforms_metrics.max_introns_num) + '\n')

            # out_handle.write('{:<100}'.format('Min. intron length') + str(basic_isoforms_metrics.min_intron_len) + '\n')
            # out_handle.write('{:<100}'.format('Max. intron length') + str(basic_isoforms_metrics.max_intron_len) + '\n')
            # out_handle.write('{:<100}'.format('Total introns length [bp]') + str(basic_isoforms_metrics.tot_introns_len) + '\n')

            logger.info('      saved to {}'.format(self.path_txt_database_metrics))


    def get_reads_coverage_report(self, reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION):
        logger.info('    Getting COVERAGE BY READS report...')

        with open(self.path_txt_reads_coverage_metrics, 'a') as out_handle:
            out_handle.write(' == COVERAGE BY READS ==\n')

            column_width_str = '{:<' + str(self.widths[0]) + '}'

            out_handle.write(column_width_str.format('Database coverage') + str(round(reads_coverage.fraction_annotation_mapped_by_reads, PRECISION)) + '\n\n')

            out_handle.write(column_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered genes') + str(reads_coverage.num_well_expressed_genes) + '\n')
            out_handle.write(column_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered genes') + str(reads_coverage.num_fully_expressed_genes) + '\n\n')

            out_handle.write(column_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered isoforms') + str(reads_coverage.num_well_expressed_isoforms) + '\n')
            out_handle.write(column_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered isoforms') + str(reads_coverage.num_fully_expressed_isoforms) + '\n\n')

            out_handle.write(column_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-covered exons') + str(reads_coverage.num_well_expressed_exons) + '\n')
            out_handle.write(column_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold * 100)) + '%-covered exons') + str(reads_coverage.num_fully_expressed_exons) + '\n\n')

        logger.info('      saved to {}'.format(self.path_txt_reads_coverage_metrics))


    def get_relative_database_coverage_report(self, transcripts_metrics, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION):
        logger.info('    Getting RELATIVE DATABASE COVERAGE report...')

        label_width_str = '{:<' + str(self.widths[0]) + '}'

        name_str = label_width_str.format('METRICS/TRANSCRIPTS')

        database_coverage_str = label_width_str.format('Relative database coverage')

        well_assembled_genes_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled genes')
        fully_assembled_genes_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled genes')

        well_covered_genes_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered genes')
        fully_covered_genes_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered genes')


        well_assembled_isoforms_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled isoforms')
        fully_assembled_isoforms_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled isoforms')

        well_covered_isoforms_str = label_width_str.format(str('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100))) + '%-covered isoforms')
        fully_covered_isoforms_str = label_width_str.format(str('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100))) + '%-covered isoforms')


        well_assembled_exons_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-assembled exons')
        fully_assembled_exons_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold * 100)) + '%-assembled exons')

        well_covered_exons_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-covered exons')
        fully_covered_exons_str = label_width_str.format('Relative ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold * 100)) + '%-covered exons')

        for i_transcripts in range(len(transcripts_metrics)):
            value_width_str = '{:<' + str(self.widths[i_transcripts + 1]) + '}'

            if transcripts_metrics[i_transcripts].assembly_completeness_metrics is not None and \
                            transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage is not None:
                relative_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.relative_database_coverage

                if relative_metrics is not None:
                    database_coverage_str += value_width_str.format(round(relative_metrics.database_coverage, PRECISION))

                    well_assembled_genes_str += value_width_str.format(round(relative_metrics.percentage_well_assembled_genes, PRECISION))
                    fully_assembled_genes_str += value_width_str.format(round(relative_metrics.percentage_fully_assembled_genes, PRECISION))

                    well_covered_genes_str += value_width_str.format(round(relative_metrics.percentage_well_covered_genes, PRECISION))
                    fully_covered_genes_str += value_width_str.format(round(relative_metrics.percentage_fully_covered_genes, PRECISION))


                    well_assembled_isoforms_str += value_width_str.format(round(relative_metrics.percentage_well_assembled_isoforms, PRECISION))
                    fully_assembled_isoforms_str += value_width_str.format(round(relative_metrics.percentage_fully_assembled_isoforms, PRECISION))

                    well_covered_isoforms_str += value_width_str.format(round(relative_metrics.percentage_well_covered_isoforms, PRECISION))
                    fully_covered_isoforms_str += value_width_str.format(round(relative_metrics.percentage_fully_covered_isoforms, PRECISION))


                    well_assembled_exons_str += value_width_str.format(round(relative_metrics.percentage_well_assembled_exons, PRECISION))
                    fully_assembled_exons_str += value_width_str.format(round(relative_metrics.percentage_fully_assembled_exons, PRECISION))

                    well_covered_exons_str += value_width_str.format(round(relative_metrics.percentage_well_covered_exons, PRECISION))
                    fully_covered_exons_str += value_width_str.format(round(relative_metrics.percentage_fully_covered_exons, PRECISION))
                else:
                    database_coverage_str += value_width_str.format('*')

                    well_assembled_genes_str += value_width_str.format('*')
                    fully_assembled_genes_str += value_width_str.format('*')

                    well_covered_genes_str += value_width_str.format('*')
                    fully_covered_genes_str += value_width_str.format('*')


                    well_assembled_isoforms_str += value_width_str.format('*')
                    fully_assembled_isoforms_str += value_width_str.format('*')

                    well_covered_isoforms_str += value_width_str.format('*')
                    fully_covered_isoforms_str += value_width_str.format('*')


                    well_assembled_exons_str += value_width_str.format('*')
                    fully_assembled_exons_str += value_width_str.format('*')

                    well_covered_exons_str += value_width_str.format('*')
                    fully_covered_exons_str += value_width_str.format('*')

        with open(self.path_txt_relative_database_coverage, 'w') as out_handle:
            out_handle.write(name_str + '\n')

            self.write_metric_str(database_coverage_str + '\n\n', out_handle, len(transcripts_metrics), first_type='RELATIVE DATABASE COVERAGE')

            self.write_metric_str(well_assembled_genes_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(fully_assembled_genes_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(well_covered_genes_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(fully_covered_genes_str + '\n\n', out_handle, len(transcripts_metrics))

            self.write_metric_str(well_assembled_isoforms_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(fully_assembled_isoforms_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(well_covered_isoforms_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(fully_covered_isoforms_str + '\n\n', out_handle, len(transcripts_metrics))

            self.write_metric_str(well_assembled_exons_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(fully_assembled_exons_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(well_covered_exons_str + '\n', out_handle, len(transcripts_metrics))
            self.write_metric_str(fully_covered_exons_str + '\n\n', out_handle, len(transcripts_metrics))

        logger.info('      saved to {}'.format(self.path_txt_relative_database_coverage))


    def write_metric_str(self, metric_str, out_handle, num_transcripts, first_type=None):
        num_values = 0
        for val in metric_str.strip().split('  '):
            if val != '':
                num_values += 1
        if metric_str.count('*') != num_transcripts and num_values == num_transcripts + 1:
            if first_type is not None:
                out_handle.write(' == ' + first_type + ' == \n')
            out_handle.write(metric_str)


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

        label_width_str = '{:<' + str(self.widths[0]) + '}'

        name_str = label_width_str.format('METRICS/TRANSCRIPTS')

        num_transcripts_str = label_width_str.format('Transcripts')
        num_transcripts_500_str = label_width_str.format('Transcripts > {} bp'.format(str(TRANSCRIPT_LENS[0])))
        num_transcripts_1000_str = label_width_str.format('Transcripts > {} bp'.format(str(TRANSCRIPT_LENS[1])))
        avg_len_transcripts_str = label_width_str.format('Average length of assembled transcripts')
        longest_transcript_str = label_width_str.format('Longest transcript')
        tot_len_str = label_width_str.format('Total length')
        n50_str = label_width_str.format('Transcript N50')

        for i_transcripts in range(len(transcripts_metrics)):
            value_width_str = '{:<' + str(self.widths[i_transcripts + 1]) + '}'

            basic_metrics = transcripts_metrics[i_transcripts].basic_metrics

            if basic_metrics is not None:
                name_str += value_width_str.format(transcripts_metrics[i_transcripts].label)
                num_transcripts_str += value_width_str.format(basic_metrics.number)

                num_transcripts_500_str += value_width_str.format(basic_metrics.num_transcripts_500)
                num_transcripts_1000_str += value_width_str.format(basic_metrics.num_transcripts_1000)

                avg_len_transcripts_str += value_width_str.format(round(basic_metrics.avg_len, PRECISION))
                longest_transcript_str += value_width_str.format(basic_metrics.max_len)
                tot_len_str += value_width_str.format(basic_metrics.tot_len)

                n50_str += value_width_str.format(basic_metrics.n50)

            # out_handle.write('{:<100}'.format('Min. length') + str(self.min_len) + '\n')
            # if args.annotation is not None:
                # number with length inside isoform length range:
                # out_handle.write('{:<100}'.format('Outliers') + str(self.num_outliers_len) + '\n')
            # out_handle.write('\n')

        out_handle = open(self.path_txt_basic, 'w')

        out_handle.write(name_str + '\n')
        self.write_metric_str(num_transcripts_str + '\n\n', out_handle, len(transcripts_metrics), first_type='BASIC TRANSCRIPTS METRICS (calculated without reference genome and gene database)')

        self.write_metric_str(num_transcripts_500_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(num_transcripts_1000_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(avg_len_transcripts_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(longest_transcript_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(tot_len_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(n50_str + '\n', out_handle, len(transcripts_metrics))

        out_handle.close()

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

        label_width_str = '{:<' + str(self.widths[0]) + '}'

        name_str = label_width_str.format('METRICS/TRANSCRIPTS')

        num_transcripts_str = label_width_str.format('Transcripts')
        num_aligned_str = label_width_str.format('Aligned')
        num_uniquely_aligned_str = label_width_str.format('Uniquely aligned')
        num_multiply_aligned_str = label_width_str.format('Multiply aligned')
        num_misassembled_by_blat_str = label_width_str.format('Misassembly candidates reported by GMAP (or BLAT)')
        num_unaligned_str = label_width_str.format('Unaligned')

        # Alignment metrics for uniquely aligned transcripts
        # genome_fraction_str = '{:<50}'.format('Genome fraction')
        # duplication_ratio_wo_database_str = '{:<50}'.format('Duplication ratio w/o database')
        avg_aligned_fraction_str = label_width_str.format('Average aligned fraction')
        avg_alignment_len_str = label_width_str.format('Average alignment length')
        avg_blocks_per_alignment = label_width_str.format('Average blocks per alignment')
        avg_block_len_str = label_width_str.format('Average block length')
        avg_mismatches_per_transcript_str = label_width_str.format('Average mismatches per transcript')
        na50_str = label_width_str.format('NA50')
        # avg_gaps_per_alignment_str = '{:<50}'.format('Average gaps per alignment')
        # avg_gap_length_str = '{:<50}'.format('Average gap length')

        for i_transcripts in range(len(transcripts_metrics)):
            value_width_str = '{:<' + str(self.widths[i_transcripts + 1]) + '}'

            basic_metrics = transcripts_metrics[i_transcripts].basic_metrics

            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics

            name_str += value_width_str.format(transcripts_metrics[i_transcripts].label)

            if basic_metrics is not None and simple_metrics is not None:
                num_transcripts_str += value_width_str.format(basic_metrics.number)

            if simple_metrics is not None:
                num_aligned_str += value_width_str.format(simple_metrics.num_aligned)

                num_uniquely_aligned_str += value_width_str.format(simple_metrics.num_unique_aligned)
                num_multiply_aligned_str += value_width_str.format(simple_metrics.num_mul_aligned)
                num_misassembled_by_blat_str += value_width_str.format(simple_metrics.num_misassembled_by_blat)

                num_unaligned_str += value_width_str.format(simple_metrics.num_unaligned)

                # number of alignments in PSL-file:
                # out_handle.write('{:<100}'.format('Alignments') + str(self.num_alignments) + '\n\n')

                # number of not misassembled aligned transcripts (best union is single union):
                # out_handle.write('{:<100}'.format('Not misassembled aligned') + str(self.num_non_misassembled) + '\n\n')

                # TEMPORARY FOR TEST BLAT ALIGNER:
                # out_handle.write('{:<100}'.format('Multiple equal query aligned') + str(simple_metrics.num_mul_equal_query_aligned) + '\n')
                # out_handle.write('{:<100}'.format('Multiple equal query aligned over misassemblies') + str(simple_metrics.num_mul_equal_query_aligned_mis) + '\n\n')

                # genome_fraction_str += '{:<25}'.format(round(simple_metrics.fraction_genome_mapped, precision))
                # duplication_ratio_wo_database_str += '{:<25}'.format(round(simple_metrics.avg_duplication_ratio, precision))

                avg_aligned_fraction_str += value_width_str.format(round(simple_metrics.avg_fraction, PRECISION))
                avg_alignment_len_str += value_width_str.format(round(simple_metrics.avg_alignment_len, PRECISION))

                avg_blocks_per_alignment += value_width_str.format(round(simple_metrics.avg_blocks_num, PRECISION))
                avg_block_len_str += value_width_str.format(round(simple_metrics.avg_block_len, PRECISION))

                avg_mismatches_per_transcript_str += value_width_str.format(round(simple_metrics.avg_mismatch_num, PRECISION))

                na50_str += value_width_str.format(simple_metrics.na50)
            else:
                num_aligned_str += value_width_str.format('*')

                num_uniquely_aligned_str += value_width_str.format('*')
                num_multiply_aligned_str += value_width_str.format('*')
                num_misassembled_by_blat_str += value_width_str.format('*')

                num_unaligned_str += value_width_str.format('*')

                avg_aligned_fraction_str += value_width_str.format('*')
                avg_alignment_len_str += value_width_str.format('*')

                avg_blocks_per_alignment += value_width_str.format('*')
                avg_block_len_str += value_width_str.format('*')

                avg_mismatches_per_transcript_str += value_width_str.format('*')

                na50_str += value_width_str.format('*')

            # avg_gaps_per_alignment_str +=  '{:<25}'.format(round(simple_metrics.avg_qgap_num, precision))
            # avg_gap_length_str += '{:<25}'.format(round(simple_metrics.avg_qgap_len, precision))

            # LENGTH:
            # out_handle.write('{:<100}'.format('Total length') + str(self.tot_len) + '\n')
            # out_handle.write('{:<100}'.format('Min. length') + str(self.min_len) + '\n')
            # out_handle.write('{:<100}'.format('Max. length') + str(self.max_len) + '\n')
            # out_handle.write('{:<100}'.format('Avg. length') + str(round(self.avg_len, precision)) + '\n')
            # if args.annotation is not None:
                # number of aligned transcripts len out isoforms without introns length:
                # out_handle.write('{:<100}'.format('Outliers') + str(self.num_outliers_len) + '\n')
            # out_handle.write('\n')

            # REFERENCE LENGTH (from begin to end of alignment):
            # out_handle.write('{:<100}'.format('Total reference length') + str(self.tot_ref_len) + '\n')
            # out_handle.write('{:<100}'.format('Min. reference length') + str(self.min_ref_len) + '\n')
            # out_handle.write('{:<100}'.format('Max. reference length') + str(self.max_ref_len) + '\n')
            # out_handle.write('{:<100}'.format('Avg. reference length') + str(round(self.avg_ref_len, precision)) + '\n')
            # if args.annotation is not None:
                # number of aligned transcripts with length more then max isoform with introns length or less then min isoform with introns length:
                # out_handle.write('{:<100}'.format('Outliers') + str(self.num_outliers_ref_len) + '\n')
            # out_handle.write('\n')

            # ALIGNMENT LENGTH (summary length of aligned blocks):
            # out_handle.write('{:<100}'.format('Total alignment length') + str(self.tot_alignment_len) + '\n')
            # out_handle.write('{:<100}'.format('Min. alignment length') + str(self.min_alignment_len) + '\n')
            # out_handle.write('{:<100}'.format('Max. alignment length') + str(self.max_alignment_len) + '\n')

            # if args.annotation is not None:
                # number of aligned transcripts with alignment length less then min isoform without introns length or more then max isoform without introns length:
                # out_handle.write('{:<100}'.format('Outliers') + str(self.num_outliers_alignment_len) + '\n')
            # out_handle.write('\n')

            # ALIGNMENT FRACTION:
            # out_handle.write('{:<100}'.format('Min. alignment fraction') + str(round(self.min_fraction, precision)) + '\n')
            # out_handle.write('{:<100}'.format('Max. alignment fraction') + str(round(self.max_fraction, precision)) + '\n')

            # BLOCKS NUMBER:
            # out_handle.write('{:<100}'.format('Total blocks number') + str(self.tot_blocks_num) + '\n')
            # out_handle.write('{:<100}'.format('Min. blocks number') + str(self.min_blocks_num) + '\n')
            # out_handle.write('{:<100}'.format('Max. blocks number') + str(self.max_blocks_num) + '\n')
            # if args.annotation is not None:
                # number and percentage of aligned transcripts with number of blocks more then max number of exons per isoform or less then min:
                # out_handle.write('{:<100}'.format('Outliers') + str(self.num_outliers_blocks_num) + '\n')
            # out_handle.write('\n')

            # BLOCKS LENGTH:
            # out_handle.write('{:<100}'.format('Min. block length') + str(self.min_block_len) + '\n')
            # out_handle.write('{:<100}'.format('Max. block length') + str(self.max_block_len) + '\n')
            # if args.annotation is not None:
                # number and percentage of blocks with length more then max length of exon and less then min length of exon:
                # out_handle.write('{:<100}'.format('Outliers') + str(self.num_outliers_blocks_len) + '\n')
            # out_handle.write('\n')

            # MISMATCHES:
            # out_handle.write('{:<100}'.format('Total mismatch number') + str(self.tot_mismatch_num) + '\n')

            # TARGET GAPS NUMBER:
            # out_handle.write('{:<100}'.format('Total target insert number') + str(self.tot_tgap_num) + '\n')
            # out_handle.write('{:<100}'.format('Min. target insert number') + str(self.min_tgap_num) + '\n')
            # out_handle.write('{:<100}'.format('Max. target insert number') + str(self.max_tgap_num) + '\n')
            # out_handle.write('{:<100}'.format('Avg. target insert number') + str(round(self.avg_tgap_num, precision)) + '\n')
            # if args.annotation is not None:
                # number and percentage of transcripts with number of target gaps more then max number of introns per isoform and less then min:
                # out_handle.write('{:<100}'.format('Outliers') + str(self.num_outliers_tgap_num) + '\n')
            # out_handle.write('\n')

            # TARGET GAPS LENGTH:
            # out_handle.write('{:<100}'.format('Total target insert length') + str(self.tot_tgap_len) + '\n')
            # out_handle.write('{:<100}'.format('Min. target insert length') + str(self.min_tgap_len) + '\n')
            # out_handle.write('{:<100}'.format('Max. target insert length') + str(self.max_tgap_len) + '\n')
            # out_handle.write('{:<100}'.format('Avg. target insert length') + str(round(self.avg_tgap_len, precision)) + '\n')
            # if args.annotation is not None:
                # number and percentage of target gaps with length more then max length of intron and less then min length of intron:
                # out_handle.write('{:<100}'.format('Outliers') + str(self.num_outliers_tgap_len) + '\n')
            # out_handle.write('\n')

            # QUERY GAPS NUMBER:
            # out_handle.write('{:<100}'.format('Total query insert number') + str(self.tot_qgap_num) + '\n')
            # QUERY GAPS LENGTH:
            # out_handle.write('{:<100}'.format('Total query insert length') + str(self.tot_qgap_len) + '\n')

        out_handle = open(self.path_txt_alignment, 'w')

        out_handle.write(name_str + '\n')
        self.write_metric_str(num_transcripts_str + '\n\n', out_handle, len(transcripts_metrics), first_type='ALIGNMENT METRICS (calculated with reference genome but without gene database)')

        self.write_metric_str(num_aligned_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(num_uniquely_aligned_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(num_multiply_aligned_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(num_misassembled_by_blat_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(num_unaligned_str + '\n\n', out_handle, len(transcripts_metrics))

        # out_handle.write(genome_fraction_str + '\n')
        # out_handle.write(duplication_ratio_wo_database_str + '\n\n')
        self.write_metric_str(avg_aligned_fraction_str + '\n', out_handle, len(transcripts_metrics), first_type='ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS')
        self.write_metric_str(avg_alignment_len_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(avg_blocks_per_alignment + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(avg_block_len_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(avg_mismatches_per_transcript_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(na50_str + '\n\n', out_handle, len(transcripts_metrics))

        # out_handle.write(avg_gaps_per_alignment_str + '\n')
        # out_handle.write(avg_gap_length_str + '\n')

        out_handle.close()

        logger.info('      saved to {}'.format(self.path_txt_alignment))


    # Alignment metrics for misassembled (chimeric) transcripts
    # Misassembly candidates reported by BLAT: transcripts are aligned to the reference genome with BLAT and discordant partial alignments are selected as misassembly candidates.
    # Misassembly candidates reported by BLASTN: transcripts are aligned to the isoform sequences (extracted from the genome using gene database) with BLASTN and discordant partial alignments are selected as misassembly candidate.
    # Misassemblies: misassembly candidates confirmed by both methods described above. Using both methods simultaneously allows to avoid considering misalignments that can be caused, for  example, by paralogous genes.
    def get_misassemblies_report(self, args_blast, transcripts_metrics, logger):
        logger.info('    Getting ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS report...')

        label_width_str = '{:<' + str(self.widths[0]) + '}'

        name_str = label_width_str.format('METRICS/TRANSCRIPTS')

        num_transcripts_str = label_width_str.format('Transcripts')

        num_misassembled_by_blat_str = label_width_str.format('Misassembly candidates reported by GMAP (or BLAT)')
        num_misassembled_by_blast_str = label_width_str.format('Misassembly candidates reported by BLASTN')
        num_misassembled_together_str = label_width_str.format('Misassemblies')

        for i_transcripts in range(len(transcripts_metrics)):
            value_width_str = '{:<' + str(self.widths[i_transcripts + 1]) + '}'

            basic_metrics = transcripts_metrics[i_transcripts].basic_metrics

            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics

            name_str += value_width_str.format(transcripts_metrics[i_transcripts].label)

            if basic_metrics is not None and simple_metrics is not None:
                num_transcripts_str += value_width_str.format(basic_metrics.number)

            if simple_metrics is not None:
                num_misassembled_by_blat_str += value_width_str.format(simple_metrics.num_misassembled_by_blat)
                num_misassembled_together_str += value_width_str.format(simple_metrics.num_misassembled_together)
            else:
                num_misassembled_by_blat_str += value_width_str.format('*')
                num_misassembled_together_str += value_width_str.format('*')

            if simple_metrics is not None and args_blast:
                num_misassembled_by_blast_str += value_width_str.format(simple_metrics.num_misassembled_by_blast)
            else:
                num_misassembled_by_blast_str += value_width_str.format('*')

        out_handle = open(self.path_txt_misassemblies, 'w')

        out_handle.write(name_str + '\n')
        self.write_metric_str(num_transcripts_str + '\n\n', out_handle, len(transcripts_metrics), first_type='ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS (calculated with reference genome or with gene database)')

        self.write_metric_str(num_misassembled_by_blat_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(num_misassembled_by_blast_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(num_misassembled_together_str + '\n\n', out_handle, len(transcripts_metrics))

        out_handle.close()

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

        label_width_str = '{:<' + str(self.widths[0]) + '}'

        name_str = label_width_str.format('METRICS/TRANSCRIPTS')

        database_coverage_str = label_width_str.format('Database coverage')

        duplication_str = label_width_str.format('Duplication ratio')
        avg_num_transcripts_in_isoform_str = label_width_str.format('Avg. number of transcripts mapped to one isoform')

        assembled_well_genes_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled genes')
        assembled_fully_genes_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled genes')
        covered_well_genes_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered genes')
        covered_fully_genes_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered genes')

        assembled_well_isoforms_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled isoforms')
        assembled_fully_isoforms_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled isoforms')
        covered_well_isoforms_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered isoforms')
        covered_fully_isoforms_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered isoforms')

        assembled_exons_well_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-assembled exons')
        assembled_exons_fully_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled exons')
        mean_isoform_assembly_str = label_width_str.format('Mean isoform assembly')
        mean_isoform_cov_str = label_width_str.format('Mean isoform coverage')
        mean_exon_cov_str = label_width_str.format('Mean exon coverage')
        isoform_well_cov_exons_str = label_width_str.format('Avg. percentage of isoform ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_exon_threshold * 100)) + '%-covered exons')
        isoform_fully_cov_exons_str = label_width_str.format('Avg. percentage of isoform ' + str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_exon_threshold * 100)) + '%-covered exons')

        busco_complete_str = label_width_str.format('Complete')
        busco_partial_str = label_width_str.format('Partial')

        geneMarkS_T_genes_str = label_width_str.format('Genes')

        for i_transcripts in range(len(transcripts_metrics)):
            value_width_str = '{:<' + str(self.widths[i_transcripts + 1]) + '}'

            if transcripts_metrics[i_transcripts].assembly_completeness_metrics is not None:
                name_str += value_width_str.format(transcripts_metrics[i_transcripts].label)

                isoforms_coverage = transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage

                if isoforms_coverage is not None:
                    database_coverage_str += value_width_str.format(round(isoforms_coverage.fraction_annotation_mapped, PRECISION))

                    duplication_str += value_width_str.format(round(isoforms_coverage.avg_duplication_ratio, PRECISION))
                    avg_num_transcripts_in_isoform_str += value_width_str.format(round(isoforms_coverage.avg_num_transcripts_mapped_to_isoform, PRECISION))

                    assembled_well_genes_str += value_width_str.format(isoforms_coverage.num_well_assembled_genes)
                    assembled_fully_genes_str += value_width_str.format(isoforms_coverage.num_fully_assembled_genes)
                    covered_well_genes_str += value_width_str.format(isoforms_coverage.num_well_covered_genes)
                    covered_fully_genes_str += value_width_str.format(isoforms_coverage.num_fully_covered_genes)

                    assembled_well_isoforms_str += value_width_str.format(isoforms_coverage.num_well_assembled_isoforms)
                    assembled_fully_isoforms_str += value_width_str.format(isoforms_coverage.num_fully_assembled_isoforms)
                    covered_well_isoforms_str += value_width_str.format(isoforms_coverage.num_well_covered_isoforms)
                    covered_fully_isoforms_str += value_width_str.format(isoforms_coverage.num_fully_covered_isoforms)

                    assembled_exons_well_str += value_width_str.format(isoforms_coverage.num_well_assembled_exons)
                    assembled_exons_fully_str += value_width_str.format(isoforms_coverage.num_fully_assembled_exons)
                    mean_isoform_cov_str += value_width_str.format(round(isoforms_coverage.avg_covered_fraction, PRECISION))
                    mean_isoform_assembly_str += value_width_str.format(round(isoforms_coverage.avg_assembled_fraction, PRECISION))
                    mean_exon_cov_str += value_width_str.format(round(isoforms_coverage.avg_covered_fraction_exons, PRECISION))
                    isoform_well_cov_exons_str += value_width_str.format(round(isoforms_coverage.avg_percentage_isoform_well_covered_exons, PRECISION))
                    isoform_fully_cov_exons_str += value_width_str.format(round(isoforms_coverage.avg_percentage_isoform_fully_covered_exons, PRECISION))
                else:
                    database_coverage_str += value_width_str.format('*')

                    duplication_str += value_width_str.format('*')
                    avg_num_transcripts_in_isoform_str += value_width_str.format('*')

                    assembled_well_genes_str += value_width_str.format('*')
                    assembled_fully_genes_str += value_width_str.format('*')
                    covered_well_genes_str += value_width_str.format('*')
                    covered_fully_genes_str += value_width_str.format('*')

                    assembled_well_isoforms_str += value_width_str.format('*')
                    assembled_fully_isoforms_str += value_width_str.format('*')
                    covered_well_isoforms_str += value_width_str.format('*')
                    covered_fully_isoforms_str += value_width_str.format('*')

                    assembled_exons_well_str += value_width_str.format('*')
                    assembled_exons_fully_str += value_width_str.format('*')
                    mean_isoform_cov_str += value_width_str.format('*')
                    mean_isoform_assembly_str += value_width_str.format('*')
                    mean_exon_cov_str += value_width_str.format('*')
                    isoform_well_cov_exons_str += value_width_str.format('*')
                    isoform_fully_cov_exons_str += value_width_str.format('*')

                busco_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.busco_metrics

                if busco_metrics is not None:
                    busco_complete_str += value_width_str.format(round(busco_metrics.complete_completeness, PRECISION))
                    busco_partial_str += value_width_str.format(round(busco_metrics.partial_completeness, PRECISION))
                else:
                    busco_complete_str += value_width_str.format('*')
                    busco_partial_str += value_width_str.format('*')

                geneMarkS_T_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.geneMarkS_T_metrics

                if geneMarkS_T_metrics is not None:
                    geneMarkS_T_genes_str += value_width_str.format(geneMarkS_T_metrics.genes)
                else:
                    geneMarkS_T_genes_str += value_width_str.format('*')

        out_handle = open(self.path_txt_sensitivity, 'w')

        out_handle.write(name_str + '\n')

        self.write_metric_str(database_coverage_str + '\n', out_handle, len(transcripts_metrics), first_type='ASSEMBLY COMPLETENESS (SENSITIVITY)')

        self.write_metric_str(duplication_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(avg_num_transcripts_in_isoform_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(assembled_well_genes_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(assembled_fully_genes_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(covered_well_genes_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(covered_fully_genes_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(assembled_well_isoforms_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(assembled_fully_isoforms_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(covered_well_isoforms_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(covered_fully_isoforms_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(assembled_exons_well_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(assembled_exons_fully_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(mean_isoform_assembly_str + '\n\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(mean_isoform_cov_str + '\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(mean_exon_cov_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(isoform_well_cov_exons_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(isoform_fully_cov_exons_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(busco_complete_str + '\n', out_handle, len(transcripts_metrics), first_type='BUSCO METRICS')
        self.write_metric_str(busco_partial_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(geneMarkS_T_genes_str + '\n\n', out_handle, len(transcripts_metrics), first_type='GeneMarkS-T METRICS')

        out_handle.close()

        logger.info('      saved to {}'.format(self.path_txt_sensitivity))


        # out_handle.write('{:<100}'.format('Assembled isoforms') + str(isoforms_coverage.num_assembled) + '\n\n')

        # out_handle.write('{:<100}'.format('Assembled exons') + str(exons_coverage.num_assembled) + '\n\n')

        # out_handle.write('COVERAGE BY EACH MAPPED TRANSCRIPT SEPARATELY:\n')
        # out_handle.write('DISTRIBUTIONS FOR ALIGNED TRANSCRIPTS:\n')

        # out_handle.write('{:<100}'.format('Avg. fraction of exon assembled in isoform') + str(round(isoforms_coverage.avg_assembled_fraction_exon, precision)) + '\n\n')

        # out_handle.write('{:<100}'.format('Avg. percentage of isoform well-assembled exons') + str(round(isoforms_coverage.avg_percentage_well_assembled_exons, precision)) + '\n')
        # out_handle.write('{:<100}'.format('Avg. percentage of isoform fully-assembled exons') + str(round(isoforms_coverage.avg_percentage_fully_assembled_exons, precision)) + '\n\n')

        # out_handle.write('DISTRIBUTIONS FOR MAPPED ISOFORMS:\n')
        # out_handle.write('{:<100}'.format('Avg. fraction of isoform assembled') + str(round(isoforms_coverage.avg_over_isoforms_assembled_fraction, precision)) + '\n\n')

        # out_handle.write('{:<100}'.format('Avg. fraction of exon assembled in isoform') + str(round(isoforms_coverage.avg_over_isoforms_avg_assembled_fraction_exon, precision)) + '\n\n')

        # out_handle.write('{:<100}'.format('Avg. percentage of isoform well-assembled exons') + str(round(isoforms_coverage.avg_over_isoforms_percentage_well_assembled_exons, precision)) + '\n')
        # out_handle.write('{:<100}'.format('Avg. percentage of isoform fully-assembled exons') + str(round(isoforms_coverage.avg_over_isoforms_percentage_fully_assembled_exons, precision)) + '\n\n\n')


        # out_handle.write('DISTRIBUTIONS FOR MAPPED EXONS:\n')
        # out_handle.write('{:<100}'.format('Avg. fraction of exon assembled') + str(round(exons_coverage.avg_over_exons_assembled_fraction, precision)) + '\n\n\n')

        # out_handle.write('FOR ALL ALIGNED TRANSCRIPTS:\n')
        # out_handle.write('{:<100}'.format('Well-assembled in average isoforms') + str(isoforms_coverage.num_well_assembled_in_average) + '\n')
        # out_handle.write('{:<100}'.format('Fully-assembled in average isoforms') + str(isoforms_coverage.num_fully_assembled_in_average) + '\n\n')

        # out_handle.write('{:<100}'.format('Well-assembled in average exons') + str(exons_coverage.num_well_assembled_in_average) + '\n')
        # out_handle.write('{:<100}'.format('Fully-assembled in average exons') + str(exons_coverage.num_fully_assembled_in_average) + '\n\n')

        # out_handle.write('{:<100}'.format('Avg. fraction of exon assembled') + str(round(exons_coverage.avg_assembled_fraction, precision)) + '\n\n\n')

        # out_handle.write('COVERAGE BY ALL MAPPED TRANSCRIPTS AT ONCE:\n')
        # out_handle.write('DISTRIBUTIONS FOR ALIGNED TRANSCRIPTS:\n')

        # out_handle.write('{:<100}'.format('Avg. fraction of exon covered in isoform') + str(round(isoforms_coverage.avg_covered_fraction_exon, precision)) + '\n\n')

        # out_handle.write('FOR ALL ALIGNED TRANSCRIPTS:\n')
        # out_handle.write('{:<100}'.format('Well-covered exons') + str(exons_coverage.num_well_covered) + '\n')
        # out_handle.write('{:<100}'.format('Fully-covered exons') + str(exons_coverage.num_fully_covered) + '\n\n')


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

        label_width_str = '{:<' + str(self.widths[0]) + '}'

        name_str = label_width_str.format('METRICS/TRANSCRIPTS')

        unannotated_str = label_width_str.format('Unannotated')

        matched_well_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_transcript_threshold * 100)) + '%-matched')
        matched_fully_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_transcript_threshold * 100)) + '%-matched')

        mean_transcript_match_str = label_width_str.format('Mean fraction of transcript matched')

        mean_block_match_str = label_width_str.format('Mean fraction of block matched')

        matched_blocks_well_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_block_threshold * 100)) + '%-matched blocks')
        matched_blocks_fully_str = label_width_str.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_block_threshold * 100)) + '%-matched blocks')

        matched_len_str = label_width_str.format('Matched length')
        unmatched_len_str = label_width_str.format('Unmatched length')

        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].assembly_correctness_metrics is not None:
                value_width_str = '{:<' + str(self.widths[i_transcripts + 1]) + '}'

                transcripts_coverage = transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage

                if transcripts_coverage is not None:
                    name_str += value_width_str.format(transcripts_metrics[i_transcripts].label)
                    unannotated_str += value_width_str.format(transcripts_coverage.num_unannotated_transcripts)

                    matched_well_str += value_width_str.format(transcripts_coverage.num_well_covered_transcripts)
                    matched_fully_str += value_width_str.format(transcripts_coverage.num_fully_covered_transcripts)

                    mean_transcript_match_str += value_width_str.format(round(transcripts_coverage.avg_covered_fraction_whole_transcript, precision))

                    mean_block_match_str += value_width_str.format(round(transcripts_coverage.avg_covered_fraction_block, precision))

                    matched_blocks_well_str += value_width_str.format(round(transcripts_coverage.percentage_well_covered_blocks, precision))
                    matched_blocks_fully_str += value_width_str.format(round(transcripts_coverage.percentage_fully_covered_blocks, precision))

                    matched_len_str += value_width_str.format(transcripts_coverage.matched_len)
                    unmatched_len_str += value_width_str.format(transcripts_coverage.unmatched_len)
                else:
                    unannotated_str += value_width_str.format('*')

                    matched_well_str += value_width_str.format('*')
                    matched_fully_str += value_width_str.format('*')

                    mean_transcript_match_str += value_width_str.format('*')

                    mean_block_match_str += value_width_str.format('*')

                    matched_blocks_well_str += value_width_str.format('*')
                    matched_blocks_fully_str += value_width_str.format('*')

                    matched_len_str += value_width_str.format('*')
                    unmatched_len_str += value_width_str.format('*')

        out_handle = open(self.path_txt_specificity, 'w')

        out_handle.write(name_str + '\n')

        self.write_metric_str(unannotated_str + '\n\n', out_handle, len(transcripts_metrics), first_type='ASSEMBLY SPECIFICITY')

        self.write_metric_str(matched_well_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(matched_fully_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(mean_transcript_match_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(mean_block_match_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(matched_blocks_well_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(matched_blocks_fully_str + '\n\n', out_handle, len(transcripts_metrics))

        self.write_metric_str(matched_len_str + '\n', out_handle, len(transcripts_metrics))
        self.write_metric_str(unmatched_len_str + '\n\n', out_handle, len(transcripts_metrics))

        out_handle.close()

        logger.info('      saved to {}'.format(self.path_txt_specificity))

        # FOR ALL TRANSCRIPTS ALIGNMENTS:
        # out_handle.write('FOR ALL TRANSCRIPTS ALIGNMENTS:\n')

        # out_handle.write('{:<100}'.format('Annotated transcripts') + str(transcripts_coverage.num_annotated_transcripts) + '\n')

        # out_handle.write('DISTRIBUTIONS FOR TRANSCRIPTS ALIGNMENTS:\n')
        # out_handle.write('{:<100}'.format('Avg. fraction of aligned part annotated') + str(round(transcripts_coverage.avg_covered_fraction_aligned_part, precision)) + '\n\n')

        # out_handle.write('{:<100}'.format('Avg. fraction of block annotated in transcript:') + str(round(transcripts_coverage.avg_covered_fraction_block_in_t, precision)) + '\n\n')

        # out_handle.write('{:<100}'.format('Avg. percentage of well-annotated blocks in transcript') + str(round(transcripts_coverage.avg_percentage_well_covered_blocks_in_t, precision)) + '\n')
        # out_handle.write('{:<100}'.format('Avg. percentage of fully-annotated blocks in transcript') + str(round(transcripts_coverage.avg_percentage_fully_covered_blocks_in_t, precision)) + '\n\n\n')

        # out_handle.write('FOR ALL TRANSCRIPTS MAPPED TO SPECIFIC ISOFORM:\n')
        # out_handle.write('{:<100}'.format('Avg. fraction of whole transcript annotated') + str(round(transcripts_coverage.avg_over_isoforms_covered_fraction_whole_transcript, precision)) + '\n')
        # out_handle.write('{:<100}'.format('Avg. fraction of aligned part annotated') + str(round(transcripts_coverage.avg_over_isoforms_covered_fraction_aligned_part, precision)) + '\n\n')

        # out_handle.write('{:<100}'.format('Avg. fraction of total transcript annotated length over total whole transcript lengths') + str(round(transcripts_coverage.avg_covered_fraction_whole_transcript_over_tot_len, precision)) + '\n')
        # out_handle.write('{:<100}'.format('Fraction of total transcript annotated length over total aligned parts') + str(round(transcripts_coverage.avg_covered_fraction_aligned_part_over_tot_len, precision)) + '\n\n\n')