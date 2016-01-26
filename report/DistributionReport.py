__author__ = 'letovesnoi'

from report import UtilsPictures


class DistributionReport():
    """Class which generate distributions"""

    def __init__(self, transcripts_metrics, db_genes_metrics, output_dir, logger, precision):
        self.output_dir = output_dir

        logger.print_timestamp('  ')
        logger.info('  Getting DISTRIBUTION report...')
        # for creating PDF file with all plots and tables:
        self.pdf_plots_figures = []

        # BASIC PLOTS:
        self.get_basic_plots(transcripts_metrics, db_genes_metrics, precision)

        # SPECIFICITY:
        if len(transcripts_metrics) != 0:
            if transcripts_metrics[0].assembly_correctness_metrics != None:
                self.get_specificity_plots(transcripts_metrics, self.output_dir, precision)

        # SENSITIVITY:
        if len(transcripts_metrics) != 0:
            if transcripts_metrics[0].assembly_completeness_metrics != None:
                self.get_sensitivity_plots(transcripts_metrics, self.output_dir, precision)

        logger.info('  Done.')


    def get_basic_plots(self, transcripts_metrics, db_genes_metrics, precision):
        # BASIC AND SIMPLE 1
        # Transcripts len distribution from aligned and assembled transcripts(avg, min, max, total)
        # Transcripts alignment len distribution (avg, min, max, total)
        # Isoforms len without introns distribution (avg, min, max, total)
        # Number of transcripts with length inside isoform  without introns length range
        self.get_transcript_length_plot(db_genes_metrics, transcripts_metrics, self.output_dir, precision, short_report_visible=True)

        # BASIC AND SIMPLE 4
        # Blocks len distribution (tot transcripts alignment length, avg, min, max)
        # Exons len distribution (tot isoforms len without introns, avg, min, max)
        self.get_block_length_plot(db_genes_metrics, transcripts_metrics, self.output_dir, precision, short_report_visible=False)

        # BASIC AND SIMPLE 0
        # Transcripts alignment fraction distribution (avg, min, max, total)
        self.get_x_aligned_plot(transcripts_metrics, self.output_dir, precision, short_report_visible=False)

        # BASIC AND SIMPLE 3
        # Number of blocks per transcript distribution  (tot, avg, min, max)
        # Number of exons per isoform distribution  (tot, avg, min, max)
        self.get_blocks_per_alignment_plot(db_genes_metrics, transcripts_metrics, self.output_dir, precision, short_report_visible=False)

        # BASIC AND SIMPLE 5
        # Multiple-aligned transcripts distribution
        # Number of unaligned / aligned / alignments / unique aligned / multiple-aligned transcripts
        self.get_alignment_multiplicity_plot(transcripts_metrics, self.output_dir, short_report_visible=False)

        # BASIC AND SIMPLE 6
        # Mismatch number per transcript distribution (tot, avg)
        self.get_mismatch_rate_plot(transcripts_metrics, self.output_dir, precision, short_report_visible=True)

        # BASIC AND SIMPLE 9
        # Query gap number per transcript distribution (tot, avg, min, max)
        # self.get_gaps_per_alignment_plot(transcripts_metrics, self.output_dir, precision, short_report_visible=False)

        # BASIC AND SIMPLE 10
        # Query gap len distribution (tot, avg, min, max)
        # self.get_gap_length_plot(transcripts_metrics, self.output_dir, precision, short_report_visible=True)

        # BASIC AND SIMPLE 11
        # Nx:
        if len(transcripts_metrics) != 0:
            if transcripts_metrics[0].basic_metrics != None:
                list_of_labels = []
                lists_of_lengths = []
                for i_transcripts in range(len(transcripts_metrics)):
                    list_of_labels.append(transcripts_metrics[i_transcripts].label)
                    lists_of_lengths.append(transcripts_metrics[i_transcripts].basic_metrics.lengths)
                UtilsPictures.Nx_plot(list_of_labels, lists_of_lengths, self.output_dir, self.pdf_plots_figures, title_str='Nx', short_report_visible=False, reference_lengths=None)

        # BASIC AND SIMPLE 12
        # NAx:
        if len(transcripts_metrics) != 0:
            if transcripts_metrics[0].simple_metrics != None:
                list_of_labels = []
                lists_of_lengths = []
                for i_transcripts in range(len(transcripts_metrics)):
                    list_of_labels.append(transcripts_metrics[i_transcripts].label)
                    lists_of_lengths.append(transcripts_metrics[i_transcripts].simple_metrics.alignment_lengths)
                UtilsPictures.Nx_plot(list_of_labels, lists_of_lengths, self.output_dir, self.pdf_plots_figures, title_str='NAx', short_report_visible=True, reference_lengths=None)


    # BASIC AND SIMPLE 1
    # Transcripts len distribution (avg, min, max, total)
    # Transcripts alignment len distribution (avg, min, max, total)
    # Isoforms len without introns distribution (avg, min, max, total)
    # Number of transcripts with length inside isoform  without introns length range
    def get_transcript_length_plot(self, db_genes_metrics, transcripts_metrics, out_dir, precision, short_report_visible):
        # set fields of plots:
        title_str = 'assembled transcripts / database isoforms length'

        if db_genes_metrics != None:
            isoforms_distribution = db_genes_metrics.isoforms_len_distribution
            str = 'Gene database'
            list_label = [('avg', round(db_genes_metrics.avg_isoform_len, precision))]
                          # ('min', basic_isoforms_metrics.min_len_wout_introns),
                          # ('max', basic_isoforms_metrics.max_len_wout_introns),
                          # ('tot', basic_isoforms_metrics.tot_len_wout_introns)]
            isoforms_label = UtilsPictures.get_label(str, list_label)
        else:
            isoforms_distribution = None
            isoforms_label = None

        transcripts_distributions = []
        transcripts_labels = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].basic_metrics != None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].basic_metrics.len_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].basic_metrics.avg_len, precision))]
                              # ('min', transcripts_metrics[i_transcripts].basic_metrics.min_len),
                              # ('max', transcripts_metrics[i_transcripts].basic_metrics.max_len),
                              # ('tot', transcripts_metrics[i_transcripts].basic_metrics.tot_len)]
                # if basic_isoforms_metrics != None:
                #     list_label.append(('outliers', round(transcripts_metrics[i_transcripts].basic_metrics.percent_outliers_len * 100, precision)))

                transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        label_x = 'length'
        label_y = 'number of transcripts / isoforms'
        name_fig = 'transcript_length'
        y_log_scale = True

        # plot:
        UtilsPictures.plot_compare_distribution(out_dir, title_str, isoforms_distribution, isoforms_label,
                                                transcripts_distributions, transcripts_labels, label_x, label_y, name_fig,
                                                short_report_visible, self.pdf_plots_figures, y_log_scale=y_log_scale)


    # BASIC AND SIMPLE 4
    # Blocks len distribution (tot transcripts alignment length, avg, min, max)
    # Exons len distribution (tot isoforms len without introns, avg, min, max)
    def get_block_length_plot(self, db_genes_metrics, transcripts_metrics, out_dir, precision, short_report_visible):
        # set fields of plots:
        title_str = 'blocks / exons length'

        if db_genes_metrics != None:
            isoforms_distribution = db_genes_metrics.exons_len_distribution
            str = 'Gene database'
            list_label = [('avg', round(db_genes_metrics.avg_exon_len, precision))]
                          # ('min', basic_isoforms_metrics.min_exon_len),
                          # ('max', basic_isoforms_metrics.max_exon_len),
                          # ('tot', basic_isoforms_metrics.tot_len_wout_introns)]
            isoforms_label = UtilsPictures.get_label(str, list_label)
        else:
            isoforms_distribution = None
            isoforms_label = None

        transcripts_distributions = []
        transcripts_labels = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].simple_metrics.blocks_len_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_block_len, precision))]
                              # ('min', transcripts_metrics[i_transcripts].simple_metrics.min_block_len),
                              # ('max', transcripts_metrics[i_transcripts].simple_metrics.max_block_len),
                              # ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_alignment_len)]
                # if basic_isoforms_metrics != None:
                #     list_label.append(('outliers', round(transcripts_metrics[i_transcripts].simple_metrics.percent_outliers_blocks_len * 100, precision)))

                transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        label_x = 'length'
        label_y = 'number of blocks / exons'
        name_fig = 'block_length'
        y_log_scale = True

        # plot:
        UtilsPictures.plot_compare_distribution(out_dir, title_str, isoforms_distribution, isoforms_label,
                                                transcripts_distributions, transcripts_labels, label_x, label_y, name_fig,
                                                short_report_visible, self.pdf_plots_figures, y_log_scale=y_log_scale)


    # # BASIC AND SIMPLE 0
    # Transcripts alignment fraction distribution (avg, min, max, total)
    def get_x_aligned_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # set fields of plots:
        title_str = 'transcript aligned fraction'

        label_x = 'aligned fraction'
        label_y = 'number of transcripts'
        name_fig = 'x-aligned'
        step = 0.1

        transcripts_distribution = []
        transcripts_labels = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distribution.append(transcripts_metrics[i_transcripts].simple_metrics.fraction_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_fraction, precision))]
                              # ('min', round(transcripts_metrics[i_transcripts].simple_metrics.min_fraction, precision)),
                              # ('max', round(transcripts_metrics[i_transcripts].simple_metrics.max_fraction, precision))]

                transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        # plot:
        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics,
                                             transcripts_labels, label_x, label_y, name_fig, short_report_visible,
                                             self.pdf_plots_figures, def_step=step)


    # BASIC AND SIMPLE 3
    # Number of blocks per transcript distribution  (tot, avg, min, max)
    # Number of exons per isoform distribution  (tot, avg, min, max)
    def get_blocks_per_alignment_plot(self, db_genes_metrics, transcripts_metrics, out_dir, precision, short_report_visible):
        # set fields of plots:
        title_str = 'number of blocks / exons per alignment / isoform'

        if db_genes_metrics != None:
            isoforms_distribution = db_genes_metrics.exons_num_distribution
            str = 'Gene database'
            list_label = [('avg', round(db_genes_metrics.avg_exons_num, precision)),
                          # ('min', basic_isoforms_metrics.min_exons_num),
                          # ('max', basic_isoforms_metrics.max_exons_num),
                          ('tot', db_genes_metrics.tot_exons_num)]
            isoforms_label = UtilsPictures.get_label(str, list_label)
        else:
            isoforms_distribution = None
            isoforms_label = None

        transcripts_distributions = []
        transcripts_labels = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].simple_metrics.blocks_num_distribution)
                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_blocks_num, precision)),
                              # ('min', transcripts_metrics[i_transcripts].simple_metrics.min_blocks_num),
                              # ('max', transcripts_metrics[i_transcripts].simple_metrics.max_blocks_num),
                              ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_blocks_num)]
                # if basic_isoforms_metrics != None:
                #     list_label.append(('outliers', round(transcripts_metrics[i_transcripts].simple_metrics.percent_outliers_blocks_num * 100, precision)))

                transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        label_x = 'number of blocks / exons'
        label_y = 'number of alignments / isoforms'
        name_fig = 'blocks_per_alignment'
        y_log_scale = True
        step = 1

        # plot:
        UtilsPictures.plot_compare_distribution(out_dir, title_str, isoforms_distribution, isoforms_label,
                                                transcripts_distributions, transcripts_labels, label_x, label_y, name_fig,
                                                short_report_visible, self.pdf_plots_figures, y_log_scale=y_log_scale,
                                                def_step=step)


    # BASIC AND SIMPLE 5
    # Multiple-aligned transcripts distribution
    # Number of unaligned / aligned / alignments / unique aligned / multiple-aligned transcripts
    def get_alignment_multiplicity_plot(self, transcripts_metrics, out_dir, short_report_visible):
        title_str = 'alignment multiplicity'
        label_x = 'number of alignments'
        label_y = 'number of transcripts'
        name_fig = 'alignment_multiplicity'
        y_log_scale = True
        step = 1

        transcripts_labels = []
        transcripts_distributions = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].simple_metrics.mult_alignment_distribution)
                str = transcripts_metrics[i_transcripts].label
                list_label = [('unaligned', transcripts_metrics[i_transcripts].simple_metrics.num_unaligned),
                              ('aligned', transcripts_metrics[i_transcripts].simple_metrics.num_aligned),
                              ('alignments', transcripts_metrics[i_transcripts].simple_metrics.num_alignments),
                              ('unique aligned', transcripts_metrics[i_transcripts].simple_metrics.num_unique_aligned),
                              ('multiple aligned', transcripts_metrics[i_transcripts].simple_metrics.num_mul_aligned)]

                transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_distribution(out_dir, title_str, None, None, transcripts_distributions,
                                                transcripts_labels, label_x, label_y, name_fig, short_report_visible,
                                                self.pdf_plots_figures, y_log_scale=y_log_scale, def_step=step)

    # BASIC AND SIMPLE 6
    # Mismatch number per transcript distribution (tot, avg)
    def get_mismatch_rate_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        title_str = 'substitution errors per alignment'
        label_x = 'number of mismatches'
        label_y = 'number of alignments'
        name_fig = 'mismatch_rate'
        y_log_scale = True
        step = 1

        transcripts_labels = []
        transcripts_distributions = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].simple_metrics.mismatch_num_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_mismatch_num, precision)),
                              ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_mismatch_num)]

                transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_distribution(out_dir, title_str, None, None, transcripts_distributions,
                                                transcripts_labels, label_x, label_y, name_fig, short_report_visible,
                                                self.pdf_plots_figures, y_log_scale=y_log_scale, def_step=step)

    # BASIC AND SIMPLE 7
    # Target gap number per transcript distribution (tot, avg, min, max)
    # Number of introns per isoform distribution  (tot, avg, min, max)
    def get_basic_simple_7(self, basic_isoforms_metrics, transcripts_metrics, out_dir, precision, short_report_visible):
        # set fields of plots:
        title_str = 'Basic and Simple 7'

        isoforms_distribution = None
        isoforms_label = None
        if basic_isoforms_metrics != None:
            isoforms_distribution = basic_isoforms_metrics.introns_num_distribution

            str = 'Number of introns per isoform'
            list_label = [('min', basic_isoforms_metrics.min_introns_num),
                          ('max', basic_isoforms_metrics.max_introns_num),
                          ('avg', round(basic_isoforms_metrics.avg_introns_num, precision)),
                          ('tot', basic_isoforms_metrics.tot_introns_num)]
            isoforms_label = UtilsPictures.get_label(str, list_label)

        transcripts_distributions = []
        transcripts_labels = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                t_distribution_i_0 = transcripts_metrics[i_transcripts].simple_metrics.tgap_num_distribution
                transcripts_distributions.append([t_distribution_i_0])

                str = 'Target gap number per transcript [{}]'.format(transcripts_metrics[i_transcripts].label)
                list_label = [('min', transcripts_metrics[i_transcripts].simple_metrics.min_tgap_num),
                              ('max', transcripts_metrics[i_transcripts].simple_metrics.max_tgap_num),
                              ('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_tgap_num, precision)),
                              ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_tgap_num)]
                if basic_isoforms_metrics != None:
                    list_label.append(('outliers', round(transcripts_metrics[i_transcripts].simple_metrics.percent_outliers_tgap_num * 100, precision)))
                t_label_i_0 = UtilsPictures.get_label(str, list_label)

                transcripts_labels.append([t_label_i_0])

        label_x = 'number of gaps/introns'
        label_y = 'number of transcripts/isoforms'
        name_fig = '7'
        y_log_scale = True
        step = 1

        # plot:
        UtilsPictures.plot_compare_distribution(out_dir, title_str, isoforms_distribution, isoforms_label, transcripts_distributions, transcripts_labels,
                                                label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, y_log_scale=y_log_scale, def_step=step)


    # BASIC AND SIMPLE 9
    # Query gap number per transcript distribution (tot, avg, min, max)
    def get_gaps_per_alignment_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        title_str = 'number of gaps per alignment'
        label_x = 'number of gaps'
        label_y = 'number of alignments'
        name_fig = 'gaps_per_alignment'
        y_log_scale = True
        step = 1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distribution.append(transcripts_metrics[i_transcripts].simple_metrics.qgap_num_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_qgap_num, precision)),
                              ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_qgap_num)]

                transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, y_log_scale=y_log_scale, def_step=step)

    # BASIC AND SIMPLE 10
    # Query gap len distribution (tot, avg, min, max)
    def get_gap_length_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        title_str = 'alignment gaps length'
        label_x = 'gap length'
        label_y = 'number of alignments'
        name_fig = 'gap_length'
        y_log_scale =True
        step = None

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distribution.append(transcripts_metrics[i_transcripts].simple_metrics.qgap_len_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_qgap_len, precision)),
                              ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_qgap_len)]

                transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, y_log_scale, step)


    def get_x_matched_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        title_str = 'fraction of transcript matched'
        label_x = 'fraction of transcript matched'
        label_y = 'number of transcripts'
        name_fig = 'x-matched'
        y_log_scale = False
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage.covered_fraction_whole_transcript_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage.avg_covered_fraction_whole_transcript, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, y_log_scale, def_step=step)


    def get_x_matched_blocks_plot(self, transcripts_metrics, out_dir, precision, short_report_visible=False):
        title_str = 'fraction of block matched'
        label_x = 'fraction of block matched'
        label_y = 'number of blocks'
        name_fig = 'x-matched_blocks'
        y_log_scale = False
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage.block_fraction_matched_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage.avg_covered_fraction_block, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, y_log_scale, def_step=step)


    def get_specificity_plots(self, transcripts_metrics, out_dir, precision):
        # DISTRIBUTIONS FOR TRANSCRIPTS ALIGNMENTS
        self.get_x_matched_plot(transcripts_metrics, out_dir, precision, short_report_visible=True)

        self.get_x_matched_blocks_plot(transcripts_metrics, out_dir, precision, short_report_visible=False)


    # CONSIDER COVERED BASES BY EACH MAPPED TRANSCRIPT SEPARATELY:
    def get_x_assembled_exons_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        title_str = 'fraction of exon captured by a single transcript'
        label_x = 'fraction of exon assembled'
        label_y = 'number of exons'
        name_fig = 'x-assembled_exons'
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.assembled_fraction_exons_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_assembled_fraction_exons, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)


    # distributions for isoforms:
    def get_x_assembled_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        title_str = 'fraction of isoform captured by a single transcript'
        label_x = 'fraction of isoform assembled'
        label_y = 'number of isoforms'
        name_fig = 'x-assembled'
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.assembled_fraction_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_assembled_fraction, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)


    # CONSIDER COVERED BASES (BY ALL MAPPED TRANSCRIPTS) WITHOUT OVERLAPS:
    # distributions for exons::
    def get_x_covered_exons_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        title_str = 'fraction of exon covered by all alignments'
        label_x = 'fraction of exon covered'
        label_y = 'number of exons'
        name_fig = 'x-covered_exons'
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.covered_fraction_exons_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_covered_fraction_exons, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)

    # distributions for isoforms:
    def get_x_covered_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        title_str = 'fraction of isoform covered by all alignments'
        label_x = 'fraction of isoform covered'
        label_y = 'number of isoforms'
        name_fig = 'x-covered'
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.covered_fraction_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_covered_fraction, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)

    def get_avg_exon_covered_fraction_over_isoforms_distribution(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        title_str = 'average fraction of exon in isoform covered'
        label_x = 'avg fraction of exon covered'
        label_y = 'number of isoforms'
        name_fig = 'avg_exon_in_isoform_fraction_covered'
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_covered_fraction_exon_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_covered_fraction_exon, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)

    def get_percentage_well_covered_exons_over_isoforms_distribution(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        title_str = 'percentage of well-covered exons in isoform '
        label_x = 'percentage of isoform well-covered exons'
        label_y = 'number of isoforms'
        name_fig = 'well_covered_exons_in_isoform'
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.percentage_well_covered_exons_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_percentage_well_covered_exons, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)

    def get_percentage_fully_covered_exons_over_isoforms_distribution(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        title_str = 'percentage of fully-covered exons in isoform'
        label_x = 'percentage of isoform fully-covered exons'
        label_y = 'number of isoforms'
        name_fig = 'fully_covered_exons_in_isoform'
        step = 0.1

        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.percentage_fully_covered_exons_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_percentage_fully_covered_exons, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_histogram(out_dir, title_str, transcripts_distribution, transcripts_metrics, transcripts_labels,
                                             label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)


    def get_alignments_per_isoform_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        title_str = 'number of alignments per isoform'
        label_x = 'number of alignments'
        label_y = 'number of isoforms'
        name_fig = 'alignments_per_isoform'
        y_log_scale = True
        x_log_scale = False
        step = 1

        transcripts_labels = []
        transcripts_distributions = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distributions.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.num_transcripts_mapped_to_isoform_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_num_transcripts_mapped_to_isoform, precision))]

            transcripts_labels.append(UtilsPictures.get_label(str, list_label))

        UtilsPictures.plot_compare_distribution(out_dir, title_str, None, None, transcripts_distributions, transcripts_labels, label_x, label_y, name_fig,
                                                short_report_visible, self.pdf_plots_figures, x_log_scale=x_log_scale, y_log_scale=y_log_scale, def_step=step)



    def get_sensitivity_plots(self, transcripts_metrics, out_dir, precision):
        self.get_x_assembled_plot(transcripts_metrics, out_dir, precision, short_report_visible=True)

        self.get_x_covered_plot(transcripts_metrics, out_dir, precision, short_report_visible=True)

        self.get_x_assembled_exons_plot(transcripts_metrics, out_dir, precision, short_report_visible=False)

        self.get_x_covered_exons_plot(transcripts_metrics, out_dir, precision, short_report_visible=False)

        self.get_alignments_per_isoform_plot(transcripts_metrics, out_dir, precision, short_report_visible=False)

        # self.get_percentage_well_covered_exons_over_isoforms_distribution(transcripts_metrics, out_dir, precision, short_report_visible=False)

        # self.get_percentage_fully_covered_exons_over_isoforms_distribution(transcripts_metrics, out_dir, precision, short_report_visible=False)