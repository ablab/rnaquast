__author__ = 'letovesnoi'

import os

from report import UtilsPictures


class DistributionReport():
    """Class which generate distributions"""

    def __init__(self, transcripts_metrics, db_genes_metrics, output_dir, logger, precision):
        self.output_dir = output_dir

        logger.print_timestamp('  ')
        logger.info('  Getting DISTRIBUTION report...')
        # for creating PDF file with all plots and tables:
        self.short_report_plots = []

        # BASIC PLOTS:
        self.get_basic_plots(transcripts_metrics, db_genes_metrics, precision)

        # SPECIFICITY:
        if len(transcripts_metrics) != 0:
            if transcripts_metrics[0].assembly_correctness_metrics is not None:
                self.get_specificity_plots(transcripts_metrics, self.output_dir, precision)

        # SENSITIVITY:
        if len(transcripts_metrics) != 0:
            if transcripts_metrics[0].assembly_completeness_metrics is not None:
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
            if transcripts_metrics[0].basic_metrics is not None:
                list_of_labels = []
                lists_of_lengths = []
                for i_transcripts in range(len(transcripts_metrics)):
                    list_of_labels.append(transcripts_metrics[i_transcripts].label)
                    lists_of_lengths.append(transcripts_metrics[i_transcripts].basic_metrics.lengths)

                title_name = 'Nx'
                label_x = 'x'
                label_y = 'Contig length'
                name_fig = 'Nx'
                short_report_visible = False

                nx_plot = \
                    UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, self.output_dir)

                nx_plot.Nx_plot(list_of_labels, lists_of_lengths, self.short_report_plots)

        # BASIC AND SIMPLE 12
        # NAx:
        if len(transcripts_metrics) != 0:
            if transcripts_metrics[0].simple_metrics is not None:
                list_of_labels = []
                lists_of_lengths = []
                for i_transcripts in range(len(transcripts_metrics)):
                    list_of_labels.append(transcripts_metrics[i_transcripts].label)
                    lists_of_lengths.append(transcripts_metrics[i_transcripts].simple_metrics.alignment_lengths)

                caption = 'Nx plot for transcripts. Nx is a maximal number $N$, such that the total length of all ' \
                          'transcripts longer than $N$ bp is at least $x%$ of the total length of all transcripts.'
                title_name = 'NAx'
                label_x = 'x'
                label_y = 'Contig length'
                name_fig = 'NAx'
                short_report_visible = True

                nax_plot = \
                    UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, self.output_dir,
                                       caption=caption)

                nax_plot.Nx_plot(list_of_labels, lists_of_lengths, self.short_report_plots)


    # BASIC AND SIMPLE 1
    # Transcripts len distribution (avg, min, max, total)
    # Transcripts alignment len distribution (avg, min, max, total)
    # Isoforms len without introns distribution (avg, min, max, total)
    # Number of transcripts with length inside isoform  without introns length range
    def get_transcript_length_plot(self, db_genes_metrics, transcripts_metrics, out_dir, precision, short_report_visible):
        if db_genes_metrics is not None:
            isoforms_distribution = db_genes_metrics.isoforms_len_distribution
            str = 'Gene database'
            list_label = [('avg', round(db_genes_metrics.avg_isoform_len, precision))]

            isoforms_label = DistributionReport.get_label(str, list_label)
        else:
            isoforms_distribution = None
            isoforms_label = None

        transcripts_distributions = []
        transcripts_labels = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].basic_metrics is not None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].basic_metrics.len_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].basic_metrics.avg_len, precision))]

                transcripts_labels.append(DistributionReport.get_label(str, list_label))

        caption = \
            'Plot showing cumulative transcript length distribution. Each point represents the ' \
            'number of transcripts in the assembly with the corresponding length or longer; black dashed line ' \
            'corresponds to the database isoforms; the plot is given in logarithmic scale.'
        title_name = 'transcript / isoform length'
        label_x = 'length'
        label_y = 'number of transcripts / isoforms'
        name_fig = 'transcript_length'

        transcript_len_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, isoforms_label, isoforms_distribution,
                               y_log_scale=True, caption=caption)

        transcript_len_plot.plot_compare_distribution(self.short_report_plots)


    # BASIC AND SIMPLE 4
    # Blocks len distribution (tot transcripts alignment length, avg, min, max)
    # Exons len distribution (tot isoforms len without introns, avg, min, max)
    def get_block_length_plot(self, db_genes_metrics, transcripts_metrics, out_dir, precision, short_report_visible):
        if db_genes_metrics is not None:
            isoforms_distribution = db_genes_metrics.exons_len_distribution
            str = 'Gene database'
            list_label = [('avg', round(db_genes_metrics.avg_exon_len, precision))]

            isoforms_label = DistributionReport.get_label(str, list_label)
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

                transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'blocks / exons length'
        label_x = 'length'
        label_y = 'number of blocks / exons'
        name_fig = 'block_length'

        block_len_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, isoforms_label, isoforms_distribution,
                               y_log_scale=True)

        block_len_plot.plot_compare_distribution(self.short_report_plots)


    # # BASIC AND SIMPLE 0
    # Transcripts alignment fraction distribution (avg, min, max, total)
    def get_x_aligned_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        transcripts_distributions = []
        transcripts_labels = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].simple_metrics.fraction_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_fraction, precision))]

                transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'transcript aligned fraction'
        label_x = 'aligned fraction'
        label_y = 'number of transcripts'
        name_fig = 'x-aligned'

        x_aligned_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, y_log_scale=True, def_step=0.1)

        x_aligned_plot.plot_compare_distribution(self.short_report_plots)


    # BASIC AND SIMPLE 3
    # Number of blocks per transcript distribution  (tot, avg, min, max)
    # Number of exons per isoform distribution  (tot, avg, min, max)
    def get_blocks_per_alignment_plot(self, db_genes_metrics, transcripts_metrics, out_dir, precision, short_report_visible):
        if db_genes_metrics is not None:
            isoforms_distribution = db_genes_metrics.exons_num_distribution
            str = 'Gene database'
            list_label = [('avg', round(db_genes_metrics.avg_exons_num, precision)),
                          # ('min', basic_isoforms_metrics.min_exons_num),
                          # ('max', basic_isoforms_metrics.max_exons_num),
                          ('tot', db_genes_metrics.tot_exons_num)]
            isoforms_label = DistributionReport.get_label(str, list_label)
        else:
            isoforms_distribution = None
            isoforms_label = None

        transcripts_distributions = []
        transcripts_labels = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics != None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].simple_metrics.blocks_num_distribution)
                str = transcripts_metrics[i_transcripts].label
                list_label = \
                    [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_blocks_num, precision)),
                     ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_blocks_num)]

                transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'number of blocks / exons per alignment / isoform'
        label_x = 'number of blocks / exons'
        label_y = 'number of alignments / isoforms'
        name_fig = 'blocks_per_alignment'

        blocks_per_alignment_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, isoforms_label, isoforms_distribution,
                               y_log_scale=True, def_step=1)

        blocks_per_alignment_plot.plot_compare_distribution(self.short_report_plots)


    # BASIC AND SIMPLE 5
    # Multiple-aligned transcripts distribution
    # Number of unaligned / aligned / alignments / unique aligned / multiple-aligned transcripts
    def get_alignment_multiplicity_plot(self, transcripts_metrics, out_dir, short_report_visible):
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

                transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'alignment multiplicity'
        label_x = 'number of alignments'
        label_y = 'number of transcripts'
        name_fig = 'alignment_multiplicity'

        alignment_multiplicity_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, y_log_scale=True, def_step=1)

        alignment_multiplicity_plot.plot_compare_distribution(self.short_report_plots)


    # BASIC AND SIMPLE 6
    # Mismatch number per transcript distribution (tot, avg)
    def get_mismatch_rate_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        transcripts_labels = []
        transcripts_distributions = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics is not None:
                transcripts_distributions.append(transcripts_metrics[i_transcripts].simple_metrics.mismatch_num_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_mismatch_num, precision)),
                              ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_mismatch_num)]

                transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'substitution errors per alignment'
        label_x = 'number of mismatches'
        label_y = 'number of alignments'
        name_fig = 'mismatch_rate'

        caption = \
            'Plot showing cumulative substitution errors per alignment distribution. Each point represents the ' \
            'number of alignments with the corresponding number of mismatches or greater; ' \
            'the plot is given in logarithmic scale.'
        mismatch_rate_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, y_log_scale=True, caption=caption,
                               def_step=1)

        mismatch_rate_plot.plot_compare_distribution(self.short_report_plots)


    '''# BASIC AND SIMPLE 7
    # Target gap number per transcript distribution (tot, avg, min, max)
    # Number of introns per isoform distribution  (tot, avg, min, max)
    def get_basic_simple_7(self, basic_isoforms_metrics, transcripts_metrics, out_dir, precision, short_report_visible):
        isoforms_distribution = None
        isoforms_label = None
        if basic_isoforms_metrics is not None:
            isoforms_distribution = basic_isoforms_metrics.introns_num_distribution

            str = 'Number of introns per isoform'
            list_label = [('min', basic_isoforms_metrics.min_introns_num),
                          ('max', basic_isoforms_metrics.max_introns_num),
                          ('avg', round(basic_isoforms_metrics.avg_introns_num, precision)),
                          ('tot', basic_isoforms_metrics.tot_introns_num)]
            isoforms_label = DistributionReport.get_label(str, list_label)

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
                t_label_i_0 = DistributionReport.get_label(str, list_label)

                transcripts_labels.append([t_label_i_0])

        title_name = 'Basic and Simple 7'
        label_x = 'number of gaps/introns'
        label_y = 'number of transcripts/isoforms'
        name_fig = '7'
        plot_7 = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, isoforms_label, isoforms_distribution,
                               y_log_scale=True, def_step=1)

        plot_7.plot_compare_distribution(self.short_report_plots)'''


    '''# BASIC AND SIMPLE 9
    # Query gap number per transcript distribution (tot, avg, min, max)
    def get_gaps_per_alignment_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics is not None:
                transcripts_distribution.append(transcripts_metrics[i_transcripts].simple_metrics.qgap_num_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_qgap_num, precision)),
                              ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_qgap_num)]

                transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'number of gaps per alignment'
        label_x = 'number of gaps'
        label_y = 'number of alignments'
        name_fig = 'gaps_per_alignment'
        gaps_per_alignment_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distribution, y_log_scale=True, def_step=1)

        gaps_per_alignment_plot.plot_compare_histogram(transcripts_metrics, self.short_report_plots)'''


    '''# BASIC AND SIMPLE 10
    # Query gap len distribution (tot, avg, min, max)
    def get_gap_length_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            if transcripts_metrics[i_transcripts].simple_metrics is not None:
                transcripts_distribution.append(transcripts_metrics[i_transcripts].simple_metrics.qgap_len_distribution)

                str = transcripts_metrics[i_transcripts].label
                list_label = [('avg', round(transcripts_metrics[i_transcripts].simple_metrics.avg_qgap_len, precision)),
                              ('tot', transcripts_metrics[i_transcripts].simple_metrics.tot_qgap_len)]

                transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'alignment gaps length'
        label_x = 'gap length'
        label_y = 'number of alignments'
        name_fig = 'gap_length'

        gap_len_plot = UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir, transcripts_labels, transcripts_distribution, y_log_scale=True)
        gap_len_plot.plot_compare_histogram(transcripts_metrics, self.short_report_plots)'''


    def get_x_matched_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        transcripts_labels = []
        transcripts_distributions = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distributions.append(transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage.covered_fraction_whole_transcript_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage.avg_covered_fraction_whole_transcript, precision))]

            transcripts_labels.append(DistributionReport.get_label(str, list_label))

        caption = 'Plot showing cumulative transcript match histogram. Each bar represents the number of transcripts ' \
                  'with matched fraction equal to or greater than the value on $x$ axis; transcript matched fraction is ' \
                  'calculated as the number of its bases covering an isoform divided by the transcript length.'
        title_name = 'transcript matched fraction'
        label_x = 'matched fraction'
        label_y = 'number of transcripts'
        name_fig = 'x-matched'

        x_matched_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, caption=caption, def_step=0.1)

        x_matched_plot.plot_compare_histogram(transcripts_metrics, self.short_report_plots)


    def get_x_matched_blocks_plot(self, transcripts_metrics, out_dir, precision, short_report_visible=False):
        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage.block_fraction_matched_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage.avg_covered_fraction_block, precision))]

            transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'block matched fraction'
        label_x = 'matched fraction'
        label_y = 'number of blocks'
        name_fig = 'x-matched_blocks'

        x_matched_blocks_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distribution, def_step=0.1)

        x_matched_blocks_plot.plot_compare_histogram(transcripts_metrics, self.short_report_plots)


    def get_specificity_plots(self, transcripts_metrics, out_dir, precision):
        # DISTRIBUTIONS FOR TRANSCRIPTS ALIGNMENTS
        if transcripts_metrics[0].assembly_correctness_metrics.transcripts_coverage is not None:
            self.get_x_matched_plot(transcripts_metrics, out_dir, precision, short_report_visible=True)

            self.get_x_matched_blocks_plot(transcripts_metrics, out_dir, precision, short_report_visible=False)


    # CONSIDER COVERED BASES BY EACH MAPPED TRANSCRIPT SEPARATELY:
    # distributions for isoforms:
    def get_x_assembled_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.assembled_fraction_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_assembled_fraction, precision))]

            transcripts_labels.append(DistributionReport.get_label(str, list_label))

        caption = 'Plot showing cumulative isoform assembly histogram. Each bar represents the number of isoforms ' \
                  'with assembled fraction equal to or greater than the value on $x$ axis; isoform assembled fraction is ' \
                  'calculated as the maximum number of captured by single assembled transcript bases divided by the ' \
                  'total isoform length.'
        title_name = 'isoform assembled fraction'
        label_x = 'assembled fraction'
        label_y = 'number of isoforms'
        name_fig = 'x-assembled'

        x_assembled_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distribution, caption=caption, def_step=0.1)

        x_assembled_plot.plot_compare_histogram(transcripts_metrics, self.short_report_plots)


    def get_x_assembled_exons_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.assembled_fraction_exons_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_assembled_fraction_exons, precision))]

            transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'exon assembled fraction'
        label_x = 'assembled fraction'
        label_y = 'number of exons'
        name_fig = 'x-assembled_exons'

        x_assembled_exons_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distribution, def_step=0.1)

        x_assembled_exons_plot.plot_compare_histogram(transcripts_metrics, self.short_report_plots)


    # CONSIDER COVERED BASES (BY ALL MAPPED TRANSCRIPTS) WITHOUT OVERLAPS:
    # distributions for isoforms:
    def get_x_covered_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.covered_fraction_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_covered_fraction, precision))]

            transcripts_labels.append(DistributionReport.get_label(str, list_label))

        caption = 'Plot showing cumulative isoform coverage histogram. Each bar represents the number of isoforms ' \
                  'with covered fraction equal to or greater than the value on $x$ axis; isoform covered fraction is ' \
                  'calculated as the number of covered bases (by all transcripts in the assembly) divided by the ' \
                  'total isoform length.'
        title_str = 'isoform covered fraction'
        label_x = 'covered fraction'
        label_y = 'number of isoforms'
        name_fig = 'x-covered'

        x_covered_plot = \
            UtilsPictures.Plot(title_str, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distribution, caption=caption, def_step=0.1)

        x_covered_plot.plot_compare_histogram(transcripts_metrics, self.short_report_plots)


    # distributions for exons:
    def get_x_covered_exons_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        transcripts_labels = []
        transcripts_distribution = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.covered_fraction_exons_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_covered_fraction_exons, precision))]

            transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_name = 'exon covered fraction'
        label_x = 'covered fraction'
        label_y = 'number of exons'
        name_fig = 'x-covered_exons'

        x_covered_exons_plot = \
            UtilsPictures.Plot(title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distribution, def_step=0.1)

        x_covered_exons_plot.plot_compare_histogram(transcripts_metrics, self.short_report_plots)


    # def get_avg_exon_covered_fraction_over_isoforms_distribution(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        # title_name = 'average fraction of exon in isoform covered'
        # label_x = 'avg fraction of exon covered'
        # label_y = 'number of isoforms'
        # name_fig = 'avg_exon_in_isoform_fraction_covered'
        # step = 0.1
        #
        # transcripts_labels = []
        # transcripts_distribution = []
        # for i_transcripts in range(len(transcripts_metrics)):
        #     transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_covered_fraction_exon_distribution)
        #
        #     str = transcripts_metrics[i_transcripts].label
        #     list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_covered_fraction_exon, precision))]
        #
        #     transcripts_labels.append(UtilsPictures.get_label(str, list_label))
        #
        # UtilsPictures.plot_compare_histogram(out_dir, title_name, transcripts_distribution, transcripts_metrics, transcripts_labels,
        #                                      label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)

    # def get_percentage_well_covered_exons_over_isoforms_distribution(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        # title_name = 'percentage of well-covered exons in isoform '
        # label_x = 'percentage of isoform well-covered exons'
        # label_y = 'number of isoforms'
        # name_fig = 'well_covered_exons_in_isoform'
        # step = 0.1
        #
        # transcripts_labels = []
        # transcripts_distribution = []
        # for i_transcripts in range(len(transcripts_metrics)):
        #     transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.percentage_well_covered_exons_distribution)
        #
        #     str = transcripts_metrics[i_transcripts].label
        #     list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_percentage_well_covered_exons, precision))]
        #
        #     transcripts_labels.append(UtilsPictures.get_label(str, list_label))
        #
        # UtilsPictures.plot_compare_histogram(out_dir, title_name, transcripts_distribution, transcripts_metrics, transcripts_labels,
        #                                      label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)

    # def get_percentage_fully_covered_exons_over_isoforms_distribution(self, transcripts_metrics, out_dir, precision, short_report_visible):
        # distribution for transcripts by all mapped transcripts at once
        # title_name = 'percentage of fully-covered exons in isoform'
        # label_x = 'percentage of isoform fully-covered exons'
        # label_y = 'number of isoforms'
        # name_fig = 'fully_covered_exons_in_isoform'
        # step = 0.1
        #
        # transcripts_labels = []
        # transcripts_distribution = []
        # for i_transcripts in range(len(transcripts_metrics)):
        #     transcripts_distribution.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.percentage_fully_covered_exons_distribution)
        #
        #     str = transcripts_metrics[i_transcripts].label
        #     list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_percentage_fully_covered_exons, precision))]
        #
        #     transcripts_labels.append(UtilsPictures.get_label(str, list_label))
        #
        # UtilsPictures.plot_compare_histogram(out_dir, title_name, transcripts_distribution, transcripts_metrics, transcripts_labels,
        #                                      label_x, label_y, name_fig, short_report_visible, self.pdf_plots_figures, def_step=step)


    def get_alignments_per_isoform_plot(self, transcripts_metrics, out_dir, precision, short_report_visible):
        transcripts_labels = []
        transcripts_distributions = []
        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_distributions.append(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.num_transcripts_mapped_to_isoform_distribution)

            str = transcripts_metrics[i_transcripts].label
            list_label = [('avg', round(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage.avg_num_transcripts_mapped_to_isoform, precision))]

            transcripts_labels.append(DistributionReport.get_label(str, list_label))

        title_str = 'number of alignments per isoform'
        label_x = 'number of alignments'
        label_y = 'number of isoforms'
        name_fig = 'alignments_per_isoform'

        alignments_per_isoform_plot = \
            UtilsPictures.Plot(title_str, label_x, label_y, name_fig, short_report_visible, out_dir,
                               transcripts_labels, transcripts_distributions, y_log_scale=True, def_step=1)

        alignments_per_isoform_plot.plot_compare_distribution(self.short_report_plots)


    def get_sensitivity_plots(self, transcripts_metrics, out_dir, precision):
        if transcripts_metrics[0].assembly_completeness_metrics.isoforms_coverage is not None:
            self.get_x_assembled_plot(transcripts_metrics, out_dir, precision, short_report_visible=True)

            self.get_x_assembled_exons_plot(transcripts_metrics, out_dir, precision, short_report_visible=False)

            self.get_x_covered_plot(transcripts_metrics, out_dir, precision, short_report_visible=True)

            self.get_x_covered_exons_plot(transcripts_metrics, out_dir, precision, short_report_visible=False)

            self.get_alignments_per_isoform_plot(transcripts_metrics, out_dir, precision, short_report_visible=False)

            # self.get_percentage_well_covered_exons_over_isoforms_distribution(transcripts_metrics, out_dir, precision, short_report_visible=False)

            # self.get_percentage_fully_covered_exons_over_isoforms_distribution(transcripts_metrics, out_dir, precision, short_report_visible=False)


    @classmethod
    def get_label(cls, str, list_labels):
        label = '{}:\n'.format(str)

        for i_label in range(len(list_labels)):
            (key, value) = list_labels[i_label]
            label += '{}={} '.format(key, value)
            if i_label % 3 == 2 and i_label != len(list_labels) - 1:
                label += '\n'

        return label