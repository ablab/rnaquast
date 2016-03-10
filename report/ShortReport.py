__author__ = 'letovesnoi'

import os

from quast23.libs import reporting

# Font of plot captions, axes labels and ticks
font = {'family': 'sans-serif', 'style': 'normal', 'weight': 'medium', 'size': 10}


class ShortReport():
    """Class which generate short report"""

    @classmethod
    def get_metrics_types(cls):
        metrics_types = ['DATABASE METRICS',
                         'BASIC TRANSCRIPTS METRICS',
                         'ALIGNMENT METRICS',
                         'ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS',
                         'ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS',
                         'ASSEMBLY COMPLETENESS (SENSITIVITY)',
                         'BUSCO METRICS',
                         'GeneMarkS-T METRICS',
                         'ASSEMBLY SPECIFICITY']
        return metrics_types


    @classmethod
    def get_metrics_labels(cls, TRANSCRIPT_LENS, WELL_FULLY_COVERAGE_THRESHOLDS):
        metrics_labels = ['Genes',
                          'Avg. number of exons per isoform',
                          'Transcripts',
                          'Transcripts > {} bp'.format(str(TRANSCRIPT_LENS[0])),
                          'Transcripts > {} bp'.format(str(TRANSCRIPT_LENS[1])),
                          'Aligned',
                          'Uniquely aligned',
                          'Multiply aligned',
                          'Unaligned',
                          'Avg. aligned fraction',
                          'Avg. alignment length',
                          'Avg. mismatches per transcript',
                          'Misassemblies',
                          'Database coverage',
                          'Relative database coverage',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled genes',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled genes',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered genes',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered genes',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled isoforms',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled isoforms',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered isoforms',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered isoforms',
                          'Mean isoform coverage',
                          'Mean isoform assembly',
                          'Complete',
                          'Partial',
                          'Predicted genes',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_transcript_threshold * 100)) + '%-matched',
                          str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_transcript_threshold * 100)) + '%-matched',
                          'Unannotated',
                          'Mean fraction of transcript matched']
        return metrics_labels


    @classmethod
    def get_metrics_type_labels_dict(cls, metrics_types, metrics_labels):
        metrics_type_labels_dict = \
            {metrics_types[0]: metrics_labels[0:2],
             metrics_types[1]: metrics_labels[2:5],
             metrics_types[2]: metrics_labels[5:9],
             metrics_types[3]: metrics_labels[9:12],
             metrics_types[4]: metrics_labels[12],
             metrics_types[5]: metrics_labels[13:25],
             metrics_types[6]: metrics_labels[25:27],
             metrics_types[7]: metrics_labels[27],
             metrics_types[8]: metrics_labels[28:]}
        return metrics_type_labels_dict


    @classmethod
    def get_best_type(cls, metrics_labels):
        # absent if there are no best values
        # 1 if the max value is the best value
        # -1 if the min value is the best value
        # 2 if the max relative value (i.e. divided by the total number of transcripts in the corresponding assembly) is the best value
        # -2 if the min relative value (i.e. divided by the total number of transcripts in the corresponding assembly) is the best value
        best_type = {metrics_labels[3]: 2, metrics_labels[4]: 2,

                     metrics_labels[5]: 2, metrics_labels[8]: -2,

                     metrics_labels[9]: 1, metrics_labels[10]: 1, metrics_labels[11]: -1,

                     metrics_labels[12]: -2,

                     metrics_labels[13]: 1, metrics_labels[14]: 1, metrics_labels[15]: 1, metrics_labels[16]: 1,
                     metrics_labels[17]: 1, metrics_labels[18]: 1, metrics_labels[19]: 1, metrics_labels[20]: 1,
                     metrics_labels[21]: 1, metrics_labels[22]: 1, metrics_labels[23]: 1, metrics_labels[24]: 1,

                     metrics_labels[25]: 1, metrics_labels[26]: 1,

                     metrics_labels[27]: 1,

                     metrics_labels[28]: 2, metrics_labels[29]: 2, metrics_labels[30]: -2, metrics_labels[31]: 1}

        return best_type

    @classmethod
    def get_i_rel_best_metrics(cls, metrics_labels, best_type):
        i_rel_best_metrics = []
        for i_metric_label in range(len(metrics_labels)):
            if metrics_labels[i_metric_label] in best_type and \
                    (best_type[metrics_labels[i_metric_label]] == 2 or best_type[metrics_labels[i_metric_label]] == -2):
                i_rel_best_metrics.append(i_metric_label + 1)
        return i_rel_best_metrics


    @classmethod
    def get_best_values(self, metrics_dict, best_type):
        best_values = {}

        transcripts_num = [int(num_str) for num_str in metrics_dict['Transcripts']]

        for metric_label in metrics_dict.keys():
            values = metrics_dict[metric_label]
            if metric_label not in best_type:
                best_values[metric_label] = None
            else:
                if best_type[metric_label] == 1:
                    max_value = - float('Inf')
                    for i_v in range(len(values)):
                        float_value = float(values[i_v])
                        if float_value > max_value:
                            max_value = float_value
                            argmax = values[i_v]
                    best_values[metric_label] = argmax
                if best_type[metric_label] == -1:
                    min_value = float('Inf')
                    for i_v in range(len(values)):
                        float_value = float(values[i_v])
                        if float_value < min_value:
                            min_value = float_value
                            argmin = values[i_v]
                    best_values[metric_label] = argmin
                if best_type[metric_label] == 2:
                    max_rel_value = - float('Inf')
                    for i_v in range(len(values)):
                        rel_value = float(values[i_v]) / transcripts_num[i_v]
                        if rel_value > max_rel_value:
                            max_rel_value = rel_value
                            argmax = values[i_v]
                    best_values[metric_label] = argmax
                if best_type[metric_label] == -2:
                    min_rel_value = float('Inf')
                    for i_v in range(len(values)):
                        rel_value = float(values[i_v]) / transcripts_num[i_v]
                        if rel_value < min_rel_value:
                            min_rel_value = rel_value
                            argmin = values[i_v]
                    best_values[metric_label] = argmin

        return best_values


    def __init__(self, args, db_genes_metrics, transcripts_metrics, outdir, separated_reports,
                 comparison_report, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS):
        self.name = 'short_report'

        # output directory:
        self.out_dir = outdir

        self.first_label = 'METRICS/TRANSCRIPTS'
        self.metrics_type = ShortReport.get_metrics_types()
        self.metrics_labels = ShortReport.get_metrics_labels(TRANSCRIPT_LENS, WELL_FULLY_COVERAGE_THRESHOLDS)
        self.metrics_type_labels_dict = ShortReport.get_metrics_type_labels_dict(self.metrics_type, self.metrics_labels)

        self.best_type = ShortReport.get_best_type(self.metrics_labels)

        # txt file:
        self.path_txt = os.path.join(outdir, '{}.txt'.format(self.name))
        # tab separated file:
        self.path_tsv = os.path.join(outdir, '{}.tsv'.format(self.name))
        # tex file:
        self.path_tex = os.path.join(outdir, '{}.tex'.format(self.name))
        # pdf file:
        self.path_pdf = os.path.join(outdir, '{}.pdf'.format(self.name))

        # get table with all short report metrics:
        self.metrics_table = self.set_metrics_table(args, db_genes_metrics, transcripts_metrics,
                                                    WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS)

        column_n = len(transcripts_metrics) + 1

        self.metrics_dict = self.get_metrics_dict(column_n)

        # if there are no transcripts, add one row for annotation metrics:
        # if len(transcripts_metrics) == 0:
        #     column_n = 2

        self.get_report(args, column_n, separated_reports, comparison_report, logger)


    def get_report(self, args, column_n, separated_reports, comparison_report, logger):
        logger.print_timestamp()
        logger.info('Getting SHORT SUMMARY report...')

        self.column_widths = self.get_column_widths(column_n)

        self.best_values = self.get_best_values(self.metrics_dict, self.best_type)

        distribution_report = None
        if not args.no_plots:
            # for several files with transcripts select only comparison report pictures:
            if comparison_report is not None:
                distribution_report = comparison_report.distribution_report
            else:
                distribution_report = separated_reports[0].distribution_report

        # TXT:
        self.print_txt()

        # TSV:
        #self.print_tsv()

        # TEX:
        self.print_tex(column_n, distribution_report)

        # PDF:
        #self.print_pdf(args, separated_reports, comparison_report, logger)

        logger.info('  saved to\n' + '    ' + '{}\n'.format(self.path_txt) + '    ' + '{}\n'.format(self.path_tex) +
                    4 * ' ' + '{}'.format(self.path_pdf))


    def get_column_widths(self, column_n):
        widths = []
        widths.append(50)
        for i in range(1, column_n):
            widths.append(25)
        return widths


    def get_metrics_dict(self, column_n):
        metrics_dict = {}
        for line in self.metrics_table:
            arr = line.split('  ')
            new_arr = []
            for substr in arr:
                if substr.strip() != '':
                    new_arr.append(substr.strip())
            if len(new_arr) == column_n:
                metrics_dict[new_arr[0]] = new_arr[1:]
        return metrics_dict


    def print_txt(self):
        fout_txt_file = open(self.path_txt, 'w')
        fout_txt_file.write('SHORT REPORT \n')
        for line in self.metrics_table:
            fout_txt_file.write(line)
        fout_txt_file.close()


    def print_tsv(self):
        fout_tsv_file = open(self.path_tsv, 'w')
        for line in self.table_to_draw:
            fout_tsv_file.write('\t'.join(line) + '\n')


    def print_tex(self, column_n, distribution_report):
        fout_tex_file = open(self.path_tex, 'w')

        print >> fout_tex_file, '\\documentclass[12pt,a4paper]{article}\n'

        print >> fout_tex_file, '\usepackage{fancyhdr}'
        print >> fout_tex_file, '\usepackage{graphicx}'
        print >> fout_tex_file, '\n'

        print >> fout_tex_file, '\\begin{document}\n'

        print >> fout_tex_file, '\pagestyle{fancy}'
        print >> fout_tex_file, r'\fancyhf{}'
        print >> fout_tex_file, '\chead{Short summary report}\n'

        print >> fout_tex_file, '\\begin{table}[!ht]'
        print >> fout_tex_file, '\\begin{center}'
        print >> fout_tex_file, '\small'

        i_rel_best_metrics = ShortReport.get_i_rel_best_metrics(self.metrics_labels, self.best_type)
        print >> fout_tex_file, \
            r'\caption {rnaQUAST metrics for assembled transcripts. In each row the best values are indicated with ' \
            r'\textbf{bold}. For the transcript metrics (rows ' + str(i_rel_best_metrics)[1:-1] + \
            r') we highlighted the best \textbf{relative} values i.e. divided by the total number of transcripts in ' \
            r'the corresponding assembly.}'

        print >> fout_tex_file, '\\begin{tabular}{|l*{' + reporting.val_to_str(column_n) + '}{|r}|}'
        print >> fout_tex_file, '\\hline'

        tex_str = '{:<50}'.format(r'\textbf{' + self.first_label + '}')
        for t_label in self.metrics_dict[self.first_label]:
            tex_str += '{:<25}'.format(' & ' + r'\textbf{' + t_label + '}')
        tex_str += 10 * ' ' + r'\\ \hline\hline' + '\n'
        for metric_type in self.metrics_type:
            type_flag = False
            tmp_tex_str_type = '{:<50}'.format(r'\multicolumn{' + str(column_n) + r'}{l}{\bf ' + metric_type + '}') + 10 * ' ' + r'\\ \hline' + '\n'
            tmp_tex_str_label = ''
            for i_metric_label in range(len(self.metrics_type_labels_dict[metric_type])):
                metric_label = self.metrics_type_labels_dict[metric_type][i_metric_label]
                if metric_label in self.metrics_dict:
                    type_flag = True
                    tmp_tex_str_label += '{:<50}'.format(metric_label.replace('>', '$>$').replace('<', '$<$').replace('%', '\%'))
                    for i_metric_value in range(len(self.metrics_dict[metric_label])):
                        metric_value = self.metrics_dict[metric_label][i_metric_value]

                        if metric_value == self.best_values[metric_label]:
                            metric_value = r'\textbf{' + metric_value + '}'

                        tmp_tex_str_label += '{:<25}'.format(' & ' + metric_value)
                    tmp_tex_str_label += r' \\' + '\n'
            if type_flag:
                tex_str += tmp_tex_str_type + tmp_tex_str_label.strip() + r' \hline' + '\n'
        fout_tex_file.write(tex_str)

        print >> fout_tex_file, '\\end{tabular}'
        print >> fout_tex_file, '\\end{center}'
        print >> fout_tex_file, '\\end{table}\n'

        print >> fout_tex_file, '\lfoot{generated by rnaQUAST}\n'


        # FIGURES:
        if distribution_report is not None:
            short_report_plots = distribution_report.short_report_plots

            for plot in short_report_plots:
                self.add_figure_to_tex(fout_tex_file, plot)

        print >> fout_tex_file, '\\end{document}'

        fout_tex_file.close()


    def add_figure_to_tex(self, fout_tex_file, plot):
        print >> fout_tex_file, r'\begin{figure}[!ht]'
        print >> fout_tex_file, r'\begin{center}'
        print >> fout_tex_file, r'\includegraphics[width = \linewidth]{' + plot.path + '}'
        print >> fout_tex_file, '\end{center}'
        print >> fout_tex_file, '\caption{' + plot.caption + '}'
        print >> fout_tex_file, '\end{figure}'
        print >> fout_tex_file, '\n'


    # full report in PDF format: all tables and plots
    def print_pdf(self, args, separated_reports, comparison_report, logger):
        pdf_tables_figures = [self.get_pdf_table_figure('Short report', 'generated by rnaQUAST', self.table_to_draw, self.column_widths, logger)]

        all_pdf_file = None
        if not args.no_plots:
            from quast23.libs import plotter  # Do not remove this line! It would lead to a warning in matplotlib.
            try:
                from matplotlib.backends.backend_pdf import PdfPages
                all_pdf_file = PdfPages(self.path_pdf)
            except:
                all_pdf_file = None
        if all_pdf_file:
            # for several files with transcripts select only comparison report pictures:
            if comparison_report is not None:
                pdf_plots_figures = comparison_report.distribution_report.pdf_plots_figures
            else:
                pdf_plots_figures = separated_reports[0].distribution_report.pdf_plots_figures

            self.fill_all_pdf_file(all_pdf_file, pdf_tables_figures, pdf_plots_figures, logger)


    # draw_report_table from quast23.libs.plotter:
    def get_pdf_table_figure(self, report_name, extra_info, table_to_draw, column_widths, logger):
        # checking if matplotlib and pylab is installed:
        matplotlib_error = False
        try:
            import matplotlib
            matplotlib.use('Agg')  # non-GUI backend
            if matplotlib.__version__.startswith('0'):
                logger.warning('matplotlib version is rather old! Please use matplotlib version 1.0 or higher for better results.')
            #from pylab import *
        except Exception:
            logger.warning('Can\'t draw plots: please install python-matplotlib and pylab.')
            matplotlib_error = True

        if matplotlib_error:
            return

        # some magic constants ..
        font_size = 12
        font_scale = 2
        external_font_scale = 10
        letter_height_coeff = 0.10
        letter_width_coeff = 0.04

        row_height = letter_height_coeff * font_scale
        nrows = len(table_to_draw)
        external_text_height = float(font["size"] * letter_height_coeff * external_font_scale) / font_size
        total_height = nrows * row_height + 2 * external_text_height
        total_width = letter_width_coeff * font_scale * sum(column_widths)

        import matplotlib.pyplot
        figure = matplotlib.pyplot.figure(figsize=(total_width, total_height))
        matplotlib.pyplot.rc('font', **font)
        matplotlib.pyplot.axis('off')
        matplotlib.pyplot.text(0.5 - float(column_widths[0]) / (2 * sum(column_widths)),
                               1. - float(2 * row_height) / total_height, report_name.replace('_', ' ').capitalize())
        matplotlib.pyplot.text(0 - float(column_widths[0]) / (2 * sum(column_widths)), 0, extra_info)

        #with_bold_table = ShortReport.get_with_bold_table(table_to_draw, self.best_values)

        colLabels=table_to_draw[0][1:]
        rowLabels=[item[0] for item in table_to_draw[1:]]
        restValues=[item[1:] for item in table_to_draw[1:]]

        matplotlib.pyplot.table(cellText=restValues, rowLabels=rowLabels, colLabels=colLabels,
            colWidths=[float(column_width) / sum(column_widths) for column_width in column_widths[1:]],
            rowLoc='left', colLoc='center', cellLoc='right', loc='center')
        return figure


    def fill_all_pdf_file(self, all_pdf, pdf_tables_figures, pdf_plots_figures, logger):
        # checking if matplotlib and pylab is installed:
        matplotlib_error = False
        try:
            import matplotlib
            matplotlib.use('Agg')  # non-GUI backend
            if matplotlib.__version__.startswith('0'):
                logger.warning('matplotlib version is rather old! Please use matplotlib version 1.0 or higher for better results.')
            #from pylab import *
        except Exception:
            logger.warning('Can\'t draw plots: please install python-matplotlib and pylab.')
            matplotlib_error = True

        if matplotlib_error or not all_pdf:
            return

        for figure in pdf_tables_figures:
            all_pdf.savefig(figure, bbox_inches='tight')
        for figure in pdf_plots_figures:
            all_pdf.savefig(figure, additional_artists='art', bbox_inches='tight')

        try:  # for matplotlib < v.1.0
            d = all_pdf.infodict()
            d['Title'] = 'rnaQUAST short report'
            d['Author'] = 'rnaQUAST'
            import datetime
            d['CreationDate'] = datetime.datetime.now()
            d['ModDate'] = datetime.datetime.now()
        except AttributeError:
            pass
        all_pdf.close()


    # generate table with all metrics for short report:
    def set_metrics_table(self, args, db_genes_metrics, transcripts_metrics, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS):
        self.metrics_table = []

        name_str = '{:<50}'.format('METRICS/TRANSCRIPTS')
        for i_transcripts in range(len(transcripts_metrics)):
            name_str += '{:<25}'.format(transcripts_metrics[i_transcripts].label)
        self.metrics_table.append(name_str + '\n')

        if db_genes_metrics is not None:
            # ======= DATABASE METRICS ===========
            self.add_database_metrics_to_table(db_genes_metrics, transcripts_metrics, PRECISION)

        if len(transcripts_metrics) >= 1:
            # ======= BASIC TRANSCRIPTS METRICS ===========
            if transcripts_metrics[0].basic_metrics is not None:
                self.add_basic_metrics_to_table(transcripts_metrics, PRECISION, TRANSCRIPT_LENS)

            # ======= ALIGNMENT METRICS ===========
            if transcripts_metrics[0].simple_metrics is not None:
                self.add_alignment_metrics_to_table(transcripts_metrics, PRECISION)

                self.add_fusion_misassemble_metrics_to_table(transcripts_metrics)

            # ======= ASSEMBLY COMPLETENESS METRICS ===========
            if transcripts_metrics[0].assembly_completeness_metrics is not None:
                self.add_assemble_completeness_metrics_to_table(transcripts_metrics, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

            # ======= ASSEMBLY CORRECTNESS METRICS ===========
            if transcripts_metrics[0].assembly_correctness_metrics is not None:
                    self.add_assemble_correctness_metrics_to_table(transcripts_metrics, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION)

        return self.metrics_table


    def add_database_metrics_to_table(self, db_genes_metrics, transcripts_metrics, PRECISION):
        # Database short report metrics:
        tot_genes_num_str = '{:<50}'.format('Genes')
        avg_exons_num_str = '{:<50}'.format('Avg. number of exons per isoform')

        # tot_protein_coding_genes_num_str = '{:<50}'.format('Protein coding genes')

        # tot_isoforms_num_str = '{:<50}'.format('Isoforms')
        # tot_protein_coding_isoforms_num_str = '{:<50}'.format('Protein coding isoforms')
        # avg_len_wout_introns_str = '{:<50}'.format('Avg. length of all isoforms')
        # avg_exon_len_str = '{:<50}'.format('Avg. exon length')

        num_row = len(transcripts_metrics)
        if len(transcripts_metrics) == 0:
            num_row = 1

        for i_transcript in range(num_row):
            if db_genes_metrics is not None:
                tot_genes_num_str += '{:<25}'.format(db_genes_metrics.genes_num)

                avg_exons_num_str += '{:<25}'.format(round(db_genes_metrics.avg_exons_num, PRECISION))

        # tot_protein_coding_genes_num_str += '{:<25}'.format(basic_genes_metrics.num_protein_coding)

        # tot_isoforms_num_str += '{:<25}'.format(basic_isoforms_metrics.number)
        # tot_protein_coding_isoforms_num_str += '{:<25}'.format(basic_isoforms_metrics.num_protein_coding)
        # avg_len_wout_introns_str += '{:<25}'.format(round(basic_isoforms_metrics.avg_len_wout_introns, precision))
        # avg_exon_len_str += '{:<25}'.format(round(basic_isoforms_metrics.avg_exon_len, precision))

        self.metrics_table.append(' == DATABASE METRICS == \n')
        if db_genes_metrics is not None:
            self.metrics_table.append(tot_genes_num_str + '\n')
            self.metrics_table.append(avg_exons_num_str + '\n')

        # gff source field unfortunately contain database name instead of transcripts type protein_coding:
        # if basic_genes_metrics.num_protein_coding != 0:
        #     self.metrics_table.append(tot_protein_coding_genes_num_str + '\n')

        # self.metrics_table.append(tot_isoforms_num_str + '\n')

        # gff source field unfortunately contain database name instead of transcripts type protein_coding:
        # if basic_isoforms_metrics.num_protein_coding != 0:
        #     self.metrics_table.append(tot_protein_coding_isoforms_num_str + '\n')

        # self.metrics_table.append(avg_len_wout_introns_str + '\n')
        # self.metrics_table.append(avg_exon_len_str + '\n')


    def add_basic_metrics_to_table(self, transcripts_metrics, PRECISION, TRANSCRIPT_LENS):
        # Basic short report metrics:
        num_transcripts_str = '{:<50}'.format('Transcripts')
        num_transcripts_500_str = '{:<50}'.format('Transcripts > {} bp'.format(str(TRANSCRIPT_LENS[0])))
        num_transcripts_1000_str = '{:<50}'.format('Transcripts > {} bp'.format(str(TRANSCRIPT_LENS[1])))

        # avg_transcript_len_str = '{:<50}'.format('Avg. length')
        # n50_str = '{:<50}'.format('N50')
        # max_transcript_len_str = '{:<50}'.format('Longest transcript')

        for i_transcripts in range(len(transcripts_metrics)):
            basic_metrics = transcripts_metrics[i_transcripts].basic_metrics

            num_transcripts_str += '{:<25}'.format(basic_metrics.number)
            num_transcripts_500_str += '{:<25}'.format(basic_metrics.num_transcripts_500)
            num_transcripts_1000_str += '{:<25}'.format(basic_metrics.num_transcripts_1000)

            # avg_transcript_len_str += '{:<25}'.format(round(basic_transcripts_metrics.avg_len, precision))
            # n50_str += '{:<25}'.format(basic_transcripts_metrics.n50)
            # max_transcript_len_str += '{:<25}'.format(basic_transcripts_metrics.max_len)

        self.metrics_table.append('\n == BASIC TRANSCRIPTS METRICS == \n')
        self.metrics_table.append(num_transcripts_str + '\n')
        self.metrics_table.append(num_transcripts_500_str + '\n')
        self.metrics_table.append(num_transcripts_1000_str + '\n')

        # self.metrics_table.append(avg_transcript_len_str + '\n')
        # self.metrics_table.append(n50_str + '\n')
        # self.metrics_table.append(max_transcript_len_str + '\n')


    def add_alignment_metrics_to_table(self, transcripts_metrics, PRECISION):
        # Alignment short report metrics:
        num_aligned_str = '{:<50}'.format('Aligned')
        num_uniquely_aligned_str = '{:<50}'.format('Uniquely aligned')
        num_multiply_aligned_str = '{:<50}'.format('Multiply aligned')
        num_unaligned_str = '{:<50}'.format('Unaligned')

        for i_transcripts in range(len(transcripts_metrics)):
            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics

            num_aligned_str += '{:<25}'.format(simple_metrics.num_aligned)
            num_uniquely_aligned_str += '{:<25}'.format(simple_metrics.num_unique_aligned)
            num_multiply_aligned_str += '{:<25}'.format(simple_metrics.num_mul_aligned)
            num_unaligned_str += '{:<25}'.format(simple_metrics.num_unaligned)

            # num_alignments_str += '{:<25}'.format(simple_metrics.num_alignments)

        self.metrics_table.append('\n == ALIGNMENT METRICS == \n')
        self.metrics_table.append(num_aligned_str + '\n')
        self.metrics_table.append(num_uniquely_aligned_str + '\n')
        self.metrics_table.append(num_multiply_aligned_str + '\n')
        self.metrics_table.append(num_unaligned_str + '\n')
        # self.metrics_table.append(num_alignments_str + '\n')

        # OTHER SECTION:
        avg_aligned_fraction_str = '{:<50}'.format('Avg. aligned fraction')
        avg_alignment_len_str = '{:<50}'.format('Avg. alignment length')
        avg_mismatches_per_transcript_str = '{:<50}'.format('Avg. mismatches per transcript')

        # na50_str = '{:<50}'.format('NA50')
        # avg_blocks_num_str = '{:<50}'.format('Avg. blocks per alignment')
        # avg_block_len_str = '{:<50}'.format('Avg. block length')
        # avg_mismatches_str = '{:<50}'.format('Avg. mismatches per transcript')

        for i_transcripts in range(len(transcripts_metrics)):
            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics

            avg_aligned_fraction_str += '{:<25}'.format(round(simple_metrics.avg_fraction, PRECISION))
            avg_alignment_len_str += '{:<25}'.format(round(simple_metrics.avg_alignment_len, PRECISION))
            avg_mismatches_per_transcript_str += '{:<25}'.format(round(simple_metrics.avg_mismatch_num, PRECISION))

            # na50_str += '{:<25}'.format(simple_metrics.na50)
            # avg_blocks_num_str += '{:<25}'.format(round(simple_metrics.avg_blocks_num, precision))
            # avg_block_len_str += '{:<25}'.format(round(simple_metrics.avg_block_len, precision))
            # avg_mismatches_str += '{:<25}'.format(round(simple_metrics.avg_mismatch_num, precision))

        self.metrics_table.append('\n == ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS == \n')
        self.metrics_table.append(avg_aligned_fraction_str + '\n')
        self.metrics_table.append(avg_alignment_len_str + '\n')
        self.metrics_table.append(avg_mismatches_per_transcript_str + '\n')

        # self.metrics_table.append(na50_str + '\n')
        # self.metrics_table.append(avg_blocks_num_str + '\n')
        # self.metrics_table.append(avg_block_len_str + '\n')
        # self.metrics_table.append(avg_mismatches_str + '\n')


    def add_fusion_misassemble_metrics_to_table(self, transcripts_metrics):
        mis_str_together = '{:<50}'.format('Misassemblies')

        # mis_str_by_blat = '{:<50}'.format('Misassembly candidates by blat')
        # mis_str_by_blast = '{:<50}'.format('Misassembly candidates by blast')
        #if args.fusion_misassemble_analyze:
        #    fusString = '{:<50}'.format('Fusion candidates')

        for i_transcripts in range(len(transcripts_metrics)):
            #if args.fusion_misassemble_analyze:
            #    misString += '{:<25}'.format(transcripts_metrics[i_transcripts].fusion_misassemble_metrics.misassemble_num)
            #    fusString += '{:<25}'.format(transcripts_metrics[i_transcripts].fusion_misassemble_metrics.fusion_num)
            #else:

            mis_str_together += '{:<25}'.format(transcripts_metrics[i_transcripts].simple_metrics.num_misassembled_together)

            # mis_str_by_blat += '{:<25}'.format(transcripts_metrics[i_transcripts].simple_metrics.num_misassembled_by_blat)
            # mis_str_by_blast += '{:<25}'.format(transcripts_metrics[i_transcripts].simple_metrics.num_misassembled_by_blast)

        self.metrics_table.append('\n == ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS == \n')
        self.metrics_table.append(mis_str_together + '\n')

        # self.metrics_table.append(mis_str_by_blat + '\n')
        # self.metrics_table.append(mis_str_by_blast + '\n')
        #if args.fusion_misassemble_analyze:
        #    self.metrics_table.append(fusString + '\n')


    def add_assemble_completeness_metrics_to_table(self, transcripts_metrics, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION):
        database_coverage_str = '{:<50}'.format('Database coverage')

        relative_database_coverage_str = '{:<50}'.format('Relative database coverage')

        assembled_well_genes_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled genes')
        assembled_fully_genes_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled genes')
        covered_well_genes_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered genes')
        covered_fully_genes_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered genes')

        assembled_well_isoforms_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-assembled isoforms')
        assembled_fully_isoforms_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-assembled isoforms')
        covered_well_isoforms_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)) + '%-covered isoforms')
        covered_fully_isoforms_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)) + '%-covered isoforms')


        mean_isoform_cov_str = '{:<50}'.format('Mean isoform coverage')
        mean_isoform_assembly_str = '{:<50}'.format('Mean isoform assembly')

        # cegma_complete_str = '{:<50}'.format('Complete')
        # cegma_partial_str = '{:<50}'.format('Partial')

        busco_complete_str = '{:<50}'.format('Complete')
        busco_partial_str = '{:<50}'.format('Partial')

        geneMarkS_T_genes_str = '{:<50}'.format('Predicted genes')

        for i_transcripts in range(len(transcripts_metrics)):
            isoforms_coverage = transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage
            if isoforms_coverage is not None:
                database_coverage_str += '{:<25}'.format(round(isoforms_coverage.fraction_annotation_mapped, PRECISION))

                relative_database_coverage = isoforms_coverage.relative_database_coverage
                if relative_database_coverage is not None:
                    relative_database_coverage_str += '{:<25}'.format(round(relative_database_coverage.database_coverage, PRECISION))

                assembled_well_genes_str += '{:<25}'.format(isoforms_coverage.num_well_assembled_genes)
                assembled_fully_genes_str += '{:<25}'.format(isoforms_coverage.num_fully_assembled_genes)
                covered_well_genes_str += '{:<25}'.format(isoforms_coverage.num_well_covered_genes)
                covered_fully_genes_str += '{:<25}'.format(isoforms_coverage.num_fully_covered_genes)

                assembled_well_isoforms_str += '{:<25}'.format(isoforms_coverage.num_well_assembled_isoforms)
                assembled_fully_isoforms_str += '{:<25}'.format(isoforms_coverage.num_fully_assembled_isoforms)
                covered_well_isoforms_str += '{:<25}'.format(isoforms_coverage.num_well_covered_isoforms)
                covered_fully_isoforms_str += '{:<25}'.format(isoforms_coverage.num_fully_covered_isoforms)

                mean_isoform_cov_str += '{:<25}'.format(round(isoforms_coverage.avg_covered_fraction, PRECISION))
                mean_isoform_assembly_str += '{:<25}'.format(round(isoforms_coverage.avg_assembled_fraction, PRECISION))

            # if transcripts_metrics[i_transcripts].assembly_completeness_metrics.cegma_metrics is not None:
            #     cegma_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.cegma_metrics
            #     if cegma_metrics != None:
            #         cegma_complete_str += '{:<25}'.format(cegma_metrics.complete_completeness)
            #         cegma_partial_str += '{:<25}'.format(cegma_metrics.partial_completeness)

            busco_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.busco_metrics
            if busco_metrics is not None:
                busco_complete_str += '{:<25}'.format(round(busco_metrics.complete_completeness, PRECISION))
                busco_partial_str += '{:<25}'.format(round(busco_metrics.partial_completeness, PRECISION))
            else:
                busco_complete_str += '{:<25}'.format('0')
                busco_partial_str += '{:<25}'.format('0')

            geneMarkS_T_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.geneMarkS_T_metrics
            if geneMarkS_T_metrics is not None:
                geneMarkS_T_genes_str += '{:<25}'.format(geneMarkS_T_metrics.genes)
            else:
                geneMarkS_T_genes_str += '{:<25}'.format('0')

        self.metrics_table.append('\n == ASSEMBLY COMPLETENESS (SENSITIVITY) == \n')

        isoforms_coverage = transcripts_metrics[0].assembly_completeness_metrics.isoforms_coverage
        if isoforms_coverage is not None:
            self.metrics_table.append(database_coverage_str + '\n')

            relative_database_coverage = isoforms_coverage.relative_database_coverage
            if relative_database_coverage is not None:
                self.metrics_table.append(relative_database_coverage_str + '\n')

            self.metrics_table.append(assembled_well_genes_str + '\n')
            self.metrics_table.append(assembled_fully_genes_str + '\n')
            self.metrics_table.append(covered_well_genes_str + '\n')
            self.metrics_table.append(covered_fully_genes_str + '\n')

            self.metrics_table.append(assembled_well_isoforms_str + '\n')
            self.metrics_table.append(assembled_fully_isoforms_str + '\n')
            self.metrics_table.append(covered_well_isoforms_str + '\n')
            self.metrics_table.append(covered_fully_isoforms_str + '\n')

            self.metrics_table.append(mean_isoform_cov_str + '\n')
            self.metrics_table.append(mean_isoform_assembly_str + '\n')

        # if transcripts_metrics[0].cegma_metrics != None:
        #     self.metrics_table.append('\n == CEGMA METRICS ==\n')
        #     self.metrics_table.append(cegma_complete_str + '\n')
        #     self.metrics_table.append(cegma_partial_str + '\n')

        if busco_complete_str.strip().split()[1:].count('0') != len(transcripts_metrics):
            self.metrics_table.append('\n == BUSCO METRICS == \n')
            self.metrics_table.append(busco_complete_str + '\n')
            self.metrics_table.append(busco_partial_str + '\n')

        if geneMarkS_T_genes_str.strip().split()[1:].count('0') != len(transcripts_metrics):
            self.metrics_table.append('\n == GeneMarkS-T METRICS == \n')
            self.metrics_table.append(geneMarkS_T_genes_str + '\n')


    def add_assemble_correctness_metrics_to_table(self, transcripts_metrics, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION):
        matched_well_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_transcript_threshold * 100)) + '%-matched')
        matched_fully_str = '{:<50}'.format(str(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_transcript_threshold * 100)) + '%-matched')
        unannotated_str = '{:<50}'.format('Unannotated')
        mean_transcript_match_str = '{:<50}'.format('Mean fraction of transcript matched')

        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_coverage = transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage
            if transcripts_coverage is not None:
                matched_well_str += '{:<25}'.format(transcripts_coverage.num_well_covered_transcripts)
                matched_fully_str += '{:<25}'.format(transcripts_coverage.num_fully_covered_transcripts)
                unannotated_str += '{:<25}'.format(transcripts_coverage.num_unannotated_transcripts)
                mean_transcript_match_str += '{:<25}'.format(round(transcripts_coverage.avg_covered_fraction_whole_transcript, PRECISION))

        if transcripts_metrics[0].assembly_correctness_metrics.transcripts_coverage is not None:
            self.metrics_table.append('\n == ASSEMBLY SPECIFICITY == \n')

            self.metrics_table.append(matched_well_str + '\n')
            self.metrics_table.append(matched_fully_str + '\n')
            self.metrics_table.append(unannotated_str + '\n')
            self.metrics_table.append(mean_transcript_match_str + '\n')