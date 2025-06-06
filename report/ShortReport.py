__author__ = 'letovesnoi'

import os

from quast_libs import reporting

# Font of plot captions, axes labels and ticks
font = {'family': 'sans-serif', 'style': 'normal', 'weight': 'medium', 'size': 10}

from general import rqconfig


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
                          'Duplication ratio',
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
             metrics_types[4]: [metrics_labels[12]],
             metrics_types[5]: metrics_labels[13:26],
             metrics_types[6]: metrics_labels[26:28],
             metrics_types[7]: [metrics_labels[28]],
             metrics_types[8]: metrics_labels[29:]}
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

                     metrics_labels[13]: 1, metrics_labels[14]: -1, metrics_labels[15]: 1, metrics_labels[16]: 1,
                     metrics_labels[17]: 1, metrics_labels[18]: 1, metrics_labels[19]: 1, metrics_labels[20]: 1,
                     metrics_labels[21]: 1, metrics_labels[22]: 1, metrics_labels[23]: 1, metrics_labels[24]: 1,
                     metrics_labels[25]: 1,

                     metrics_labels[26]: 1, metrics_labels[27]: 1,

                     metrics_labels[28]: 1,

                     metrics_labels[29]: 2, metrics_labels[30]: 2, metrics_labels[31]: -2, metrics_labels[32]: 1}

        return best_type

    @classmethod
    def get_i_rel_best_metrics(cls, metrics_labels, metrics_dict, best_type):
        i_rel_best_metrics = []
        num_absent = 0
        for i_metric_label in range(len(metrics_labels)):
            metric_label = metrics_labels[i_metric_label]
            if metric_label not in metrics_dict:
                num_absent += 1
            else:
                if metric_label in best_type and (best_type[metric_label] == 2 or best_type[metric_label] == -2):
                    i_rel_best_metrics.append(i_metric_label + 1 - num_absent)
        return i_rel_best_metrics


    @classmethod
    def get_best_values(self, metrics_dict, best_type):
        best_values = {}

        for metric_label in metrics_dict.keys():
            values = metrics_dict[metric_label]
            if metric_label not in best_type:
                best_values[metric_label] = None
            else:
                if best_type[metric_label] == 1 or best_type[metric_label] == 2:
                    max_value = - float('Inf')
                    for i_v in range(len(values)):
                        if values[i_v] == '*':
                            continue
                        float_value = float(values[i_v])
                        if float_value > max_value:
                            max_value = float_value
                            argmax = values[i_v]
                    best_values[metric_label] = argmax
                if best_type[metric_label] == -1 or best_type[metric_label] == -2:
                    min_value = float('Inf')
                    for i_v in range(len(values)):
                        if values[i_v] == '*':
                            continue
                        float_value = float(values[i_v])
                        if float_value < min_value:
                            min_value = float_value
                            argmin = values[i_v]
                    best_values[metric_label] = argmin
                
        return best_values


    def get_table_to_draw(self):
        table_to_draw = []
        table_to_draw.append([self.first_label] + self.metrics_dict[self.first_label])

        for metric_type in self.metrics_type:
            for metric_label in self.metrics_type_labels_dict[metric_type]:
                if metric_label in self.metrics_dict:
                    table_to_draw.append([metric_label] + self.metrics_dict[metric_label])

        return table_to_draw


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
        self.metrics_dict = self.set_metrics_dict(db_genes_metrics, transcripts_metrics, PRECISION, TRANSCRIPT_LENS)

        column_n = len(transcripts_metrics) + 1

        # if there are no transcripts, add one row for annotation metrics:
        # if len(transcripts_metrics) == 0:
        #     column_n = 2

        self.get_report(args, column_n, separated_reports, comparison_report, logger)


    def get_report(self, args, column_n, separated_reports, comparison_report, logger):
        logger.print_timestamp()
        logger.info('Getting SHORT SUMMARY report...')

        if len(separated_reports) != 0:
            self.column_widths = ShortReport.get_column_widths(self.metrics_dict, self.first_label)
        else:
            self.column_widths = [rqconfig.space_label, rqconfig.space_value]

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
        self.print_tsv()

        # TEX:
        self.print_tex(column_n, distribution_report)

        # PDF:
        table_to_draw = self.get_table_to_draw()
        self.print_pdf(args, table_to_draw, separated_reports, comparison_report, logger)

        logger.info('  saved to\n' + '    ' + '{}\n'.format(self.path_txt) + '    ' + '{}\n'.format(self.path_tex) +
                    4 * ' ' + '{}'.format(self.path_pdf))

    @classmethod
    def get_column_widths(cls, metrics_dict, first_label):
        widths = []

        widths.append(rqconfig.space_label)
        for t_label in metrics_dict[first_label]:
            curr_width = len(t_label) + 2
            widths.append(max(curr_width, rqconfig.space_value))

        return widths


    # def get_metrics_dict(self, column_n):
    #     metrics_dict = {}
    #     for line in self.metrics_table:
    #         arr = line.split('  ')
    #         new_arr = []
    #         for substr in arr:
    #             if substr.strip() != '':
    #                 new_arr.append(substr.strip())
    #         if len(new_arr) == column_n:
    #             metrics_dict[new_arr[0]] = new_arr[1:]
    #     return metrics_dict


    def print_txt(self):
        fout_txt_file = open(self.path_txt, 'w')

        fout_txt_file.write('SHORT SUMMARY REPORT \n\n')

        column_width_str = '{:<' + str(self.column_widths[0]) + '}'
        txt_str = column_width_str.format(self.first_label)
        for i_t_label in range(len(self.metrics_dict[self.first_label])):
            t_label = self.metrics_dict[self.first_label][i_t_label]
            column_width_str = '{:<' + str(self.column_widths[i_t_label + 1]) + '}'
            txt_str += column_width_str.format(t_label)
        txt_str += '\n'
        for metric_type in self.metrics_type:
            type_flag = False
            tmp_txt_str_type = '\n == ' + metric_type + ' == \n'
            tmp_txt_str_label = ''
            for i_metric_label in range(len(self.metrics_type_labels_dict[metric_type])):
                metric_label = self.metrics_type_labels_dict[metric_type][i_metric_label]
                if metric_label in self.metrics_dict:
                    type_flag = True

                    column_width_str = '{:<' + str(self.column_widths[0]) + '}'
                    tmp_txt_str_label += column_width_str.format(metric_label)
                    for i_metric_value in range(len(self.metrics_dict[metric_label])):
                        metric_value_str = str(self.metrics_dict[metric_label][i_metric_value])

                        column_width_str = '{:<' + str(self.column_widths[i_metric_value + 1]) + '}'
                        tmp_txt_str_label += column_width_str.format(metric_value_str)
                    tmp_txt_str_label += '\n'
            if type_flag:
                txt_str += tmp_txt_str_type + tmp_txt_str_label.strip() + '\n'
        fout_txt_file.write(txt_str)

        fout_txt_file.close()


    def print_tsv(self):
        fout_tsv_file = open(self.path_tsv, 'w')

        tsv_str = self.first_label
        for i_t_label in range(len(self.metrics_dict[self.first_label])):
            tsv_str += '\t' + self.metrics_dict[self.first_label][i_t_label]
        tsv_str += '\n'

        for metric_type in self.metrics_type:

            for i_metric_label in range(len(self.metrics_type_labels_dict[metric_type])):
                metric_label = self.metrics_type_labels_dict[metric_type][i_metric_label]
                if metric_label in self.metrics_dict:
                    tsv_str += metric_label
                    for i_metric_value in range(len(self.metrics_dict[metric_label])):
                        metric_value = self.metrics_dict[metric_label][i_metric_value]

                        tsv_str += '\t' + str(metric_value)
                    tsv_str += '\n'

        fout_tsv_file.write(tsv_str)


    def print_tex(self, column_n, distribution_report):
        with open(self.path_tex, 'w') as fout_tex_file:
            fout_tex_file.write('\\documentclass[12pt,a4paper]{article}\n')
    
            fout_tex_file.write('\\usepackage{fancyhdr}\n')
            fout_tex_file.write('\\usepackage{graphicx}\n')
            fout_tex_file.write('\\usepackage{placeins}\n')
            fout_tex_file.write('\\usepackage{adjustbox}\n')
            fout_tex_file.write('\n')
    
            fout_tex_file.write('\\begin{document}\n')
    
            fout_tex_file.write('\\pagestyle{fancy}\n')
            fout_tex_file.write('\\fancyhf{}\\n')
            fout_tex_file.write('\\chead{Short summary report}\n')
    
            # TABLE:
            self.add_table_to_tex(fout_tex_file, column_n)
    
            fout_tex_file.write('\\lfoot{generated by rnaQUAST}\n')
    
            # FIGURES:
            if distribution_report is not None:
                short_report_plots = distribution_report.short_report_plots
    
                for plot in short_report_plots:
                    self.add_figure_to_tex(fout_tex_file, plot)
    
            fout_tex_file.write('\\end{document}\n')


    def add_table_to_tex(self, fout_tex_file, column_n):
        fout_tex_file.write('\\begin{table}[t]\n')
        fout_tex_file.write('\\centering\n')

        i_rel_best_metrics = ShortReport.get_i_rel_best_metrics(self.metrics_labels, self.metrics_dict, self.best_type)
        fout_tex_file.write(
            '\\caption {rnaQUAST metrics for assembled transcripts. In each row the best values are indicated with ' \
            '\\textbf{bold}. For the transcript metrics (rows ' + str(i_rel_best_metrics)[1:-1] + \
            ') we highlighted the best \\textbf{relative} values i.e. divided by the total number of transcripts in ' \
            'the corresponding assembly.}')

        fout_tex_file.write(r'\begin{adjustbox}{width=1\textwidth}')
        fout_tex_file.write('\\small')


        fout_tex_file.write('\\begin{tabular}{|l*{' + reporting.val_to_str(column_n) + '}{|r}|}')
        fout_tex_file.write('\\hline')

        column_width_str = '{:<' + str(self.column_widths[0]) + '}'
        tex_str = column_width_str.format(r'\textbf{' + self.first_label + '}')
        for i_t_label in range(len(self.metrics_dict[self.first_label])):
            t_label = self.metrics_dict[self.first_label][i_t_label]
            column_width_str = '{:<' + str(self.column_widths[i_t_label + 1]) + '}'
            tex_str += column_width_str.format(' & ' + r'\textbf{' + t_label.replace('_', r'\_') + '}')
        tex_str += r' \\ \hline\hline' + '\n'
        for metric_type in self.metrics_type:
            type_flag = False
            column_width_str = '{:<' + str(rqconfig.space_type) + '}'
            tmp_tex_str_type = column_width_str.format(r'\multicolumn{' + str(column_n) + r'}{l}{\bf ' + metric_type + '}') + \
                               10 * ' ' + r'\\ \hline' + '\n'
            tmp_tex_str_label = ''
            for i_metric_label in range(len(self.metrics_type_labels_dict[metric_type])):
                metric_label = self.metrics_type_labels_dict[metric_type][i_metric_label]
                if metric_label in self.metrics_dict:
                    type_flag = True

                    column_width_str = '{:<' + str(self.column_widths[0]) + '}'
                    tmp_tex_str_label += column_width_str.format(metric_label.replace('>', '$>$').replace('<', '$<$').replace('%', r'\%'))
                    for i_metric_value in range(len(self.metrics_dict[metric_label])):
                        metric_value_str = str(self.metrics_dict[metric_label][i_metric_value])

                        if metric_value_str == str(self.best_values[metric_label]):
                            metric_value_str = r'\textbf{' + metric_value_str + '}'

                        column_width_str = '{:<' + str(self.column_widths[i_metric_value + 1]) + '}'
                        tmp_tex_str_label += column_width_str.format(' & ' + metric_value_str)
                    tmp_tex_str_label += r' \\' + '\n'
            if type_flag:
                tex_str += tmp_tex_str_type + tmp_tex_str_label.strip() + r' \hline' + '\n'
        fout_tex_file.write(tex_str)

        fout_tex_file.write('\\end{tabular}\n')
        fout_tex_file.write('\\end{adjustbox}\n')
        fout_tex_file.write('\\end{table}\n')

        fout_tex_file.write('\\FloatBarrier\n')
        fout_tex_file.write('\\clearpage\n')


    def add_figure_to_tex(self, fout_tex_file, plot):
        fout_tex_file.write('\\begin{figure}[t]\n')
        fout_tex_file.write('\\centering\n')
        fout_tex_file.write('\\includegraphics[width = \\linewidth]{' + plot.path + '}\n')
        fout_tex_file.write('\\caption{' + plot.caption + '}\n')
        fout_tex_file.write('\\end{figure}\n')
        fout_tex_file.write('\\FloatBarrier\n')
        fout_tex_file.write('\\clearpage\n')
        fout_tex_file.write('\n')


    # full report in PDF format: all tables and plots
    def print_pdf(self, args, table_to_draw, separated_reports, comparison_report, logger):
        pdf_tables_figures = [self.get_pdf_table_figure('Short report', 'generated by rnaQUAST', table_to_draw, self.column_widths, logger)]

        all_pdf_file = None
        if not args.no_plots:
            from quast_libs import plotter  # Do not remove this line! It would lead to a warning in matplotlib.
            try:
                from matplotlib.backends.backend_pdf import PdfPages
                all_pdf_file = PdfPages(self.path_pdf)
            except:
                all_pdf_file = None
        if all_pdf_file:
            # for several files with transcripts select only comparison report pictures:
            if comparison_report is not None:
                short_report_plots = comparison_report.distribution_report.short_report_plots
            else:
                short_report_plots = separated_reports[0].distribution_report.short_report_plots

            self.fill_all_pdf_file(all_pdf_file, pdf_tables_figures, short_report_plots, logger)


    # draw_report_table from quast_libs.plotter:
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
        nrows = len(self.metrics_dict)
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

        colLabels = table_to_draw[0][1:]
        if len(colLabels) == 0:
            colLabels = ['']
        rowLabels = [item[0] for item in table_to_draw[1:]]
        restValues = [item[1:] for item in table_to_draw[1:]]

        matplotlib.pyplot.table(cellText=restValues, rowLabels=rowLabels, colLabels=colLabels,
            colWidths = [float(column_width) / sum(column_widths) for column_width in column_widths[1:]],
            rowLoc = 'left', colLoc='center', cellLoc='right', loc='center')
        return figure


    def fill_all_pdf_file(self, all_pdf, pdf_tables_figures, short_report_plots, logger):
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
        for plot in short_report_plots:
            all_pdf.savefig(plot.fig, bbox_inches='tight')

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
    def set_metrics_dict(self, db_genes_metrics, transcripts_metrics, PRECISION, TRANSCRIPT_LENS):
        self.metrics_dict = {}

        self.metrics_dict[self.first_label] = []
        for i_transcripts in range(len(transcripts_metrics)):
            self.metrics_dict[self.first_label].append(transcripts_metrics[i_transcripts].label)

        if db_genes_metrics is not None:
            # ======= DATABASE METRICS ===========
            self.add_database_metrics_to_table(db_genes_metrics, transcripts_metrics, PRECISION)

        if len(transcripts_metrics) >= 1:
            # ======= BASIC TRANSCRIPTS METRICS ===========
            self.add_basic_metrics_to_table(transcripts_metrics)

            # ======= ALIGNMENT METRICS ===========
            self.add_alignment_metrics_to_table(transcripts_metrics, PRECISION)

            self.add_fusion_misassemble_metrics_to_table(transcripts_metrics)

            # ======= ASSEMBLY COMPLETENESS METRICS ===========
            self.add_assemble_completeness_metrics_to_table(transcripts_metrics, PRECISION)

            # ======= ASSEMBLY CORRECTNESS METRICS ===========
            self.add_assemble_correctness_metrics_to_table(transcripts_metrics, PRECISION)

        return self.metrics_dict


    def add_database_metrics_to_table(self, db_genes_metrics, transcripts_metrics, PRECISION):
        # Database short report metrics:
        num_row = len(transcripts_metrics)
        if len(transcripts_metrics) == 0:
            num_row = 1

        for i_label in range(0, 2):
            self.metrics_dict[self.metrics_labels[i_label]] = []

        for i_transcript in range(num_row):
                self.metrics_dict[self.metrics_labels[0]].append(db_genes_metrics.genes_num)

                self.metrics_dict[self.metrics_labels[1]].append(round(db_genes_metrics.avg_exons_num, PRECISION))

        # self.metrics_table.append(' == DATABASE METRICS == \n')


    def add_basic_metrics_to_table(self, transcripts_metrics):
        # Basic short report metrics:
        for i_label in range(2, 5):
            self.metrics_dict[self.metrics_labels[i_label]] = []

        for i_transcripts in range(len(transcripts_metrics)):
            basic_metrics = transcripts_metrics[i_transcripts].basic_metrics

            if basic_metrics is not None:
                self.metrics_dict[self.metrics_labels[2]].append(basic_metrics.number)
                self.metrics_dict[self.metrics_labels[3]].append(basic_metrics.num_transcripts_500)
                self.metrics_dict[self.metrics_labels[4]].append(basic_metrics.num_transcripts_1000)
            else:
                for i_label in range(2, 5):
                    self.metrics_dict[self.metrics_labels[i_label]].append('*')

        for i_label in range(2, 5):
            if self.metrics_dict[self.metrics_labels[i_label]].count('*') == \
                    len(self.metrics_dict[self.metrics_labels[i_label]]):
                del self.metrics_dict[self.metrics_labels[i_label]]

        # self.metrics_table.append('\n == BASIC TRANSCRIPTS METRICS == \n')


    def add_alignment_metrics_to_table(self, transcripts_metrics, PRECISION):
        # Alignment short report metrics:
        for i_label in range(5, 9):
            self.metrics_dict[self.metrics_labels[i_label]] = []

        for i_transcripts in range(len(transcripts_metrics)):
            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics

            if simple_metrics is not None:
                self.metrics_dict[self.metrics_labels[5]].append(simple_metrics.num_aligned)
                self.metrics_dict[self.metrics_labels[6]].append(simple_metrics.num_unique_aligned)
                self.metrics_dict[self.metrics_labels[7]].append(simple_metrics.num_mul_aligned)
                self.metrics_dict[self.metrics_labels[8]].append(simple_metrics.num_unaligned)
            else:
                for i_label in range(5, 9):
                    self.metrics_dict[self.metrics_labels[i_label]].append('*')

        for i_label in range(5, 9):
            if self.metrics_dict[self.metrics_labels[i_label]].count('*') == \
                    len(self.metrics_dict[self.metrics_labels[i_label]]):
                del self.metrics_dict[self.metrics_labels[i_label]]

        # self.metrics_table.append('\n == ALIGNMENT METRICS == \n')


        # OTHER SECTION:
        for i_label in range(9, 12):
            self.metrics_dict[self.metrics_labels[i_label]] = []

        for i_transcripts in range(len(transcripts_metrics)):
            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics

            if simple_metrics is not None:
                self.metrics_dict[self.metrics_labels[9]].append(round(simple_metrics.avg_fraction, PRECISION))
                self.metrics_dict[self.metrics_labels[10]].append(round(simple_metrics.avg_alignment_len, PRECISION))
                self.metrics_dict[self.metrics_labels[11]].append(round(simple_metrics.avg_mismatch_num, PRECISION))
            else:
                for i_label in range(9, 12):
                    self.metrics_dict[self.metrics_labels[i_label]].append('*')

        for i_label in range(9, 12):
            if self.metrics_dict[self.metrics_labels[i_label]].count('*') == \
                    len(self.metrics_dict[self.metrics_labels[i_label]]):
                del self.metrics_dict[self.metrics_labels[i_label]]

        # self.metrics_table.append('\n == ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS == \n')


    def add_fusion_misassemble_metrics_to_table(self, transcripts_metrics):
        self.metrics_dict[self.metrics_labels[12]] = []

        for i_transcripts in range(len(transcripts_metrics)):
            simple_metrics = transcripts_metrics[i_transcripts].simple_metrics
            if simple_metrics is not None:
                self.metrics_dict[self.metrics_labels[12]].append(simple_metrics.num_misassembled_together)
            else:
                self.metrics_dict[self.metrics_labels[12]].append('*')

        if self.metrics_dict[self.metrics_labels[12]].count('*') == len(self.metrics_dict[self.metrics_labels[12]]):
            del self.metrics_dict[self.metrics_labels[12]]

        # self.metrics_table.append('\n == ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS == \n')


    def add_assemble_completeness_metrics_to_table(self, transcripts_metrics, PRECISION):
        for i_label in range(13, 29):
            self.metrics_dict[self.metrics_labels[i_label]] = []

        for i_transcripts in range(len(transcripts_metrics)):
            isoforms_coverage = transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage
            if isoforms_coverage is not None:
                self.metrics_dict[self.metrics_labels[13]].append(round(isoforms_coverage.fraction_annotation_mapped, PRECISION))

                self.metrics_dict[self.metrics_labels[14]].append(round(isoforms_coverage.avg_duplication_ratio, PRECISION))

                relative_database_coverage = isoforms_coverage.relative_database_coverage
                if relative_database_coverage is not None:
                    self.metrics_dict[self.metrics_labels[15]].append(round(relative_database_coverage.database_coverage, PRECISION))
                else:
                    self.metrics_dict[self.metrics_labels[15]].append('*')

                self.metrics_dict[self.metrics_labels[16]].append(isoforms_coverage.num_well_assembled_genes)
                self.metrics_dict[self.metrics_labels[17]].append(isoforms_coverage.num_fully_assembled_genes)
                self.metrics_dict[self.metrics_labels[18]].append(isoforms_coverage.num_well_covered_genes)
                self.metrics_dict[self.metrics_labels[19]].append(isoforms_coverage.num_fully_covered_genes)

                self.metrics_dict[self.metrics_labels[20]].append(isoforms_coverage.num_well_assembled_isoforms)
                self.metrics_dict[self.metrics_labels[21]].append(isoforms_coverage.num_fully_assembled_isoforms)
                self.metrics_dict[self.metrics_labels[22]].append(isoforms_coverage.num_well_covered_isoforms)
                self.metrics_dict[self.metrics_labels[23]].append(isoforms_coverage.num_fully_covered_isoforms)

                self.metrics_dict[self.metrics_labels[24]].append(round(isoforms_coverage.avg_covered_fraction, PRECISION))
                self.metrics_dict[self.metrics_labels[25]].append(round(isoforms_coverage.avg_assembled_fraction, PRECISION))
            else:
                self.metrics_dict[self.metrics_labels[13]].append('*')
                self.metrics_dict[self.metrics_labels[14]].append('*')
                for i_label in range(16, 26):
                    self.metrics_dict[self.metrics_labels[i_label]].append('*')

            busco_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.busco_metrics
            if busco_metrics is not None:
                self.metrics_dict[self.metrics_labels[26]].append(round(busco_metrics.complete_completeness, PRECISION))
                self.metrics_dict[self.metrics_labels[27]].append(round(busco_metrics.partial_completeness, PRECISION))
            else:
                self.metrics_dict[self.metrics_labels[26]].append('*')
                self.metrics_dict[self.metrics_labels[27]].append('*')

            geneMarkS_T_metrics = transcripts_metrics[i_transcripts].assembly_completeness_metrics.geneMarkS_T_metrics
            if geneMarkS_T_metrics is not None:
                self.metrics_dict[self.metrics_labels[28]].append(geneMarkS_T_metrics.genes)
            else:
                self.metrics_dict[self.metrics_labels[28]].append('*')

        for i_label in range(13, 29):
            if self.metrics_dict[self.metrics_labels[i_label]].count('*') == \
                    len(self.metrics_dict[self.metrics_labels[i_label]]):
                del self.metrics_dict[self.metrics_labels[i_label]]

        # self.metrics_table.append('\n == ASSEMBLY COMPLETENESS (SENSITIVITY) == \n')


    def add_assemble_correctness_metrics_to_table(self, transcripts_metrics, PRECISION):
        for i_label in range(29, 33):
            self.metrics_dict[self.metrics_labels[i_label]] = []

        for i_transcripts in range(len(transcripts_metrics)):
            transcripts_coverage = transcripts_metrics[i_transcripts].assembly_correctness_metrics.transcripts_coverage

            if transcripts_coverage is not None:
                self.metrics_dict[self.metrics_labels[29]].append(transcripts_coverage.num_well_covered_transcripts)
                self.metrics_dict[self.metrics_labels[30]].append(transcripts_coverage.num_fully_covered_transcripts)
                self.metrics_dict[self.metrics_labels[31]].append(transcripts_coverage.num_unannotated_transcripts)
                self.metrics_dict[self.metrics_labels[32]].append(round(transcripts_coverage.avg_covered_fraction_whole_transcript, PRECISION))
            else:
                for i_label in range(29, 33):
                    self.metrics_dict[self.metrics_labels[i_label]].append('*')

        for i_label in range(29, 33):
            if self.metrics_dict[self.metrics_labels[i_label]].count('*') == \
                    len(self.metrics_dict[self.metrics_labels[i_label]]):
                del self.metrics_dict[self.metrics_labels[i_label]]
