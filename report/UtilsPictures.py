__author__ = 'letovesnoi'

import subprocess
import os
import collections
import itertools

from general.log import get_logger
from general import rqconfig

logger = get_logger(rqconfig.LOGGER_DEFAULT_NAME)

ext_plots = 'png'

# Font of plot captions, axes labels and ticks
font = {'family': 'sans-serif', 'style': 'normal', 'weight': 'medium', 'size': 10}

list_colors=['blue', 'red', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'black']


# checking if matplotlib and pylab is installed:
matplotlib_error = False
try:
    import matplotlib
    matplotlib.use('Agg')  # non-GUI backend
    if matplotlib.__version__.startswith('0'):
        logger.warning('matplotlib version is rather old! Please use matplotlib version 1.0 or higher for better results.')
    from pylab import *
except Exception:
    logger.warning('Can\'t draw plots: please install python-matplotlib and pylab.')
    matplotlib_error = True


class Plot():

    def __init__(self, title_name, label_x, label_y, name_fig, short_report_visible, out_dir,
                 transcripts_labels, transcripts_distribution, isoforms_label=None, isoforms_distribution=None,
                 x_log_scale=False, y_log_scale=False, caption='', def_step=None, ext_plots='png'):
        self.def_step = def_step

        self.path = os.path.join(out_dir, name_fig + '.' + ext_plots)

        self.caption = caption
        self.title_name = title_name
        self.label_x = label_x
        self.label_y = label_y
        self.name_fig = name_fig
        self.x_log_scale = x_log_scale
        self.y_log_scale = y_log_scale

        self.short_report_visible = short_report_visible

        self.transcripts_labels = transcripts_labels
        self.transcripts_distributions = transcripts_distribution

        self.isoforms_label = isoforms_label
        self.isoforms_distribution = isoforms_distribution

        self.num_points = None

        self.title_str = None

        self.fig = None


    # common routine for Nx-plot and NGx-plot (and probably for others Nyx-plots in the future)
    def Nx_plot(self, list_of_labels, lists_of_lengths, out_dir, short_report_plots, title_str, short_report_visible,
                reference_lengths=None):
        if matplotlib_error:
            return

        logger.info('    Drawing ' + title_str.replace('\n', ' ') + ' plot...')

        self.fig = figure()

        for id, (contigs_fpath, lengths) in enumerate(itertools.izip(list_of_labels, lists_of_lengths)):
            if len(lengths) == 0:
                continue
            lengths.sort(reverse=True)
            # calculate values for the plot
            vals_Nx = [0.0]
            vals_l = [lengths[0]]
            lcur = 0
            # if Nx-plot then we just use sum of contigs len, else use reference_length
            lsum = sum(lengths)
            if reference_lengths:
                lsum = reference_lengths[id]
            for l in lengths:
                lcur += l
                x = lcur * 100.0 / lsum
                vals_Nx.append(vals_Nx[-1] + 1e-10) # eps
                vals_l.append(l)
                vals_Nx.append(x)
                vals_l.append(l)
                # add to plot
            vals_Nx.append(vals_Nx[-1] + 1e-10) # eps
            vals_l.append(0.0)

            plot(vals_Nx, vals_l, '-', label=list_of_labels[id], color=list_colors[id % len(list_colors)])

        title(title_str)

        legend(fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))

        xlim(0, 100)
        xlabel('x')
        ylabel('Contig length')

        savefig(os.path.join(out_dir, title_str + '.' + ext_plots), additional_artists='art', bbox_inches='tight')
        logger.debug('      saved to ' + os.path.join(out_dir, title_str + '.' + ext_plots))

        if short_report_visible:
           short_report_plots.append(self)

        close()

        return self.fig


    def plot_compare_distribution(self, short_report_plots, num_points=100):
        self.num_points = num_points

        if matplotlib_error:
            return

        if self.transcripts_distributions == [] and (self.isoforms_distribution is None or self.isoforms_distribution == {}):
            return

        logger.info('    Drawing cumulative ' + self.title_name.replace('\n', ' ') + ' plot...')

        self.fig = figure()

        self.title_str = 'Cumulative ' + self.title_name + ' plot'
        title(self.title_str)

        if self.isoforms_distribution is not None and len(self.isoforms_distribution.keys()) > 1:
            step = Plot.get_step(self.def_step, [self.isoforms_distribution], self.num_points)
        else:
            step = Plot.get_step(self.def_step, self.transcripts_distributions, self.num_points)

        show_transcripts_distributions = []
        for i_transcripts in range(len(self.transcripts_distributions)):
            show_transcripts_distributions.append(Plot.show_distribution(self.transcripts_distributions[i_transcripts], step))

        show_isoforms_distribution = None
        if self.isoforms_distribution is not None and len(self.isoforms_distribution) > 1:
            show_isoforms_distribution = Plot.show_distribution(self.isoforms_distribution, step)

            x_begin, x_end = Plot.get_x_begins_ends_plot([show_isoforms_distribution], self.x_log_scale)

            show_isoforms_distribution = Plot.add_null_in_distribution(show_isoforms_distribution, step, x_begin, x_end)

            cumulate_ordered_isoforms_distribution = Plot.cumulate(collections.OrderedDict(sorted(show_isoforms_distribution.items())))

            plot(cumulate_ordered_isoforms_distribution.keys(), cumulate_ordered_isoforms_distribution.values(), '--',
                 label=self.isoforms_label, color=list_colors[-1])

            y_begin, y_end = Plot.get_y_begins_ends_plot([cumulate_ordered_isoforms_distribution], self.y_log_scale)
        else:
            x_begin, x_end = Plot.get_x_begins_ends_plot(show_transcripts_distributions, self.x_log_scale)

        cumulate_ordered_transcripts_distribution = []
        for i_transcripts in range(len(self.transcripts_distributions)):
            if self.transcripts_distributions[i_transcripts] == {}:
                continue

            show_transcripts_distributions[i_transcripts] = \
                Plot.add_null_in_distribution(show_transcripts_distributions[i_transcripts], step, x_begin, x_end)

            cumulate_ordered_transcripts_distribution.append(
                Plot.cumulate(collections.OrderedDict(sorted(show_transcripts_distributions[i_transcripts].items()))))

            plot(cumulate_ordered_transcripts_distribution[i_transcripts].keys(),
                 cumulate_ordered_transcripts_distribution[i_transcripts].values(), '-',
                 label=self.transcripts_labels[i_transcripts], color=list_colors[i_transcripts % len(list_colors)])

        if show_isoforms_distribution is None:
            y_begin, y_end = Plot.get_y_begins_ends_plot(cumulate_ordered_transcripts_distribution, self.y_log_scale)

        if self.x_log_scale:
            # the interval near 0 will be on a linear scale, so 0 can be displayed
            xscale('symlog', linthreshx=0.1)
        if self.y_log_scale:
            # the interval near 0 will be on a linear scale, so 0 can be displayed
            yscale('symlog')

        xlim(x_begin, x_end)
        ylim(y_begin, y_end)

        legend(fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))

        xlabel(self.label_x)
        ylabel(self.label_y)

        savefig(self.path, additional_artists='art', bbox_inches='tight')

        logger.debug('      saved to ' + self.path)

        if self.short_report_visible:
            short_report_plots.append(self)

        close()

        return self.fig


    def plot_compare_histogram(self, transcripts_metrics, short_report_plots, num_points=10):
        self.num_points = num_points

        if matplotlib_error:
            return

        if self.transcripts_distributions == []:
            return

        logger.info('    Drawing cumulative ' + self.title_name.replace('\n', ' ') + ' histogram...')

        self.fig = figure()

        self.title_str = 'Cumulative ' + self.title_name + ' histogram'
        title(self.title_str)

        step = Plot.get_step(self.def_step, self.transcripts_distributions, self.num_points)
        space = step * 2.0 / 3
        shift = space / len(transcripts_metrics)

        show_transcripts_distributions = {}
        for i_transcripts in range(len(self.transcripts_distributions)):
            if self.transcripts_distributions[i_transcripts] == {}:
                continue

            show_transcripts_distributions[i_transcripts] = Plot.show_distribution(self.transcripts_distributions[i_transcripts], step)

        x_begin, x_end = Plot.get_x_begins_ends_plot(show_transcripts_distributions, False)

        legend_text = []
        cumulate_ordered_transcripts_distributions = {}
        shift_cumulate_transcripts_distributions = {}
        for i_transcripts in range(len(self.transcripts_distributions)):
            if self.transcripts_distributions[i_transcripts] == {}:
                continue

            show_transcripts_distributions[i_transcripts] = Plot.add_null_in_distribution(show_transcripts_distributions[i_transcripts], step, x_begin, x_end)

            cumulate_ordered_transcripts_distributions[i_transcripts] = Plot.cumulate(collections.OrderedDict(sorted(show_transcripts_distributions[i_transcripts].items())))

            shift_cumulate_transcripts_distributions[i_transcripts] = {}
            for key in show_transcripts_distributions[i_transcripts]:
                shift_cumulate_transcripts_distributions[i_transcripts][key + i_transcripts * shift - space / 2] = cumulate_ordered_transcripts_distributions[i_transcripts][key]

            bar(shift_cumulate_transcripts_distributions[i_transcripts].keys(),
                shift_cumulate_transcripts_distributions[i_transcripts].values(),
                width=shift, color=list_colors[i_transcripts % len(list_colors)])

            legend_text.append(self.transcripts_labels[i_transcripts])

        y_begin, y_end = Plot.get_y_begins_ends_plot(shift_cumulate_transcripts_distributions, self.y_log_scale)

        if self.y_log_scale:
            # the interval near 0 will be on a linear scale, so 0 can be displayed
            yscale('symlog', linthreshx=0.1)

        xlim(x_begin - space / 2, x_end + space / 2)
        #ylim(y_begin, y_end)

        legend(legend_text, fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))

        xlabel(self.label_x)
        ylabel(self.label_y)

        savefig(self.path, additional_artists='art', bbox_inches='tight')

        logger.debug('      saved to ' + self.path)

        if self.short_report_visible:
            short_report_plots.append(self)

        close()

        return self.fig


    @classmethod
    def get_x_begins_ends_plot(cls, distributions, x_log_scale):
        x_begin = +inf
        x_end = -inf
        for i_distribution in range(len(distributions)):
            if distributions[i_distribution] == {}:
                continue
            x_begin = min(x_begin, min(distributions[i_distribution].keys()))
            x_end = max(x_end, max(distributions[i_distribution].keys()))

        if x_log_scale:
            x_end = round(x_end * pow(x_end, 0.1))
        else:
            x_end = round(x_end * 1.1)

        return x_begin, x_end


    @classmethod
    def get_y_begins_ends_plot(cls, distributions, y_log_scale):
        y_begin = +inf
        y_end = -inf

        for i_distribution in range(len(distributions)):
            if distributions[i_distribution] == {}:
                continue
            y_begin = min(y_begin, min(distributions[i_distribution].values()))
            y_end = max(y_end, max(distributions[i_distribution].values()))

        if y_log_scale:
            y_end = round(y_end * pow(y_end, 0.1))
        else:
            y_end = round(y_end * 1.1)

        return y_begin, y_end


    @classmethod
    def get_step(cls, def_step, distributions, num_points):
        step = def_step

        if def_step is None:
            x_begin, x_end = cls.get_x_begins_ends_plot(distributions, False)
            step = round((x_end - x_begin + 1) * 1.0 / num_points, 1)

        return step


    @classmethod
    def show_distribution(cls, distribution, step):
        new_distribution = {}
        for key in distribution:
            if key - int(key / step) * step > step / 2:
                new_key = round(int(key / step) * step + step, 1)
            else:
                new_key = round(int(key / step) * step, 1)
            if new_key not in new_distribution:
                new_distribution[new_key] = 0
            new_distribution[new_key] += distribution[key]

        return new_distribution


    @classmethod
    def add_null_in_distribution(cls, distribution, step, min_key_distr, max_key_distr):
        for i in range(int((max_key_distr - min_key_distr) / step) + 1):
            if round(min_key_distr + i * step, 1) not in distribution:
                distribution[round(min_key_distr + i * step, 1)] = 0
        return distribution


    @classmethod
    def cumulate(cls, sorted_distribution):
        cur_value = 0
        for key in sorted_distribution:
            cur_value += sorted_distribution[key]
        new_distribution = collections.OrderedDict()
        for key in sorted_distribution:
            new_distribution[key] = cur_value
            cur_value -= sorted_distribution[key]
        return new_distribution