__author__ = 'letovesnoi'

import subprocess
import os

import collections

from general.log import get_logger
from general import rqconfig

logger = get_logger(rqconfig.LOGGER_DEFAULT_NAME)

ext_plots = 'png'

# Font of plot captions, axes labels and ticks
font = {'family': 'sans-serif', 'style': 'normal', 'weight': 'medium', 'size': 10}


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


def show_distribution(distribution, step):
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


def add_null_in_distribution(distribution, step, min_key_distr, max_key_distr):
    for i in range(int((max_key_distr - min_key_distr) / step) + 1):
        if round(min_key_distr + i * step, 1) not in distribution:
            distribution[round(min_key_distr + i * step, 1)] = 0
    return distribution


def cumulate(sorted_distribution):
    cur_value = 0
    for key in sorted_distribution:
        cur_value += sorted_distribution[key]
    new_distribution = collections.OrderedDict()
    for key in sorted_distribution:
        new_distribution[key] = cur_value
        cur_value -= sorted_distribution[key]
    return new_distribution


# common routine for Nx-plot and NGx-plot (and probably for others Nyx-plots in the future)
def Nx_plot(list_of_labels, lists_of_lengths, out_dir, pdf_plots_figures, title_str, short_report_visible,
            reference_lengths=None, list_colors=['blue', 'red', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'black']):
    if matplotlib_error:
        return

    if not os.path.exists(out_dir):
        command = 'mkdir {}'.format(out_dir)
        subprocess.call(command, shell=True)

    logger.info('    Drawing ' + title_str.replace('\n', ' ') + ' plot...')

    fig = figure()

    import itertools
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
       pdf_plots_figures.append(fig)

    close()


def get_x_begins_ends_plot(distributions, x_log_scale):
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


def get_y_begins_ends_plot(distributions, y_log_scale):
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


def get_step(def_step, distributions, num_points):
    step = def_step

    if def_step is None:
        x_begin, x_end = get_x_begins_ends_plot(distributions, False)
        step = round((x_end - x_begin + 1) * 1.0 / num_points, 1)

    return step


def plot_compare_histogram(out_dir, title_str, transcripts_distributions, transcripts_metrics, transcripts_labels,
                           label_x, label_y, name_fig, short_report_visible, pdf_plots_figures, y_log_scale=False, def_step=None,
                           list_colors=['blue', 'red', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'black'], num_points=10):
    if matplotlib_error == True:
        return

    if transcripts_distributions == []:
        return

    logger.info('    Drawing cumulative ' + title_str.replace('\n', ' ') + ' histogram...')

    fig = figure()

    title('Cumulative ' + title_str + ' histogram')
    step = get_step(def_step, transcripts_distributions, num_points)
    space = step * 2.0 / 3
    shift = space / len(transcripts_metrics)

    show_transcripts_distributions = {}
    for i_transcripts in range(len(transcripts_distributions)):
        if transcripts_distributions[i_transcripts] == {}:
            continue

        show_transcripts_distributions[i_transcripts] = show_distribution(transcripts_distributions[i_transcripts], step)

    x_begin, x_end = get_x_begins_ends_plot(show_transcripts_distributions, False)

    legend_text = []
    cumulate_ordered_transcripts_distributions = {}
    shift_cumulate_transcripts_distributions = {}
    for i_transcripts in range(len(transcripts_distributions)):
        if transcripts_distributions[i_transcripts] == {}:
            continue

        show_transcripts_distributions[i_transcripts] = add_null_in_distribution(show_transcripts_distributions[i_transcripts], step, x_begin, x_end)

        cumulate_ordered_transcripts_distributions[i_transcripts] = cumulate(collections.OrderedDict(sorted(show_transcripts_distributions[i_transcripts].items())))

        shift_cumulate_transcripts_distributions[i_transcripts] = {}
        for key in show_transcripts_distributions[i_transcripts]:
            shift_cumulate_transcripts_distributions[i_transcripts][key + i_transcripts * shift - space / 2] = cumulate_ordered_transcripts_distributions[i_transcripts][key]

        bar(shift_cumulate_transcripts_distributions[i_transcripts].keys(),
            shift_cumulate_transcripts_distributions[i_transcripts].values(),
            width=shift, color=list_colors[i_transcripts % len(list_colors)])

        legend_text.append(transcripts_labels[i_transcripts])

    y_begin, y_end = get_y_begins_ends_plot(shift_cumulate_transcripts_distributions, y_log_scale)

    if y_log_scale:
        # the interval near 0 will be on a linear scale, so 0 can be displayed
        yscale('symlog', linthreshx=0.1)

    xlim(x_begin - space / 2, x_end + space / 2)
    #ylim(y_begin, y_end)

    legend(legend_text, fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))

    xlabel(label_x)
    ylabel(label_y)

    savefig(os.path.join(out_dir, name_fig + '.' + ext_plots), additional_artists='art', bbox_inches='tight')

    logger.debug('      saved to ' + os.path.join(out_dir, name_fig + '.' + ext_plots))

    if short_report_visible:
       pdf_plots_figures.append(fig)

    close()


def plot_compare_distribution(out_dir, title_str, isoforms_distribution, isoforms_label, transcripts_distributions,
                              transcripts_labels, label_x, label_y, name_fig, short_report_visible, pdf_plots_figures,
                              x_log_scale=False, y_log_scale=False, def_step=None, num_points=100,
                              list_colors=['blue', 'red', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'black']):
    if matplotlib_error == True:
        return

    if transcripts_distributions == [] and (isoforms_distribution is None or isoforms_distribution == {}):
        return

    logger.info('    Drawing cumulative ' + title_str.replace('\n', ' ') + ' plot...')

    fig = figure()

    title('Cumulative ' + title_str + ' plot')

    if isoforms_distribution is not None and len(isoforms_distribution.keys()) > 1:
        step = get_step(def_step, [isoforms_distribution], num_points)
    else:
        step = get_step(def_step, transcripts_distributions, num_points)

    show_transcripts_distributions = []
    for i_transcripts in range(len(transcripts_distributions)):
        show_transcripts_distributions.append(show_distribution(transcripts_distributions[i_transcripts], step))

    show_isoforms_distribution = None
    if isoforms_distribution is not None and len(isoforms_distribution) > 1:
        show_isoforms_distribution = show_distribution(isoforms_distribution, step)

        x_begin, x_end = get_x_begins_ends_plot([show_isoforms_distribution], x_log_scale)

        show_isoforms_distribution = add_null_in_distribution(show_isoforms_distribution, step, x_begin, x_end)

        cumulate_ordered_isoforms_distribution = cumulate(collections.OrderedDict(sorted(show_isoforms_distribution.items())))

        plot(cumulate_ordered_isoforms_distribution.keys(), cumulate_ordered_isoforms_distribution.values(), '--',
             label=isoforms_label, color=list_colors[-1])

        y_begin, y_end = get_y_begins_ends_plot([cumulate_ordered_isoforms_distribution], y_log_scale)
    else:
        x_begin, x_end = get_x_begins_ends_plot(show_transcripts_distributions, x_log_scale)

    cumulate_ordered_transcripts_distribution = []
    for i_transcripts in range(len(transcripts_distributions)):
        if transcripts_distributions[i_transcripts] == {}:
            continue

        show_transcripts_distributions[i_transcripts] = \
            add_null_in_distribution(show_transcripts_distributions[i_transcripts], step, x_begin, x_end)

        cumulate_ordered_transcripts_distribution.append(
            cumulate(collections.OrderedDict(sorted(show_transcripts_distributions[i_transcripts].items()))))

        plot(cumulate_ordered_transcripts_distribution[i_transcripts].keys(),
             cumulate_ordered_transcripts_distribution[i_transcripts].values(), '-',
             label=transcripts_labels[i_transcripts], color=list_colors[i_transcripts % len(list_colors)])

    if show_isoforms_distribution is None:
        y_begin, y_end = get_y_begins_ends_plot(cumulate_ordered_transcripts_distribution, y_log_scale)

    if x_log_scale:
        # the interval near 0 will be on a linear scale, so 0 can be displayed
        xscale('symlog', linthreshx=0.1)
    if y_log_scale:
        # the interval near 0 will be on a linear scale, so 0 can be displayed
        yscale('symlog')

    xlim(x_begin, x_end)
    ylim(y_begin, y_end)

    legend(fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))

    xlabel(label_x)
    ylabel(label_y)

    savefig(os.path.join(out_dir, name_fig + '.' + ext_plots), additional_artists='art', bbox_inches='tight')

    logger.debug('      saved to ' + os.path.join(out_dir, name_fig + '.' + ext_plots))

    if short_report_visible:
       pdf_plots_figures.append(fig)

    close()


def get_label(str, list_labels):
    label = '{}:\n'.format(str)

    for i_label in range(len(list_labels)):
        (key, value) = list_labels[i_label]
        label += '{}={} '.format(key, value)
        if i_label % 3 == 2 and i_label != len(list_labels) - 1:
            label += '\n'

    return label