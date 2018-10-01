import sys
import os

import numpy as np
import pandas
import matplotlib.pyplot as plt

max_degree = 2
postfix_i = "_transcript"
postfix_g = "_gene"
list_colors = ['blue', 'red', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'black']

# read simulated data file
df_g = pandas.read_csv(sys.argv[1], delim_whitespace=True)
df_t = pandas.read_csv(sys.argv[2], delim_whitespace=True)

def read_paths(txt_path):
    with open(txt_path, 'r') as fin:
        paths = fin.readlines()
    paths = [path.strip() for path in paths]
    return paths

paths_g50 = read_paths(sys.argv[3])
paths_g95 = read_paths(sys.argv[4])
paths_i50 = read_paths(sys.argv[5])
paths_i95 = read_paths(sys.argv[6])

def filter_by_assembled(path_assembled, df, postfix):
    # read 50 / 95%-assembled genes / isoforms
    good_transcripts_df = pandas.read_csv(path_assembled, header=None)
    good_transcripts = set(good_transcripts_df[0].str[0:-len(postfix)])
    df_filtered = df[df.transcript_id.isin(good_transcripts)]
    return df_filtered


def plot_coverage_plot(paths, df, postfix, name):
    plt.title('Cumulative plot')
    legend_text = []
    for i_path in range(len(paths)):
        path = paths[i_path]
        i_ext = os.path.basename(path).find('.')
        label = os.path.basename(path)[:i_ext]
        legend_text.append(label)

        bins = np.logspace(0, max_degree, max_degree + 1)

        df_filtered = filter_by_assembled(path, df, postfix)

        groups_by_TPM = df_filtered.groupby(np.digitize(df_filtered.TPM, bins))

        x = np.insert(bins, len(bins), 10 ** (max_degree + 1))
        plt.plot(x, np.cumsum(groups_by_TPM.size()), '.-', color=list_colors[i_path % len(list_colors)])
        # plt.ylim(0, len(df_filtered) + 100)
        plt.xscale('symlog')
        plt.yscale('log')

    plt.xlabel('TPM')
    plt.ylabel('Number')
    plt.legend(legend_text, fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))
    plt.savefig(name, additional_artists='art', bbox_inches='tight')
    plt.show()

def get_label(path):
    i_ext = os.path.basename(path).find('.')
    label = os.path.basename(path)[:i_ext]
    return label

def plot_FP(paths_50, paths_95, df, postfix, name):
    zero_cov_num = df[df.TPM == 0].shape[0]
    plt.title('FP histogram\n (zero coverage genes {})'.format(zero_cov_num))
    legend_text = []
    x = np.array([0, 1])
    for i_path in range(len(paths_50)):
        path_50 = paths_50[i_path]
        path_95 = paths_95[i_path]

        legend_text.append(get_label(path_50))

        df_filtered_50 = filter_by_assembled(path_50, df, postfix)
        df_filtered_95 = filter_by_assembled(path_95, df, postfix)
        FP_50 = df_filtered_50[df_filtered_50.TPM == 0].shape[0]
        FP_95 = df_filtered_95[df_filtered_95.TPM == 0].shape[0]
        plt.plot(x, [FP_50, FP_95], '.', color=list_colors[i_path % len(list_colors)])
    plt.xticks(x, ['50%-assembled', '95%-assembled'])
    plt.ylabel('FP')
    plt.legend(legend_text, fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))
    plt.savefig(name, additional_artists='art', bbox_inches='tight')
    plt.show()

plot_coverage_plot(paths_g50, df_g, postfix_g, '50%-behavior.genes.png')
plot_coverage_plot(paths_g95, df_g, postfix_g, '95%-behavior.genes.png')
plot_coverage_plot(paths_i50, df_t, postfix_i, '50%-behavior.isoforms.png')
plot_coverage_plot(paths_i95, df_t, postfix_i, '95%-behavior.isoforms.png')

plot_FP(paths_g50, paths_g95, df_g, postfix_g, 'FP.genes.png')
plot_FP(paths_i50, paths_g95, df_t, postfix_i, 'FP.isoforms.png')