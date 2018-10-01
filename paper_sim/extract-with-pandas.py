import sys
import os

import numpy as np
import pandas
import matplotlib.pyplot as plt

max_degree = 2
postfix = "_transcript"
list_colors = ['blue', 'red', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'black']

# read simulated data file
df = pandas.read_csv(sys.argv[1], delim_whitespace=True)

def read_paths(txt_path):
    with open(txt_path, 'r') as fin:
        paths = fin.readlines()
    paths = [path.strip() for path in paths]
    return paths

def filter_by_assembled(path_assembled, df):
    # read 50 / 95%-assembled genes / isoforms
    good_transcripts_df = pandas.read_csv(path_assembled, header=None)
    good_transcripts = set(good_transcripts_df[0].str[0:-len(postfix)])
    df_filtered = df[df.transcript_id.isin(good_transcripts)]
    return df_filtered


def plot_coverage_plot(paths, df):
    plt.title('Cumulative plot')
    legend_text = []
    for i_path in range(len(paths)):
        path = paths[i_path]
        i_ext = os.path.basename(path).find('.')
        label = os.path.basename(path)[:i_ext]
        legend_text.append(label)

        bins = np.logspace(0, max_degree, max_degree + 1)

        df_filtered = filter_by_assembled(path, df)

        groups_by_TPM = df_filtered.groupby(np.digitize(df_filtered.TPM, bins))

        x = np.insert(bins, len(bins), 10 ** (max_degree + 1))
        plt.plot(x, np.cumsum(groups_by_TPM.size()), '.-', color=list_colors[i_path % len(list_colors)])
        # plt.ylim(0, len(df_filtered) + 100)
        plt.xscale('symlog')
        plt.yscale('log')

    plt.xlabel('TPM')
    plt.ylabel('Number')
    plt.legend(legend_text, fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))
    plt.savefig('cov_plot.png', additional_artists='art', bbox_inches='tight')
    plt.show()

def get_label(path):
    i_ext = os.path.basename(path).find('.')
    label = os.path.basename(path)[:i_ext]
    return label

def plot_FP(paths_50, paths_95, df):
    zero_cov_num = df[df.TPM == 0].shape[0]
    plt.title('FP histogram (zero coverage {})'.format(zero_cov_num))
    legend_text = []
    x = np.array([0, 1])
    for i_path in range(len(paths_50)):
        path_50 = paths_50[i_path]
        path_95 = paths_95[i_path]

        legend_text.append(get_label(path_50))

        df_filtered_50 = filter_by_assembled(path_50, df)
        df_filtered_95 = filter_by_assembled(path_95, df)
        FP_50 = df_filtered_50[df_filtered_50.TPM == 0].shape[0]
        FP_95 = df_filtered_95[df_filtered_95.TPM == 0].shape[0]
        plt.plot(x, [FP_50, FP_95], '.', color=list_colors[i_path % len(list_colors)])
    plt.xticks(x, ['50%-assembled', '95%-assembled'])
    plt.ylabel('FP')
    plt.legend(legend_text, fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))
    plt.savefig('FP.png', additional_artists='art', bbox_inches='tight')
    plt.show()

plot_coverage_plot(read_paths(sys.argv[2]), df)
plot_coverage_plot(read_paths(sys.argv[3]), df)
plot_FP(read_paths(sys.argv[2]), read_paths(sys.argv[3]), df)