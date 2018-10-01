import sys
import os

import numpy as np
import pandas
import matplotlib.pyplot as plt

max_degree = 2
postfix = "_transcript"
list_colors=['blue', 'red', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'black']

df = pandas.read_csv(sys.argv[1], delim_whitespace=True)

with open(sys.argv[2], 'r') as fin:
    paths = fin.readlines()
paths = [path.strip() for path in paths]

legend_text = []
for i_path in range(len(paths)):
    path = paths[i_path]
    i_ext = os.path.basename(path).rfind('.')
    label = os.path.basename(path)[:i_ext]
    legend_text.append(label)

    good_transcripts_df = pandas.read_csv(path, header=None)
    good_transcripts = set(good_transcripts_df[0].str[0:-len(postfix)])
    df_filtered = df[df.transcript_id.isin(good_transcripts)]
    bins = np.logspace(0, max_degree, max_degree + 1)
    groups = df_filtered.groupby(np.digitize(df_filtered.TPM, bins))

    x = np.insert(bins, len(bins), 10 ** (max_degree + 1))
    plt.plot(x, np.cumsum(groups.size()), 'ro', color=list_colors[i_path % len(list_colors)])
    # plt.ylim(0, len(df_filtered) + 100)
    plt.xscale('symlog')
    plt.yscale('log')

plt.legend(legend_text, fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))
plt.savefig('plot.png')
plt.show()
