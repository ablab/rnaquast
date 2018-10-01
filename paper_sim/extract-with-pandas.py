import sys

import numpy as np
import pandas
import matplotlib.pyplot as plt

max_degree = 3
postfix = "_transcript"
good_transcripts_df = pandas.read_csv(sys.argv[2], header=None)
good_transcripts = set(good_transcripts_df[0].str[0:-len(postfix)])
df = pandas.read_csv(sys.argv[1], delim_whitespace=True)
df_filtered = df[df.transcript_id.isin(good_transcripts)]
bins = np.logspace(0, max_degree, max_degree + 1)
groups = df_filtered.groupby(np.digitize(df_filtered.TPM, bins))

x = np.insert(bins, 0, 0)
plt.plot(x, np.cumsum(groups.size()), 'ro')
# plt.ylim(0, len(df_filtered) + 100)
plt.xscale('symlog')
plt.show()