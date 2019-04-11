import sys
import pandas

from IPython.display import display

kallisto_path = sys.argv[1]
RL = int(sys.argv[2])

df_kallisto = pandas.read_csv(sys.argv[1], delim_whitespace=True)

for cutoff in [2, 3, 4]:
    str_cutoff = 'cutoff ' + str(cutoff)

    df_kallisto[str_cutoff] = df_kallisto.est_counts * RL / df_kallisto.eff_length > cutoff

    print('{}: {}'.format(str_cutoff, sum(df_kallisto[str_cutoff])))

display(df_kallisto.head(10))

