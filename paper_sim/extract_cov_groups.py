import sys
import os
import math

# checking if matplotlib and pylab is installed:
matplotlib_error = False
try:
    import matplotlib
    matplotlib.use('Agg')  # non-GUI backend
    if matplotlib.__version__.startswith('0'):
        print('matplotlib version is rather old! Please use matplotlib version 1.0 or higher for better results.')
    from pylab import *
except Exception:
    print('Can\'t draw plots: please install python-matplotlib and pylab.')
    matplotlib_error = True

list_colors=['blue', 'red', 'green', 'yellow', 'magenta', 'orange', 'cyan', 'black']
ext_plots = 'png'

def extract_genes_info(sim_results):
    genes_info = {}
    with open(sim_results, 'r') as fin:
        for line in fin:
            curr_list = line.strip().split()
            if curr_list[0] == 'gene_id':
                continue
            genes_info[curr_list[0]] = \
                {'transcript_id(s)': curr_list[1].split(','),
                 'length': float(curr_list[2]),
                 'effective_length': float(curr_list[3]),
                 'count': float(curr_list[4]),
                 'TPM': float(curr_list[5]),
                 'FPKM': float(curr_list[6])}
    return genes_info


def split_genes_by_counts(sim_results, key='TPM', log_scale=True, cumulative=True, step=1):
    genes_info = extract_genes_info(sim_results)

    genes_bins = {}
    for gene_id in genes_info:
        # print genes_info[gene_id]

        curr_value = genes_info[gene_id][key]
        # print 'TPM ', curr_value

        if curr_value == 0:
            curr_bin = 0
        else:
            curr_key = math.log10(curr_value)
            curr_bin = int(curr_key / step)
            # print 'curr_key ', curr_key

        # print 'curr_bin ', curr_bin
        # print '\n\n'

        if curr_bin not in genes_bins:
            genes_bins[curr_bin] = set()
        genes_bins[curr_bin].add(gene_id)

    return genes_bins


def get_bins_values(assembled_path, genes_bins):
    ass_genes = set()
    bins_values = {}
    with open(assembled_path, 'r') as fin:
        end_len = len('_transcript')
        for line in fin:
            if line[0] == '>':
                gene_id = line.strip().split()[1:-end_len]
                ass_genes.add(gene_id)
    for bin in genes_bins:
        value = len(genes_bins[bin].intersection(ass_genes))
        bins_values[bin] = value
    return bins_values


def get_pathes_list(paths_to_ass):
    assembled_pathes = set()
    with open(paths_to_ass, 'r') as fin:
        for line in fin:
            assembled_pathes.add(line.strip())
    return assembled_pathes


def get_shifted_distribution(distr, i_num, shift, space):
    new_distr = {}
    for key in distr:
        new_distr[key + i_num * shift - space / 2] = distr[key]
    return new_distr


def get_hist(paths_to_ass, sim_results, step=1):
    if matplotlib_error:
        return

    assembled_pathes = get_pathes_list(paths_to_ass)

    space = step ** 10 * 2.0 / 3
    shift = space / len(assembled_pathes)

    figure()

    title('Cumulative histogram')

    genes_bins = split_genes_by_counts(sim_results)

    legend_text = []
    i_assembly = 0

    for path in assembled_pathes:
        ext = '.fa'
        i_ext = os.path.basename(path).find(ext)
        label = os.path.basename(path)[:i_ext]
        legend_text.append(label)

        bins_values = get_bins_values(path, genes_bins)

        bar(bins_values.keys(),
            bins_values.values(),
            width=shift, color=list_colors[i_assembly % len(list_colors)])

        i_assembly += 1

    legend(legend_text, fontsize='x-small', loc='center left', bbox_to_anchor=(1.01, 0.5))

    xlabel('TPM')
    xscale('log')
    ylabel('Number')

    fig_path = 'hist.png'
    savefig(fig_path, additional_artists='art', bbox_inches='tight')

    close()

sim_results = sys.argv[1]
# genes_info = extract_genes_info(sim_results)
# genes_bins = split_genes_by_counts(genes_info)

# print genes_bins

path_to_ass = sys.argv[2]

get_hist(path_to_ass, sim_results)

