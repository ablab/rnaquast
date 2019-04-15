import sys
import pandas

from IPython.display import display

kallisto_path = sys.argv[1]
RL = int(sys.argv[2])
gff = sys.argv[3]

def get_transcript_to_gene(gff):
    transcript_to_gene = {}
    a_t = 0
    for l in open(gff):
        tokens = l.split()
        if len(tokens) < 3:
            continue
        t = ""
        g = ""
        for i in range(len(tokens)):
            if tokens[i] == "transcript_id":
                t = tokens[i + 1][1:-2]
            if tokens[i] == "gene_id":
                g = tokens[i + 1][1:-2]

        if t and g:
            if t in transcript_to_gene and transcript_to_gene[t] != g:
                a_t += 1
                print("Ambiguous transcript " + t)
            else:
                transcript_to_gene[t] = g
    print("Number of Ambiguous transcript " + str(a_t))
    # print(transcript_to_gene)
    return transcript_to_gene


def get_covered_genes_num(t_to_g, transcript_ids, covered_t):
    covered_genes = set()
    for t_i in range(len(transcript_ids)):
        if covered_t[t_i]:
            covered_genes.add(t_to_g[transcript_ids[t_i]])
    return len(covered_genes)


t_to_g = get_transcript_to_gene(sys.argv[3])

# GET COVERED GENES BY COVERED TRANSCRIPTS WITH KALLISTO
df_kallisto = pandas.read_csv(sys.argv[1], delim_whitespace=True)

df_kallisto['target_id'] = df_kallisto['target_id'].str.replace('_transcript', '')

genes_num = []
transcripts_num = []
for cutoff in [2, 3, 4, 5, 10, 20, 30]:
    str_cutoff = 'cutoff ' + str(cutoff)

    df_kallisto[str_cutoff] = df_kallisto.est_counts * RL * 2 / df_kallisto.length > cutoff

    genes_num.append(get_covered_genes_num(t_to_g, df_kallisto.target_id.values, df_kallisto[str_cutoff].values))
    transcripts_num.append(sum(df_kallisto[str_cutoff]))

print('#transcripts')
print(*transcripts_num, sep='\n')
print('\n#genes')
print(*genes_num, sep='\n')

display(df_kallisto.head(10))


