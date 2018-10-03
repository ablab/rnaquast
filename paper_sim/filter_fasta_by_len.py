from Bio import SeqIO

import sys

records = list(SeqIO.parse(sys.argv[1], "fasta"))
new_records = []

for record in records:
    if len(record.seq) >= int(sys.argv[3]):
        new_records.append(record)

SeqIO.write(new_records, sys.argv[2], "fasta")