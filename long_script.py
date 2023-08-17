#!usr/bin/env python

import sys

from Bio import SeqIO

myList = []

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    myList.append([seq_record.id, str(seq_record.seq), len(seq_record)])

myList.sort(key=lambda x: x[2])

print(">", myList[-1][0], sep='')
print(myList[-1][1])
