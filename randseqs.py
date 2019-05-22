# -*- coding: utf-8 -*-
''' 
python code to retrieve the species code from a list 
of InterPro entries

Daniel Martinez, March - 2019
'''

from Bio import SeqIO
import sys
import random

# input and output
inFile = sys.argv[1]
outFile = sys.argv[2]
prop = sys.argv[3]

# read fasta file
records = list(SeqIO.parse(inFile, "fasta"))

n = round(len(records) * float(prop), 0)

print('\nFiltering sequences, now you have a total of', int(n),'\n')

subset = random.sample(records, int(n))

SeqIO.write(subset, outFile, "fasta")