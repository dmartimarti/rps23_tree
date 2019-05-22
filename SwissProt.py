# parse list with fasta

from Bio import SeqIO
import sys
import numpy as np

def txt_parser(file):
	f = open(file, "r")
	swiss = []
	for line in f:
		line = line.rstrip()
		swiss.append(line)
	return(swiss)

# read the ids from the fasta file
def read_ids(seqs):
	ids = []
	for seq in seqs:
		ids.append(seq.id)
	return ids

# read fasta file
arch = list(SeqIO.parse('Archaea.fasta', "fasta"))
euk = list(SeqIO.parse('Eukarya.fasta', "fasta"))

swiss = np.asarray(txt_parser("swissprot_entry_list.txt"))
euk_id = np.asarray(read_ids(euk))
arch_id = np.asarray(read_ids(arch))

euk_id_hq = np.intersect1d(euk_id, swiss)

arch_id_hq = np.intersect1d(arch_id, swiss)

file = open('euk_HQ.txt','w') 
for i in range(len(euk_id_hq)):
	file.write('>' + euk_id_hq[i] + '\n')

file.close()



file = open('arch_HQ.txt','w') 
for i in range(len(arch_id_hq)):
	file.write('>' + arch_id_hq[i] + '\n')

file.close()



#######

grep -w -A 2 -f  arch_HQ.txt Archaea.fasta | grep -v -- "^--$" > Arch_HQ.fasta 
grep -w -A 2 -f  euk_HQ.txt Eukarya.fasta | grep -v -- "^--$" > Euk_HQ.fasta 
cat Arch_HQ.fasta Euk_HQ.fasta  > All_seqs.fasta

