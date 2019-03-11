# -*- coding: utf-8 -*-
''' 
python code to retrieve the species code from a list 
of InterPro entries

Daniel Martinez, May - 2019
'''

from Bio import SeqIO
import sys
from lxml import html
import requests
from time import sleep
import re
import subprocess
from termcolor import colored

# input and output
inFile = sys.argv[1]
outFile = sys.argv[2]

# read fasta file
records = list(SeqIO.parse(inFile, "fasta"))

### functions
# filter by length 
def len_filt(seqs):
	empty_list = []
	for seq in seqs:
		if len(seq) < 200:
			empty_list.append(seq)
	return empty_list

# read the ids from the fasta file
def read_ids(seqs):
	ids = []
	for seq in seqs:
		ids.append(seq.id)
	return ids

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

# retrieve species names, and replace spaces for underscores
def retrieve_html(seqs):
	species = []
	l = len(seqs)
	for i in range(len(seqs)):
		name = seqs[i]
		page = requests.get('http://www.ebi.ac.uk/interpro/protein/' + name)
		tree = html.fromstring(page.content)
		name = tree.xpath('//div[@class="prot_gal_desc"]/text()')[2]
		species.append(name.replace(" ", "_"))
		sleep(0.1)
		# Update Progress Bar
		printProgressBar(i + 1, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
	return species

def id_change(seqs, ids):
	for i in range(len(seqs)):
		seqs[i].id = ids[i]
		seqs[i].name = ''
		seqs[i].description = ''

print("Filtering sequences shorter than 200 \n")
short_sequences = len_filt(records)

# number of sequences filtered
filtered = round(((len(records) - len(short_sequences))/len(records)) * 100, 2)

print("Original number of sequences: " + str(len(records)))
print("Final number of sequences: " + str(len(short_sequences)))
print(str(filtered) + "%" + " of sequences filtered \n")

# extracting species and ids
seq_id = read_ids(short_sequences)
print("Getting species names from InterPro website, it may take a while\n")
species = retrieve_html(read_ids(short_sequences))


# filtering bad characters in species
for i in range(len(species)):
	species[i] = re.sub('[/\:(){}<>\'.\[\]-]', '', species[i])

if (len(seq_id) == len(species)) == True:
	pass
else:
	print("vectors are not the same!")

new_id = []
for i in range(len(species)):
	temp = species[i] + "_" + seq_id[i] 
	new_id.append(temp)

id_change(short_sequences, new_id)


# last check
check = len(set(read_ids(short_sequences))) == len(read_ids(short_sequences))
if check == True:
	print("\nProcess finished successfully!\n")
else:
	print("\nYou may have repeated ids, be careful!\n")

SeqIO.write(short_sequences, outFile, "fasta")

mafft_order = 'mafft --auto ' + str(outFile) + ' > ' + outFile.split('.')[0] + '.aln' 


print(colored('Launching mafft!\n', 'green'))
subprocess.run(mafft_order, shell = True)
