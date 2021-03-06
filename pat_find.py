# -*- coding: utf-8 -*-
''' 
python code to retrieve the species code from a list 
of InterPro entries

Daniel Martinez, March - 2019
'''

from Bio import SeqIO
import sys
from termcolor import colored

# input and output
inFile = sys.argv[1]

# read fasta file
records = list(SeqIO.parse(inFile, "fasta"))

### functions
def pattern_find(seqs):
	temp, pos1, pos2 = 0, [], []
	# find patterns in sequences
	for i in range(len(records)):
		if ('GVEA' in records[i].seq)  == True:
			temp = temp + 1
			pos1.append(records[i].seq.find('GVEA') + 4)
		elif ('GIEA' in records[i].seq) == True:
			temp = temp + 1
			pos2.append(records[i].seq.find('GIEA') + 4)
		else:
			pass
	if temp > 0:
		print('Patterns found in ', str(temp), 'sequences!\n')
	elif temp == 0:
		print('Pattern not found...\n')
	else:
		print('Something went wrong. I dont feel good, Mr. Stark...\n')
	print(colored(('For GVEA pattern, the aa is at position: ' + str(list(set(pos1))[0] + 1)), 'green'))
	print(colored(('For GIEA pattern, the aa is at position: ' + str(list(set(pos2))[0] + 1)+ '\n'), 'blue'))
	pos1 = list(set(pos1))[0] + 1
	return(pos1)

pattern_find(records)

