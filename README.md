# Phylogenetic tree for RPS23 protein

This folder contain a few scripts that can be linked as a pipeline to create a phylogenetic tree of InterPro protein sequences. 

Short description of scripts:
 - short_reads.py: given a fasta file, deletes sequences longer than 200 aa, and retreives the species names of each fasta identifier. After this, it launches _mafft_ to run the alignment.
 - pat_find.py: looks for the flanking sequences of the aa of interest for this project, and captures its position. The input is the alignment file from before.
 - drawtree.R: it takes the alignment, the tree (in Newick format) and the position to draw the phylogenetic tree.
 
Python requirements: _Biopython_, _lxml_, _termcolor_
R requirements: _optparse_, _seqinr_, _ggtree_

