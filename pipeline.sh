# example of pipeline run

python ../../scripts/short_seqs.py -i Archaea_Eukarya.fasta -o Archaea_Eukarya_names_short.fasta -n -s 200 
mafft --auto Archaea_Eukarya.fasta > Archaea_Eukarya.aln
iqtree -s Archaea_Eukarya.aln -st AA -nt 6 -m LG+R8 -n 0 -pre AE

python ../../scripts/pat_find.py Archaea_Eukarya.aln
Rscript ../../scripts/drawtree.R -a Archaea_Eukarya.aln -t AE.treefile -p 543 

# for subtrees


python ../../scripts/randseqs.py Archaea_Eukarya_names_short.fasta Archaea_Eukarya_names_short_subset.fasta 0.35
mafft Archaea_Eukarya_names_short_subset.fasta > Archaea_Eukarya_names_short_subset.aln
python ../../scripts/pat_find.py Archaea_Eukarya_names_short_subset.aln


iqtree -s Archaea_Eukarya_names_short_subset.aln -st AA -nt 6 -m LG+R8 -n 0 -pre subAE
Rscript ../../scripts/drawtree.R -a Archaea_Eukarya_names_short_subset.aln -t subAE.treefile -p 327 
