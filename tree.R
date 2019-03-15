## R script to read alignment, and generate the phylogenetic tree

# libraries

#library(tidyverse)
library(seqinr) # necessary to read alignments
library(ggtree)
library(openxlsx)
# session options
options(width = 220)

# read the alignment
aln = read.alignment('Archaea_Eukarya_short.aln', "fasta", forceToLower = TRUE)
aln2 = read.alignment('alignment_mafft_full.aln', "fasta", forceToLower = TRUE)
aln$nam[7] = 'Caenorhabditis_elegans'
aln2$nam[6] = 'Caenorhabditis_elegans'

# load trees
tree = read.tree('AE.treefile')
tree2 = read.tree('rsp23_newick_tree_full.tre')

### tree 1: the short tree
# simple loop to extract the K or R in position 70
mut = c()
for (i in 1:69){ mut = c(mut, substr(aln$seq[[i]][1], 70,70)) }

# creates a dataframe with names and mutation
df = data.frame(aln$nam, mut)

# split names into two different groups 
groupInfo = split(df$aln.nam, df$mut)

# generate a tree
tree1 = groupOTU(tree, groupInfo)
ggtree(tree1, aes(color = group), layout = 'circular', branch.length = "none") + 
	geom_tiplab(size = 2.5, aes(angle = angle)) + 
	theme(legend.position = c(0.54,0.455))

ggsave(file = 'summary_tree.pdf', width = 120, height = 120, units = 'mm', scale = 2, device = 'pdf')



### tree 2: the large tree
# simple loop to extract the K or R in position 70
mut = c()
for (i in 1:300){ mut = c(mut, substr(aln2$seq[[i]][1], 82,82)) }

# creates a dataframe with names and mutation
df2 = data.frame(aln2$nam, mut)

# split names into two different groups 
groupInfo = split(df2$aln2.nam, df2$mut)

# generate a tree
tree2 = groupOTU(tree2, groupInfo)
ggtree(tree2, aes(color = group), layout = 'circular', branch.length = "none") + 
	geom_tiplab(size = 1.5, aes(angle = angle)) + 
	theme(legend.position = c(0.5,0.4))

ggsave(file = 'full_tree.pdf', width = 120, height = 120, units = 'mm', scale = 2, device = 'pdf')





### analysis of all Opisthokonta organisms (>1500 seqs)
opis_aln = read.alignment('sub_trees/Opis_align.aln', "fasta", forceToLower = TRUE)

mut = c()
for (i in 1:length(opis_aln$nam)){ mut = c(mut, substr(opis_aln$seq[[i]][1], 299,299)) }
# creates a dataframe with names and mutation
df_opis = data.frame(opis_aln$nam, mut)

df_opis %>%
	filter(mut != 'k')


### analysis of all Opisthokonta organisms (>1500 seqs)


plant_aln = read.alignment('sub_trees/Plantae_align.aln', "fasta", forceToLower = TRUE)

mut = c()
for (i in 1:length(plant_aln$nam)){ mut = c(mut, substr(plant_aln$seq[[i]][1], 311,311)) }
# creates a dataframe with names and mutation
df_plant = data.frame(plant_aln$nam, mut)

df_plant %>%
	filter(mut != 'k')

# no hits :(


### analysis of all Opisthokonta organisms (>1500 seqs)


arch_aln = read.alignment('sub_trees/Archaea_align.aln', "fasta", forceToLower = TRUE)

mut = c()
for (i in 1:length(arch_aln$nam)){ mut = c(mut, substr(arch_aln$seq[[i]][1], 116,116)) }
# creates a dataframe with names and mutation
df_arch = data.frame(arch_aln$nam, mut)

df_arch %>%
	filter(mut != 'k') %>% select(arch_aln.nam) %>% write.table('archaea_mutants.txt', quote = F, col.names = F, row.names = F)



##########################################
### analysis of Eukaryotes and Archaea ###
##########################################


# read the alignment
aln = read.alignment('Archaea_Eukarya_names_short.aln', "fasta", forceToLower = TRUE)

# load trees
tree = read.tree('AE.treefile')

mut = c()
for (i in 1:length(aln$nam)){ mut = c(mut, substr(aln$seq[[i]][1], 543,543)) }
# creates a dataframe with names and mutation
df = data.frame(aln$nam, mut)

mutations = df %>% group_by(mut) %>% summarise(Count = n())
deletion = df %>% filter(mut == '-') %>% select(aln.nam)
a_mut = df %>% filter(mut == 'a') %>% select(aln.nam)
e_mut = df %>% filter(mut == 'e') %>% select(aln.nam)
k_mut = df %>% filter(mut == 'k') %>% select(aln.nam)
r_mut = df %>% filter(mut == 'r') %>% select(aln.nam)
t_mut = df %>% filter(mut == 't') %>% select(aln.nam)

l = list('Mutations' = mutations, 'Deletions' = deletion, 'A mutation' = a_mut, 'E mutation' = e_mut, 'K mutation' = k_mut, 'R mutation' = r_mut, 'T mutation' = t_mut)

write.xlsx(l, 'summary.xlsx', colNames = T, rowNames = F)



# split names into two different groups 
groupInfo = split(df$aln.nam, df$mut)

# generate a tree
tree1 = groupOTU(tree, groupInfo)
ggtree(tree1, aes(color = group), layout = 'circular', branch.length = "none") + 
	geom_tiplab(size = 0.5, aes(angle = angle)) + 
	theme(legend.position = c(0.54,0.455))

ggsave(file = 'summary_tree.pdf', width = 120, height = 120, units = 'mm', scale = 2, device = 'pdf')


























