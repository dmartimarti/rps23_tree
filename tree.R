## R script to read alignment, and generate the phylogenetic tree

# libraries

#library(tidyverse)
library(seqinr) # necessary to read alignments
library(ggtree)
library(openxlsx)
# session options
options(width = 220)

# my library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')


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







####################################
### Archaea - Eukarya - Bacteria ###
####################################

# working directory: "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Filipe/rsp23/sub_trees/final_trees"

require(dplyr)
library(tidytree)
library(treeio)

odir <- 'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

# let's create a dataframe with all species within this tree

a.aln = read.alignment('Archaea_n.aln', "fasta", forceToLower = TRUE)
b.aln = read.alignment('Bacteria_n.aln', "fasta", forceToLower = TRUE)
e.aln = read.alignment('Eukarya_n.aln', "fasta", forceToLower = TRUE)
t.aln = read.alignment('tree_of_life.aln', "fasta", forceToLower = TRUE)

sp_names = c(a.aln$nam, b.aln$nam, e.aln$nam)
df = data.frame(sp_names)
colnames(df) = 'IDs'

df = df %>% 
	mutate(Domain = ifelse(IDs %in% a.aln$nam, 'Archaea', 
		ifelse(IDs %in% b.aln$nam, 'Bacteria', 'Eukarya')))


mut = c()
for (i in 1:length(t.aln$nam)){ mut = c(mut, substr(t.aln$seq[[i]][1], 539,539)) }
# creates a dataframe with names and mutation
df2 = data.frame(t.aln$nam, mut)
colnames(df2) = c('IDs', 'mut')

df = df %>% left_join(df2)

df3 = df %>% rename(label = IDs, group = Domain) %>% as_tibble



# load trees
tree = read.tree('tree_files/tree_of_life.aln.treefile')


# split names into two different groups 
groupInfo = split(df$IDs, df$Domain)

# generate a tree
tree1 = groupOTU(tree, groupInfo)

ggtree(tree1,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 0.5, aes(angle = angle)) + geom_hilight(node = 3, fill='darkgreen') +
	theme(legend.position = c(0.5,0.36))

ggsave(file = paste0(odir,'/summary_domain_tree.pdf'), width = 300, height = 300, units = 'mm', scale = 2, device = 'pdf')

# generate a tree
tree2 = groupOTU(tree, groupInfo)
tree2 = reroot(tree2, 3135)
ggtree(tree2, layout = 'circular', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 0.5, aes(angle = angle))  +
	theme(legend.position = c(0.5,0.36))

ggsave(file = paste0(odir,'/summary_root_tree.pdf'), width = 300, height = 300, units = 'mm', scale = 2, device = 'pdf')






#################
## complete tree with colours


new_tree = full_join(as_tibble(tree1), df3, by ="label") %>% 
select(-group.y) %>%
rename(group = group.x) %>%
as.treedata

ggtree(new_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 0.5, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	# scale_color_brewer("group", palette = "Set1", type = "seq") +
	# scale_colour_manual(values = colsel(9, palette = 'sat1')) +
	scale_colour_manual(values = c('red',            # deletion
								   'coral1',         # root
								   'dodgerblue3',    # archaea
								   '#64912D',        # bacteria
								   '#D5D815',        # mutation E
								   'coral1',         # eukarya
								   'gray50',         # mutation K
								   'black',          # mutation R
								   '#EA3B1C')) +
	theme(legend.position = "bottom",
		legend.title = element_blank()) 

ggsave(file = paste0(odir,'/summary_root_tree_colors.pdf'), width = 250, height = 250, units = 'mm', scale = 2, device = 'pdf')











#################
## short final tree 

library(dplyr)
library(seqinr) # necessary to read alignments
library(ggtree)
library(openxlsx)
library(treeio)
library(tidytree)
library(ggplot2)

# session options
options(width = 220)





odir <- 'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

# let's create a dataframe with all species within this tree

a.aln = read.alignment('Archaea_n_short.aln', "fasta", forceToLower = TRUE)
b.aln = read.alignment('Bacteria_n.aln', "fasta", forceToLower = TRUE)
e.aln = read.alignment('Eukarya_n_short.aln', "fasta", forceToLower = TRUE)
t.aln = read.alignment('tree_of_life_short.aln', "fasta", forceToLower = TRUE)

sp_names = c(a.aln$nam, b.aln$nam, e.aln$nam)
df = data.frame(sp_names)
colnames(df) = 'IDs'

df = df %>% 
	mutate(Domain = ifelse(IDs %in% a.aln$nam, 'Archaea', 
		ifelse(IDs %in% b.aln$nam, 'Bacteria', 'Eukarya')))


mut = c()
for (i in 1:length(t.aln$nam)){ mut = c(mut, substr(t.aln$seq[[i]][1], 198,198)) }
# creates a dataframe with names and mutation
df2 = data.frame(t.aln$nam, mut)
colnames(df2) = c('IDs', 'mut')

df = df %>% left_join(df2)

df3 = df %>% rename(label = IDs, group = Domain) %>% as_tibble



# load trees
tree = read.tree('tree_files/tree_of_life_short.aln.treefile')

root_tree = phytools::reroot(as.phylo(tree), node = 1)


# split names into two different groups 
groupInfo = split(df$IDs, df$Domain)

# generate a tree
tree1 = groupOTU(root_tree, groupInfo)



new_tree = full_join(as_tibble(tree1), df3, by ="label") %>% 
select(-group.y) %>%
rename(group = group.x) %>%
as.treedata

ggtree(new_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 0.9, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	# scale_color_brewer("group", palette = "Set1", type = "seq") +
	# scale_colour_manual(values = colsel(9, palette = 'pastel1')) +
	scale_colour_manual(values = c('red',            # deletion
								   'coral1',         # root
								   'dodgerblue3',    # archaea
								   '#64912D',        # bacteria
								   '#D5D815',        # mutation E
								   'coral1',         # eukarya
								   'gray50',         # mutation K
								   'black',          # mutation R
								   '#EA3B1C')) +     
	theme(legend.position = "bottom",
		legend.title = element_blank()) 

ggsave(file = paste0(odir,'/summary_root_tree_colors.pdf'), width = 150, height = 150, units = 'mm', scale = 2, device = 'pdf')



####
# this loooooooong string corresponds to the fixed species names

species = c("Escherichia coli","Pseudomonas putida","Natrialba hulunbeirensis","Haloterrigena daqingensis","Natrinema gari",
"Natrialba asiatica","Natronorubrum tibetense","Natronolimnobius sp","Natronolimnobius aegyptiacus","Natronorubrum texcoconense",
"Haloterrigena turkmenica","Natronolimnobius baerhuensis","Natronococcus amylolyticus","Halovivax ruber",
"Halobacteriales archaeon QH 6 64 20","Halobacteriales archaeon QH 8 64 26","Halococcus hamelinensis",
"Halorientalis sp","Halococcus saccharolyticus","Halobacteriales archaeon SW 9 67 24","Halobacteriales archaeon SW 7 65 23",
"Halobacteriales archaeon SW 8 66 22","Haloferax mucosum","Haloferax lucentense","Haloferax sulfurifontis","Haloferax sp Atlit6N",
"Haloferax sp Atlit48N","Haloferax sp Atlit24N","Haloquadratum walsbyi","Halobacterium jilantaiense","Halorubrum lacusprofundi",
"Halorubrum vacuolatum","Halorubrum sp 48 1 W","Halopenitus malekzadehii","Halorubrum sp C191","Halorubrum tropicale","Halorubrum",
"Halobaculum gomorrense","Halobacteriales archaeon SW 12 71 31","Halarchaeum acidiphilum","Halobacteriales archaeon SW 5 68 122",
"Haloarcula rubripromontorii","Haloarcula amylolytica","Haloarcula argentinensis","Halobacteriales archaeon QH 8 67 27",
"Halobacteriales archaeon QS 1 67 19","Haladaptatus litoreus","Halobacteriales archaeon SW 12 69 24","Halobacteriales archaeon QH 3 68 24",
"Nanohaloarchaea archaeon SW 7 43 1","Nanosalina sp strain J07AB43","Nanohaloarchaea archaeon SG9","Nanohaloarchaea archaeon SW 7 46 7",
"Haloredivivus sp strain G17","archaeon SCGAAA382B04","archaeon 13 1 20CM 2 54 9","candidate divison MSBL1 archaeon SCGCAAA261G05",
"Hadesarchaea archaeon YNP N21","Hadesarchaea archaeon DG331","Candidatus Aenigmarchaeota archaeon 1","Candidatus Aenigmarchaeota archaeon 2",
"archaeon","archaeon Candidatus Huberarchaea CG18","Candidatus Micrarchaeota archaeon","Candidatus Micrarchaeota archaeon","archaeon CG07 land ",
"Candidatus Diapherotrites archaeon","ANME1 cluster archaeon","uncultured archaeon","Methanosarcina barkeri","Methanosarcina sp Ant1",
"Methanosarcina mazei","Methanosarcina sp 2HT1A15","Methanosarcina sp 2HT1A6","Methanosarcina sp WH1","Methanosarcina sp 795",
"Methanosarcina sp DSM","Methanomethylovorans sp PtaU1Bin073","Methanomethylovorans sp PtaU1Bin093","Methanohalophilus portucalensis",
"Methanohalobium evestigatum","Methanosalsum zhilinae","ANME2 cluster archaeon","Candidatus Methanoperedens nitroreducens",
"Candidatus Methanoperedenaceae archaeon","Methanocella arvoryzae","Methanocella paludicola","Methanocella conradii","Methanosphaera sp rholeuAM6",
"Methanosphaera stadtmanae","Methanobacterium sp MZA1","Methanobacterium formicicum","Methanobrevibacter gottschalkii","Methanobrevibacter sp YE315",
"Methanobrevibacter smithii","Methanobrevibacter curvatus","Methanobrevibacter millerae","Methanobacterium lacus","Methanothermobacter tenebrarum","Methanoculleus bourgensis",
"Methanoregula sp PtaU1Bin051","Methanoregula formicica","Euryarchaeota archaeon ADurbBin294","Methanospirillum hungatei","Methanofollis liminatans",
"Euryarchaeota archaeon ADurbBin009","Methanoculleus thermophilus","Euryarchaeota archaeon ADurbBin023","Arc I group archaeon B15fssc0709","Theionarchaea archaeon",
"Euryarchaeota archaeon","Marine Group II euryarchaeote","Marine Group II euryarchaeote","Euryarchaeota archaeon 1","Euryarchaeota archaeon 2","Euryarchaeota archaeon 3",
"Euryarchaeota archaeon 4","Euryarchaeota archaeon 5","Euryarchaeota archaeon 6","Marine Group II euryarchaeote MEDG38","Euryarchaeota archaeon TMED97","Euryarchaeota archaeon 7",
"Euryarchaeota archaeon TMED215","Euryarchaeota archaeon 8","Euryarchaeota archaeon 9","Euryarchaeota archaeon 10","Euryarchaeota archaeon 11",
"uncultured marine group IIIII euryarchaeote KM3 54 D07","uncultured marine group IIIII euryarchaeote SAT1000 18 B12","Euryarchaeota archaeon 12",
"Euryarchaeota archaeon 13","Euryarchaeota archaeon 14","Thermoplasmata archaeon","Marine Group II euryarchaeote","Euryarchaeota archaeon 15","Marine Group II euryarchaeote","Euryarchaeota archaeon 16",
"Euryarchaeota archaeon 17","Marine Group II euryarchaeote","Euryarchaeota archaeon 18","Euryarchaeota archaeon 19","Euryarchaeota archaeon 20","Euryarchaeota archaeon 21",
"Euryarchaeota archaeon 22","Euryarchaeota archaeon 23","Euryarchaeota archaeon 24","Candidatus Methanoplasma termitum","Methanomassiliicoccales archaeon",
"Euryarchaeota archaeon RBG 13 61 15","Euryarchaeota archaeon RBG 16 67 27","Euryarchaeota archaeon 25","Marine Group III euryarchaeote CGEpi1","Euryarchaeota archaeon 26",
"Marine Group III euryarchaeote CGEpi3","Marine Group III euryarchaeote CGEpi2","Euryarchaeota archaeon 27","Euryarchaeota archaeon 28","Thermoplasmatales archaeon SM150",
"Thermoplasmata archaeon M11B2D","Thermoplasmatales archaeon SG8524","Euryarchaeota archaeon  29","Acidiplasma cupricumulans","Acidiplasma aeolicum","Picrophilus torridus",
"Candidatus Altiarchaeales archaeon HGWAltiarchaeales1","Candidatus Altiarchaeum","Candidatus Altiarchaeum","Archaeoglobus veneficus","Archaeoglobus fulgidus","Thermococcus sibiricus",
"Thermococcus sibiricus","Thermococcus sp 2319x1","Thermococcales archaeon 44 46","Thermococcus cleftensis","Thermococcus pacificus","Pyrococcus furiosus","Pyrococcus furiosus",
"Pyrococcus sp ST04","Methanotorris formicicus","Methanocaldococcus villosus","Methanococcus maripaludis OS7","Ignicoccus hospitalis","Metallosphaera yellowstonensis MK1",
"Acidianus brierleyi","Acidianus hospitalis","Sulfolobus acidocaldarius","Sulfolobus sp AB777 L09","Sulfolobus sp AB777 G06","Sulfolobus islandicus","Sulfolobus islandicus",
"Sulfolobus islandicus","uncultured Acidilobus sp OSP8","Acidilobus sp SCGC AC742 E15","archaeon HR02","archaeon HR01","Thaumarchaeota archaeon ex4484 121",
"Candidatus Marsarchaeota G1 archaeon BE D","Candidatus Marsarchaeota G1 archaeon OSP B","Candidatus Marsarchaeota G2 archaeon OSP D",
"Candidatus Marsarchaeota G2 archaeon ECH B 2","Staphylothermus hellenicus","Desulfurococcales archaeon ex4484 58","Thermoproteus tenax",
"Thermoproteus sp CP80","Thermoproteus sp CIS 19","Pyrobaculum aerophilum","Pyrobaculum neutrophilum","Pyrobaculum sp WP30","Pyrobaculum ferrireducens",
"Pyrobaculum oguniense","Thermocladium sp ECH B","Candidatus Heimdallarchaeota archaeon LC 3","Candidatus Thorarchaeota archaeon SMTZ45",
"uncultured marine thaumarchaeote KM3 90 G11","uncultured marine thaumarchaeote KM3 79 D07","uncultured marine crenarchaeote HF4000 APKG5C13",
"Thaumarchaeota archaeon SCGC AC337 F14","uncultured marine crenarchaeote HF4000 APKG3E18","uncultured marine thaumarchaeote KM3 35 E05","uncultured marine thaumarchaeote KM3 144 G01",
"Candidatus Nitrosopumilus sp NM25","Candidatus Nitrosopumilus sediminis","Marine Group I thaumarchaeote SCGC RSA3","Candidatus Nitrosopumilus salaria BD31",
"Candidatus Nitrosopumilus adriaticus","Nitrosopumilus sp","Nitrosopumilus sp","Nitrosarchaeum koreense MY1","Candidatus Nitrosoarchaeum sp",
"Candidatus Nitrosoarchaeum limnia SFB1","Candidatus Nitrosoarchaeum limnia BG20","Candidatus Nitrososphaera sp 13 1 40CM 48 12","Candidatus Nitrocosmicus oleophilus",
"archaeon HR05","Candidatus Bathyarchaeota archaeon CG07","Candidatus Bathyarchaeota archaeon ex4484","Candidatus Bathyarchaeota archaeon B261","Lokiarchaeum sp GC14",
"miscellaneous Crenarchaeota group archaeon SMTZ155","miscellaneous Crenarchaeota group15 archaeon DG45","uncultured korarchaeote","Candidatus Korarchaeota archaeon",
"uncultured korarchaeote","Candidatus Nanoclepta minutus","Candidatus Pacearchaeota archaeon","Archaeon GW2011 AR13","Candidatus Diapherotrites archaeon ADurbBin253",
"Candidatus Pacearchaeota archaeon 1","Candidatus Pacearchaeota archaeon 2","Candidatus Pacearchaeota archaeon 3","Candidatus Pacearchaeota archaeon 4","Candidatus Pacearchaeota archaeon 5",
"Candidatus Pacearchaeota archaeon 6","Candidatus Pacearchaeota archaeon 7","Candidatus Pacearchaeota archaeon 8","Candidatus Woesearchaeota archaeon","archaeon D22",
"Candidatus Woesearchaeota archaeon","Candidatus Woesearchaeota archaeon","Candidatus Woesearchaeota archaeon","Candidatus Pacearchaeota archaeon","Nicotiana attenuata","Juglans regia",
"Citrus clementina","Eucalyptus grandis","Medicago truncatula","Brassica rapa","Arabidopsis thaliana","Triticum aestivum","Oryza sativa","Hordeum vulgare subsp vulgare",
"Hordeum vulgare subsp vulgare","Oryza sativa subsp japonica Rice","Hordeum vulgare subsp vulgare","Zea mays","Panicum hallii var hallii","Cuscuta australis",
"Apostasia shenzhenica","Wollemia nobilis","Glycine max","Lotus japonicus","Glycine max","Cicer arietinum","Populus trichocarpa","Aquilegia coerulea","Nicotiana sylvestris",
"Coffea canephora","Eucalyptus grandis","Anthurium amnicola","Dorcoceras hygrometricum","Gossypium barbadense","Gossypium raimondii","Chlamydomonas reinhardtii","Noccaea caerulescens",
"Auxenochlorella protothecoides","Auxenochlorella protothecoides","Ostreococcus tauri","Plasmodium yoelii","Plasmodium yoelii","Plasmodium sp gorilla clade G2","Plasmodium malariae",
"Plasmodium chabaudi adami","Karenia brevis Red tide dinoflagellate","Daphnia magna","Cryptosporidium ubiquitum","Toxoplasma gondii VAND","Toxoplasma gondii COUG",
"Toxoplasma gondii strain ATCC 50611  Me49","Corethrella appendiculata","Aedes albopictus","Anopheles gambiae","Anopheles christyi","Cuerna arida","Maconellicoccus hirsutus",
"Drosophila yakuba","Glossina palpalis gambiensis","Glossina morsitans morsitans","Glossina morsitans morsitans","Drosophila erecta","Haematobia irritans","Agrilus planipennis",
"Periplaneta americana","Heliconius melpomene","Papilio polytes","Rhodnius neglectus","Daphnia magna","Daphnia magna","Daphnia magna","Daphnia magna","Daphnia magna","Daphnia magna",
"Apis cerana","Habropoda laboriosa","Trachymyrmex zeteki","Trichomalopsis sarcophagae","Orchesella cincta","Pelinobius muticus","Tropilaelaps mercedesae","Dermacentor variabilis",
"Ictidomys tridecemlineatus","Leptonychotes weddellii","Rhinopithecus roxellana","Pongo abelii","Pan troglodytes","Papio anubis","Cebus capucinus imitator","Tarsius syrichta","Macaca fascicularis",
"Lipotes vexillifer","Rhinopithecus bieti","Rhinopithecus bieti","Macaca fascicularis","Sus scrofa","Rhinopithecus bieti","Macaca nemestrina","Macaca nemestrina","Rhinopithecus bieti",
"Rhinopithecus roxellana","Macaca mulatta","Pan troglodytes","Pan paniscus","Pan paniscus","Rhinopithecus bieti","Saimiri boliviensis","Aotus nancymaae","Stegastes partitus","Nothobranchius furzeri",
"Nothobranchius korthausae","Lithobates catesbeiana","Nothobranchius kuhntae","Ictalurus punctatus","Iconisemion striatum","Lonchura striata domestica","Ailuropoda melanoleuca","Limosa lapponica baueri",
"Antrostomus carolinensis","Mus musculus","Micrurus surinamensis","Poecilia latipinna","Caligus rogercresseyi","Stichopus japonicus","Diploscapter pachys","Heligmosomoides polygyrus bakeri","Toxocara canis",
"Trichinella patagoniensis","Trichinella pseudospiralis","Macrostomum lignano","Strongyloides ratti","Rhabditophanes sp KR3021","Biomphalaria glabrata","Arion vulgaris","Pomacea canaliculata","Daphnia magna",
"Capitella teleta","Echinostoma caproni","Schistosoma margrebowiei","Schistosoma mansoni","Absidia repens","Rhizopus microsporus","Basidiobolus meristosporus","Spizellomyces punctatus DAOM BR117",
"Batrachochytrium dendrobatidis JEL423","Batrachochytrium dendrobatidis JEL423","Trichoderma gamsii","Trichoderma arundinaceum","Hypocrea jecorina strain QM6a","Trichoderma harzianum CBS 22695",
"Hypocrea jecorina","Fusarium oxysporum","Gibberella zeae","Fusarium oxysporum","Fusarium proliferatum","Fusarium oxysporum Fusarium vascular wilt","Trichoderma parareesei","Verticillium dahliae",
"Verticillium dahliae","Colletotrichum tofieldiae","Colletotrichum orchidophilum","Colletotrichum chlorophyti","Ustilaginoidea virens","Lomentospora prolificans","Ophiocordyceps polyrhachisfurcata",
"Cordyceps confragosa","Cordyceps fumosorosea","Stachybotrys elegans","Trichoderma longibrachiatum","Madurella mycetomatis","Thielavia terrestris","Chaetomium thermophilum","Sordaria macrospora strain",
"Valsa mali var pyri","Togninia minima","Diaporthe ampelina","Magnaporthe grisea","Aspergillus turcosus","Aspergillus carbonarius","Aspergillus luchuensis","Aspergillus calidoustus","Aspergillus oryzae",
"Aspergillus flavus","Aspergillus niger","Penicillium digitatum","Penicillium coprophilum","Penicillium nordicum","Penicilliopsis zonata","Rasamsonia emersonii","Schizosaccharomyces octosporus","Endocarpon pusillum",
"Fonsecaea multimorphosa","Araucaria cunninghamii","Cladophialophora carrionii","Exophiala mesophila","Talaromyces stipitatus","Talaromyces islandicus","Sphaceloma murrayae","Diplodia corticola",
"Alternaria alternata","Emmonsia crescens","Emmonsia sp CAC2015a","Ajellomyces dermatitidis","Blastomyces gilchristii","Trichophyton verrucosum","Trichophyton rubrum","Trichophyton rubrum",
"Baudoinia panamericana","Zymoseptoria tritici","Zymoseptoria tritici","Pseudogymnoascus sp","Pseudogymnoascus sp","Pseudogymnoascus sp","Pseudogymnoascus destructans","Phialophora cf",
"Erysiphe pulchra","Glarea lozoyensis","Rutstroemia sp NJR2017a BBW","Phialocephala scopiformis","Meliniomyces bicolor E","Cyberlindnera jadinii",
"Lachancea thermotolerans","Metschnikowia bicuspidata","Meyerozyma guilliermondii","Pichia sorbitophila strain","Debaryomyces hansenii strain","Pichia membranifaciens",
"Pichia kudriavzevii","Candida maltosa strain","Saccharomyces cerevisiae","Saccharomyces cerevisiae","Zygosaccharomyces bailii","Zygosaccharomyces parabailii","Zygosaccharomyces rouxii",
"Zygosaccharomyces rouxii","Hanseniaspora guilliermondii","Lachancea lanzarotensis","Wickerhamiella sorbophila","Lipomyces starkeyi","Medicago truncatula","Tuber aestivum","Uncinocarpus reesii",
"Pestalotiopsis fici","Microbotryum intermedium","Puccinia triticina","Melampsora laricipopulina","Lentinula edodes","Pycnoporus cinnabarinus","Laetiporus sulphureus","Kwoniella dejecticola",
"Cryptococcus depauperatus","Kockovaella imperatae","Cryptococcus neoformans","Amanita muscaria","Sistotremastrum niveocremeum","Hypholoma sublateritium","Malassezia globosa","Fopius arisanus",
"Coprinopsis cinerea","Rhizoctonia solani","Pochonia chlamydosporia","Torrubiella hemipterigena","Rhizophagus irregularis","Nematocida sp 1 ERTm6")  





new_tree = full_join(as_tibble(tree1), df3, by ="label") %>% 
select(-group.y) %>%
rename(group = group.x) 

# puting the correct species names
new_tree$label[1:521] = species

new_tree = new_tree %>%
as.treedata



ggtree(new_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 1.1, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	# scale_color_brewer("group", palette = "Set1", type = "seq") +
	# scale_colour_manual(values = colsel(9, palette = 'pastel1')) +
	scale_colour_manual(values = c('red',            # deletion
								   'coral1',         # root
								   'dodgerblue3',    # archaea
								   '#64912D',        # bacteria
								   '#D5D815',        # mutation E
								   'coral1',         # eukarya
								   'gray50',         # mutation K
								   'black',          # mutation R
								   '#EA3B1C')) +     
	theme(legend.position = "bottom",
		legend.title = element_blank()) 

ggsave(file = paste0(odir,'/summary_root_tree_colors_species.pdf'), width = 120, height = 120, units = 'mm', scale = 2, device = 'pdf')






ggtree(new_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 0.9, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	# scale_color_brewer("group", palette = "Set1", type = "seq") +
	# scale_colour_manual(values = colsel(9, palette = 'pastel1')) +
	scale_colour_manual(values = c('red',            # deletion
								   'coral1',         # root
								   'dodgerblue3',    # archaea
								   '#64912D',        # bacteria
								   '#D5D815',        # mutation E
								   'coral1',         # eukarya
								   'gray50',         # mutation K
								   'black',          # mutation R
								   '#EA3B1C')) +     
	geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
	theme(legend.position = "bottom",
		legend.title = element_blank()) 

ggsave(file = paste0(odir,'/summary_root_tree_colors_nodes.pdf'), width = 300, height = 300, units = 'mm', scale = 2, device = 'pdf')









t1 = ggtree(new_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 0.9, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	# scale_color_brewer("group", palette = "Set1", type = "seq") +
	# scale_colour_manual(values = colsel(9, palette = 'pastel1')) +
	scale_colour_manual(values = c('red',            # deletion
								   'coral1',         # root
								   'dodgerblue3',    # archaea
								   '#64912D',        # bacteria
								   '#D5D815',        # mutation E
								   'coral1',         # eukarya
								   'gray50',         # mutation K
								   'black',          # mutation R
								   '#EA3B1C')) +     
	theme(legend.position = "bottom",
		legend.title = element_blank()) 




collapse(t1, node = c(880, 908, 843, 854, 862, 824, 809, 815, 820, 922, 983, 924, 976, 989, 965, 991, 993, 939))











