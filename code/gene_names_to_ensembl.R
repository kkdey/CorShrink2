

#########################  Convert gene names to Ensembl names   ###########################

gene_pos = read.table("../data/glist-hg19.nick.txt", header=TRUE)
gene_names = as.character(gene_pos[,4])
library(mygene)

