

#########################  Convert gene names to Ensembl names   ###########################

gene_pos = read.table("../data/glist-hg19.nick.txt", header=TRUE)
gene_names = as.character(gene_pos[,4])
library(mygene)
df = queryMany(gene_names, scopes="symbol", fields="ensembl.gene", species="human")

library(ensembldb)
library(EnsDb.Hsapiens.v79)
genes=c("CEACAM8", "CXCL1","CXCL2","CXCL5","CXCL8","MPO","ARG1","FCGR3B",
        "CCL20","CCL2","IL1RN","IL10","TGFB1","IL6","CCL3","CCL4")
genes(EnsDb.Hsapiens.v79, filter=list(GenenameFilter(genes),GeneIdFilter("ENSG", "startsWith")),
      return.type="data.frame", columns=c("gene_id"))

df = genes(EnsDb.Hsapiens.v79, filter=list(GenenameFilter(gene_names),GeneIdFilter("ENSG", "startsWith")),
      return.type="data.frame", columns=c("gene_id"))

new_pos_file = cbind.data.frame(gene_pos[match(df[,2], as.character(gene_pos[,4])),], df[,1])
