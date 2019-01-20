

##########################   Make gene lists from CorShrink  ################################

library(corrplot)
corshrink_mat = get(load("../data/CorShrink2Sparse_all_genes.rda"))
person_tissue_genes = get(load("../data/person_tissue_genes_voom.rda"))
tissue_names = names(person_tissue_genes[1,,1])
gene_names_gtex = as.character(read.table("../data/gene_names_GTEX_V6.txt")[,1])
gene_names_gtex = as.character(sapply(gene_names_gtex, function(z) return(strsplit(z, "[.]")[[1]][1])))
list_genes_qc = read.table("../data/list_genes_qc.txt", header = TRUE)
ensg_names_list_genes = as.character(sapply(as.character(list_genes_qc$ENSG),
                                        function(z) return(strsplit(z, "[.]")[[1]][1])))


tab <- array(0, dim(corshrink_mat)[3])
for(m in 1:dim(corshrink_mat)[3]){
  temp <- corshrink_mat[(8:20), (8:20), m]
  temp1 <- corshrink_mat[-(8:20),-(8:20), m]
  tab[m] <-  (quantile(temp[row(temp) > col(temp)], 0.25)) - (quantile(temp1[row(temp1) > col(temp1)], 0.75))
  cat("We are at gene", m, "\n")
}
imp_ids = order(tab, decreasing = FALSE)[1:1000]
imp_genes = gene_names_gtex[imp_ids]
matched_ids = match(imp_genes, ensg_names_list_genes)
newlist = list_genes_qc[matched_ids[!is.na(matched_ids)],]
table(newlist$chr)
df = newlist[, c("chr", "start.ensembl", "end.ensembl")]
write.table(df, file = "../data/BEDFILES/just_non_brain_corr_genes_ENSEMBL_top_1000.bed",
           quote=FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

col2 <- c("blue", "white", "red")
corrplot(as.matrix(corshrink_mat[,,imp_ids[10]]), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

list_genes_qc = read.table("../data/list_genes_qc.txt", header = TRUE)

match(gene_names_gtex, list_genes_qc$ENSG)
