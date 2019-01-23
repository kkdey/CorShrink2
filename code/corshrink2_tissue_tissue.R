

#########################  Run CorShrink2 on the tissue-tissue correlation matrix  ####################

dat = get(load("../data/person_tissue_genes_voom.rda"))
num_samples_per_tissue = apply(dat[,,1], 2, function(x) return(length(which(!is.na(x)))))
gene_names_gtex = as.character(read.table("../data/gene_names_GTEX_V6.txt")[,1])
gene_names_gtex = as.character(sapply(gene_names_gtex, function(z) return(strsplit(z, "[.]")[[1]][1])))
# devtools::install_github("kkdey/CorShrink")
library(CorShrink)
library(CVXR)
corshrink_mat = array(0, c(dim(dat)[2], dim(dat)[2], 100))
for(n in 1:100){
  out = suppressMessages(suppressWarnings(CorShrink2Data(dat[,,n])))
  corshrink_mat[,,n] = out
  cat("We are at gene", n, "\n")
}
save(corshrink_mat, file = "../data/CorShrink2Sparse_all_genes.rda")


col2 <- c("blue", "white", "red")
corrplot(as.matrix(corshrink_mat[,,95]), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

corrplot(as.matrix(cor(dat[,,95], use = "pairwise.complete.obs")), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")
