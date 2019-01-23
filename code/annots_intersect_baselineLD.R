
#########################  intersect with Baseline LD annotations   ########################

library(data.table)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

cell_path <- toString(args[1])
folder <- toString(args[2])
baseline_path <- toString(args[3])

#cell_path = "/n/groups/price/kushal/CorShrink2/data/ANNOTATIONS"
#folder="all_corr_genes_ENSEMBL_top_1000"
#baseline_path = "/n/groups/price/kushal/LDSC/baselineLD_v2.1"

for(numchr in 1:22){
  annot = data.frame(fread(paste0("zcat ", cell_path, "/", folder, "/", folder, ".", numchr, ".annot.gz")))
  base = data.frame(fread(paste0("zcat ", baseline_path, "/baselineLD.", numchr, ".annot.gz")))
  annot1 = base$Coding_UCSC*annot[,5]
  annot2 = base$Promoter_UCSC*annot[,5]
  annot3 = base$TFBS_ENCODE*annot[,5]
  sum(annot1)
  sum(annot2)
  sum(annot3)
  df1 = cbind.data.frame(annot[,1:4], annot1)
  df2 = cbind.data.frame(annot[,1:4], annot2)
  df3 = cbind.data.frame(annot[,1:4], annot3)
  colnames(df1) = c("CHR", "BP", "SNP", "CM", "AN")
  colnames(df2) = c("CHR", "BP", "SNP", "CM", "AN")
  colnames(df3) = c("CHR", "BP", "SNP", "CM", "AN")

  if(!dir.exists(paste0(cell_path, "/", folder, "_Coding"))){
    dir.create(paste0(cell_path, "/", folder, "_Coding"))}
  if(!dir.exists(paste0(cell_path, "/", folder, "_Promoter"))){
    dir.create(paste0(cell_path, "/", folder, "_Promoter"))}
  if(!dir.exists(paste0(cell_path, "/", folder, "_TFBS"))){
    dir.create(paste0(cell_path, "/", folder, "_TFBS"))}

  write.table(df1, file = gzfile(paste0(cell_path, "/", folder, "_Coding/", folder, "_Coding.",
                                       numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  write.table(df2, file = gzfile(paste0(cell_path, "/", folder, "_Promoter/", folder, "_Promoter.",
                                        numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  write.table(df3, file = gzfile(paste0(cell_path, "/", folder, "_TFBS/", folder, "_TFBS.",
                                        numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  cat("We are at chr", numchr, "\n")
}





