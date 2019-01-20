

######################  Single annot column to a proper .annot file  ##############################

library(data.table)

args <- eval(parse(text=commandArgs(T)))
print(args)

cell_path <- as.string(args[1])
output_path <- as.string(args[2])
folder <- as.string(args[3])
bimpath <- as.string(args[4])

#cell_path = "/n/groups/price/kushal/CorShrink2/data/BEDANNOTATIONS"
#output_path = "/n/groups/price/kushal/CorShrink2/data/ANNOTATIONS"
#folder="all_corr_genes_ENSEMBL_top_1000"
#bimpath = "/n/groups/price/kushal/DATA/BIMS"

for(numchr in 1:22){
  bimfile = data.frame(fread(paste0(bimpath, "/", "1000G.EUR.QC.", numchr, ".bim")))
  annot = data.frame(fread(paste0("zcat ", cell_path, "/", folder, "/", folder, ".", numchr, ".annot.gz")))
  df = cbind.data.frame(bimfile[,c(1,4,2,3)], annot)
  write.table(df, file = gzfile(paste0(output_path, "/", folder, "/", folder, ".",
                                       numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)
  colnames(df) = c("CHR", "BP", "SNP", "CM", "AN")
  cat("We are at chr", numchr, "\n")
}
