#!/usr/bin/env Rscript

library(data.table)

download.file("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz", destfile = "data/ncbi/gene_info.gz")

gene.info <- fread("gzcat data/ncbi/gene_info.gz", header = T, sep = "\t", na.strings = "-")
gene.info <- gene.info[`#tax_id` == 9606 & !is.na(chromosome) & type_of_gene == "protein-coding"]
gene.info <- gene.info[, .(GeneID, chromosome, map_location)]
fwrite(gene.info, file = "data/ncbi/gene_info.tsv", col.names = T, sep = "\t", row.names = F, quote = T, na = "None")
if(file.exists("data/ncbi/gene_info.tsv.gz")) {
  file.remove("data/ncbi/gene_info.tsv.gz")
}
system("gzip -9v data/ncbi/gene_info.tsv")