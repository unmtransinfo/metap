#!/usr/bin/env Rscript

library(data.table)
library(readxl)

pheno.list <- read_xlsx("data/output/pheno_list.xlsx", sheet = 1)
setDT(pheno.list)



for(i in 1:nrow(pheno.list)) {
  dt <- readRDS(paste0("data/input/",pheno.list[i, phenotype_id], ".rds"))
  pheno.list[i, y_pos := dt[Y == "pos" & subset == "train", .N]]
  cat(pheno.list[i, phenotype_id], "processed\n")
}

fwrite(pheno.list[, .(phenotype_id, y_pos)], file = "y_pos_count.tsv", col.names = T, sep = "\t", quote = T, row.names = F, na = "NA")

