#!/usr/bin/env Rscript

library(readxl)
library(data.table)

pheno.list <- read_xlsx("data/output/pheno_list.xlsx", sheet = 1)
setDT(pheno.list)

for(id in pheno.list[, phenotype_id]) {
  fn.list <- list.files(path = paste0("data/output/", id), pattern = "\\d+\\.model")
  cat(id,"has", length(fn.list), "models","\n") 
}