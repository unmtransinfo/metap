#!/usr/bin/env Rscript

library(readxl)
library(data.table)

fn.missing <- "test.auc.tsv"

dt <- read_xlsx("data/output/pheno_list.xlsx", sheet = 1)
setDT(dt)

cat("#!/usr/bin/bash", file = "run1.sh", append = F, sep = "\n")
cat("#!/usr/bin/bash", file = "run2.sh", append = F, sep = "\n")
cat("#!/usr/bin/bash", file = "run3.sh", append = F, sep = "\n")
cat("#!/usr/bin/bash", file = "run4.sh", append = F, sep = "\n")

for(i in 1:nrow(dt)) {
  phenotype_id <- dt[i, phenotype_id]
  if(!file.exists(paste0("data/output/",phenotype_id,"/",fn.missing))) {
    cat(dt[i, run1], file = "run1.sh", append = T, sep = "\n")
    cat(dt[i, run2], file = "run2.sh", append = T, sep = "\n")
    cat(dt[i, run3], file = "run3.sh", append = T, sep = "\n")
    cat(dt[i, run4], file = "run4.sh", append = T, sep = "\n")
  }
}
