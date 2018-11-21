#!/usr/bin/env Rscript
###
# Inventory the XGBOOST models generated from OMIM or MP phenotype based datasets.
# For each dataset, multiple models are generated and compared.
###
#
library(xgboost)

t0 <- proc.time()

wd <- getwd()

ddir <- sprintf("%s/data/output", wd)

mdirs <- list.dirs(ddir, recursive=F)
message(sprintf("Model dirs: %d", length(mdirs)))

for (mdir in mdirs) {
  mfiles <- list.files(mdir, pattern="*\\.model$")
  message(sprintf("%s model files: %d", basename(mdir), length(mfiles)))
  for (mfile in mfiles) {
    writeLines(sprintf("%s::%s", basename(mdir), mfile))
    model <- xgb.load(sprintf("%s/%s", mdir, mfile))
    attrs <- xgb.attributes(model)
    for (tag in names(attrs)) {
      if (length(attrs[[tag]])==1 & (typeof(attrs[[tag]]) %in% c("character","integer","double"))) {
        writeLines(sprintf("\t\"%s\":\"%s\",", tag, as.character(attrs[[tag]])))
      } else {
        writeLines(sprintf("\t\"%s_type\":\"%s\",", tag, typeof(attrs[[tag]])))
        writeLines(sprintf("\t\"%s_length\":%d,", tag, length(attrs[[tag]])))
      }
    }
  }
}


message(sprintf("elapsed time (total): %.2fs",(proc.time()-t0)[3]))

