#!/usr/bin/env Rscript
###
# Inventory the ML-ready datasets, each specific to either OMIM ID or MP ID.
# JSON2TSV: pandas.read_json(intxt, orient='records').to_csv(sep='\t', index=False)
###
#
library(data.table)

t0 <- proc.time()

ddir <- "data/input"

dfiles <- list.files(ddir, pattern="*\\.rds")
message(sprintf("Dataset files: %d", length(dfiles)))

writeLines("[")
for (f in dfiles) {
  dt <- readRDS(sprintf("%s/%s", ddir, f))
  setDT(dt)
  writeLines(sprintf("{\"filename\":\"%s\",", f))
  writeLines(sprintf("\t\"rows\": %d,", nrow(dt)))
  writeLines(sprintf("\t\"cols\": %d,", ncol(dt)))
  attrs <- attributes(dt)
  for (tag in names(attrs)) {
    if (length(attrs[[tag]])==1 & (typeof(attrs[[tag]]) %in% c("character","integer","double"))) {
      writeLines(sprintf("\t\"%s\":\"%s\",", tag, as.character(attrs[[tag]])))
    } else {
      writeLines(sprintf("\t\"%s_type\":\"%s\",", tag, typeof(attrs[[tag]])))
      writeLines(sprintf("\t\"%s_length\":%d,", tag, length(attrs[[tag]])))
    }
  }
  writeLines("}")
  rm(dt)
}
writeLines("]")

message(sprintf("elapsed time (total): %.2fs",(proc.time()-t0)[3]))

