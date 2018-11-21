#!/usr/bin/env Rscript
###
# Inventory the ML-ready datasets, each specific to either OMIM ID or MP ID.
###
#
library(data.table)

t0 <- proc.time()

ddir <- "data/input"

dfiles <- list.files(ddir, pattern="*\\.rds")
message(sprintf("Dataset files: %d", length(dfiles)))

for (f in dfiles) {
  dt <- readRDS(sprintf("%s/%s", ddir, f))
  setDT(dt)
  writeLines(sprintf("{\"filename\":\"%s\",", f))
  writeLines(sprintf("\t\"rows\": %d,", nrow(dt)))
  writeLines(sprintf("\t\"cols\": %d,", ncol(dt)))
  attrs <- attributes(dt)
  for (key in names(attrs)) {
    if (length(attrs[[key]])==1 & (typeof(attrs[[key]]) %in% c("character","integer","double"))) {
      writeLines(sprintf("\t\"%s\":\"%s\",", key, as.character(attrs[[key]])))
    } else {
      writeLines(sprintf("\t\"%s_type\":\"%s\",", key, length(attrs[[key]])))
      writeLines(sprintf("\t\"%s_length\":%d,", key, length(attrs[[key]])))
    }
  }
  rm(dt)
}

message(sprintf("elapsed time (total): %.2fs",(proc.time()-t0)[3]))

