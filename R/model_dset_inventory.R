#!/usr/bin/env Rscript
###
# Inventory the ML-ready datasets, each specific to either OMIM ID or MP ID.
###
library(data.table)

ddir <- "data/input"

dfiles <- list.files(ddir, pattern="^[0-9]+.rds")
message(sprintf("Dataset files: %d", length(dfiles)))

i <- 0
for (f in dfiles) {
  i <- i + 1
  dt <- readRDS(sprintf("%s/%s", ddir, f))
  writeLines(sprintf("%d. %s:", i, f))
  writeLines(sprintf("\trows: %d", nrow(dt)))
  writeLines(sprintf("\tcols: %d", ncol(dt)))
  attrs <- attributes(dt)
  for (key in c("Date", "mim_id", "omim_title")) {
    writeLines(sprintf("\t%s: %s", key, as.character(attrs[[key]])))
  }
}

