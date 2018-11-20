#!/usr/bin/env Rscript
###
# Inventory the ML-ready datasets, each specific to either OMIM ID or MP ID.
###
#
library(data.table)

t0 <- proc.time()

ddir <- "data/input"

dfiles <- list.files(ddir, pattern="MP.*\\.rds")
writeLines(sprintf("Dataset files: %d", length(dfiles)))

i <- 0
for (f in dfiles) {
  i <- i + 1
  dt <- readRDS(sprintf("%s/%s", ddir, f))
  setDT(dt)
  writeLines(sprintf("%d. %s:", i, f))
  writeLines(sprintf("\trows: %d", nrow(dt)))
  writeLines(sprintf("\tcols: %d", ncol(dt)))
  attrs <- attributes(dt)
  for (key in names(attrs)) {
    if (length(attrs[[key]])==1 & (typeof(attrs[[key]]) %in% c("character","integer","double"))) {
      writeLines(sprintf("\t%s: %s", key, as.character(attrs[[key]])))
    } else {
      writeLines(sprintf("\t%s: length = %d; type = %s", key, length(attrs[[key]]),
	typeof(attrs[[key]])))
    }
  }
  rm(dt)
}

print(sprintf("elapsed time (total): %.2fs",(proc.time()-t0)[3]))

