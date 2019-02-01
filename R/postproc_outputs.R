#!/usr/bin/env Rscript
###
# Postprocess output: aggregate predictions
###
library(readr)
library(data.table)

N_MAX_HITS <- 500L

#indir <- "/home/oleg/workspace/metap/data/output/" #seaborgium
indir <- "data/output"

outdir <- "data/output_post"
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

mods <- list.dirs(indir, recursive=F, full.names=F)
mods <- mods[order(mods)]

for (mod in mods) {
  pred_file <- sprintf("%s/%s/blind.pred.tsv", indir, mod)
  if (!file.exists(pred_file)) {
    message(sprintf("ERROR: missing %s; skipping.", pred_file))
    next
  }
  pred_this <- read_delim(pred_file, "\t")
  setDT(pred_this)
  setorder(pred_this, -pred.prob, zscore)
  pred_this <- pred_this[1:N_MAX_HITS, ]
  ofile_this <- sprintf("%s/blind.pred_top%d_%s.tsv", outdir, N_MAX_HITS, mod)
  write_delim(pred_this, ofile_this, "\t")
}
