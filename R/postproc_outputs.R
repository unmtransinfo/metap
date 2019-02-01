#!/usr/bin/env Rscript
###
# Postprocess output: aggregate predictions
###
library(readr)
library(data.table)

outdir <- "data/output_post"
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

mods <- dir("data/output", pattern = "(^MP_|^PS|^[0-9]+$)")
mods <- mods[order(mods)]

for (mod in mods) {
  pred_file <- sprintf("%s/%s/blind.pred.tsv", "data/output", mod)
  if (!file.exists(pred_file)) {
    message(sprintf("ERROR: missing %s; skipping.", pred_file))
    next
  }
  pred_this <- read_delim(pred_file, "\t")
  setDT(pred_this)
  setorder(pred_this, -pred.prob, zscore)
  writeLines(sprintf("%s: predictions: %d", mod, nrow(pred_this)))
  ofile_this <- sprintf("%s/blind.pred_%s.tsv", outdir, mod)
  write_delim(pred_this, ofile_this, "\t")
}
