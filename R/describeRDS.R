#!/usr/bin/env Rscript
###
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0) {
  ifile <- args[1]
} else {
  message("Error: Syntax: describeRDS RDSFILE")
  quit()
}

message(sprintf("Reading %s...", ifile))
dt <- readRDS(ifile)
setnames(dt, "id1", "protein_id")
message(sprintf("Rows: %d; Cols: %d", nrow(dt), length(dt)))

message(sprintf("Columns, 1:3 and samples:"))
jj <- sort(c(1:3, sample(4:nrow(dt), 12)))
message(sprintf("%5d.\t%s\n", jj, names(dt)[jj]))

dt_counts <- dt[, .(.N), by=c("subset", "Y")]
writeLines(sprintf("subset:%s\tY:%s\tN:%5d", dt_counts$subset, dt_counts$Y, dt_counts$N))

