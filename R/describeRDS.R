#!/usr/bin/env Rscript
###
library(readr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0) {
  ifile <- args[1]
} else {
  message("Error: Syntax: describeRDS RDSFILE [TRAIN_IDS_LABELLED]")
  quit()
}

writeLines(sprintf("%s:", ifile))
dt <- readRDS(ifile)
setnames(dt, "id1", "protein_id")
writeLines(sprintf("\tRows: %d; Cols: %d", nrow(dt), length(dt)))

writeLines(sprintf("\tColumns, 1:3 and samples:"))
jj <- sort(c(1:3, sample(4:nrow(dt), 12)))
writeLines(sprintf("\t\t%5d.\t%s\n", jj, names(dt)[jj]))

dt_counts <- dt[, .(.N), by=c("subset")]
writeLines(sprintf("\tsubset:%s\tN:%5d", dt_counts$subset, dt_counts$N))
dt_counts <- dt[, .(.N), by=c("subset", "Y")]
writeLines(sprintf("\tsubset:%s\tY:%s\tN:%5d", dt_counts$subset, dt_counts$Y, dt_counts$N))

# Save training set IDs and labels, no data:

dt_train <- dt[subset=="train", .(protein_id, Y)]

#write_delim(dt_train[order(-Y)], "125853_train.tsv", "\t")
#write_delim(dt_train[Y=="neg", .(protein_id)], "125853_train_neg.pid", "\t", col_names=F)
