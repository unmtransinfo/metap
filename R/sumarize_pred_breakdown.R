#!/usr/bin/env Rscript

library(data.table)
library(fst)
library(RPostgreSQL)

fn <- "data/input/nature_schizophrenia.rds"

#fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)

nthread <- 82

conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)


test.file.list <- list.files(path = fn.base, pattern = "*.test.pred.breakdown.fst", full.names = T)
train.file.list <- list.files(path = fn.base, pattern = "*.train.pred.breakdown.fst", full.names = T)

all <- data.table()
for(filename in test.file.list) {
  dt <- read_fst(filename, as.data.table = T)
  dt <- dt[feature.long != "intercept"]
  all <- rbindlist(list(all, dt[, .(protein_id, log.odds, feature.long)]))
}

all <- all[, .(median.log.odds = median(log.odds, na.rm = T), mean.log.odds = mean(log.odds, na.rm = T), sd.log.odds = sd(log.odds, na.rm = T), count = .N), by = .(protein_id, feature.long)]
all <- merge(all, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)

fwrite(all, paste0(fn.base, "/", "pred.breakdown.summary.tsv"), col.names = T, row.names = F, sep = "\t", quote = T, na = "NA")

all <- data.table()
for(filename in train.file.list) {
  dt <- read_fst(filename, as.data.table = T)
  dt <- dt[feature.long != "intercept"]
  all <- rbindlist(list(all, dt[, .(protein_id, log.odds, feature.long)]))
}

all <- all[, .(median.log.odds = median(log.odds, na.rm = T), mean.log.odds = mean(log.odds, na.rm = T), sd.log.odds = sd(log.odds, na.rm = T), count = .N), by = .(protein_id, feature.long)]
all <- merge(all, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
fwrite(all, paste0(fn.base, "/", "pred.breakdown.summary.tsv"), col.names = F, row.names = F, sep = "\t", quote = T, na = "NA", append = T)