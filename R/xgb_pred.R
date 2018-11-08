#!/usr/bin/env Rscript

library(data.table)
library(xgboost)
library(pROC)
library(caret)

#fn <- commandArgs(T)[1]
fn <- "data/input/nature_schizophrenia.rds"
fn.base <- sub(".rds", "", fn, fixed = T)

#pred.subset <- commandArgs(T)[2]
pred.subset <- "test"

model <- xgb.load(paste0(fn.base, ".xgb.model"))

dt <- readRDS(fn)
X <- data.matrix(dt[subset == pred.subset, -c("Y", "id1", "subset")])
y <- ifelse(dt[subset == pred.subset, Y] == "pos", 1, 0)
dmatrix <- xgb.DMatrix(X, label = y)
sumpos <- sum(y == 1)
sumneg <- sum(y == 0)
rm(X, y)
meta <- dt[subset == pred.subset, .(id1, Y, subset)]

pred <- predict(model, dmatrix)

meta[, pred.prob := pred]
  
rocobj <- roc(meta$Y, meta$pred.prob)
cutoff <- coords(rocobj, "b", ret="t", best.method="closest.topleft")
meta[, Y.pred := ifelse(meta$pred.prob >= cutoff, "pos", "neg")]

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

setnames(meta, "id1", "protein_id")

meta <- merge(meta, protein, by.x = "protein_id", by.y = "protein_id", sort = F)

fwrite(meta, file = paste0(fn.base, ".xgb.pred.", pred.subset, ".tsv"), row.names = F, quote = T, col.names = T, sep = "\t")

