#!/usr/bin/env Rscript

library(data.table)
library(xgboost)

fn <- "data/input/nature_schizophrenia.rds"
#fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)

test.auc <- fread(paste0(fn.base, ".test.auc.tsv"), header = T, sep = "\t", quote = "\"")

dt <- readRDS(fn)
test <- dt[subset == "test"]
train <- dt[subset == "train"]
rm(dt)

X <- data.matrix(test[, -c("Y", "id1", "subset")])
y <- ifelse(test[, Y] == "pos", 1, 0)
dtest <- xgb.DMatrix(X, label = y)
rm(X, y)

test.pred <- test[, .(id1, Y)]
library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)
setnames(test.pred, "id1", "protein_id")

pred.total <- data.table()

for(i in 1:nrow(test.auc)) {
  seed <- test.auc[i, seed]
  auc <- test.auc[i, auc]
  model <- xgb.load(paste0(fn.base, "/", seed, ".model"))
  pred <- predict(model, dtest)
  c.pred <- copy(test.pred)
  c.pred[, pred.prob := pred]
  c.pred[, weight := auc]
  c.pred[, seed := seed]
  pred.total <- rbindlist(list(pred.total, c.pred), use.names = T)
  c.pred <- merge(c.pred, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
  fwrite(c.pred, paste0(fn.base, "/", seed, ".test.pred.tsv"), col.names = T, row.names = F, quote = T, sep = "\t", na = "NA")
}

pred.total.weighted <- pred.total[, .(pred.weighted = weighted.mean(pred.prob, weight)), by = protein_id]
pred.total.weighted <- merge(pred.total.weighted, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)

fwrite(pred.total.weighted[order(-pred.weighted)], file = paste0(fn.base, ".test.ext.pred.weighted.tsv"), sep = "\t", row.names = F, col.names = T, quote = T)

X <- data.matrix(train[, -c("Y", "id1", "subset")])
y <- ifelse(train[, Y] == "pos", 1, 0)
dtrain <- xgb.DMatrix(X, label = y)
rm(X, y)

test.pred <- train[, .(id1, Y)]
setnames(test.pred, "id1", "protein_id")

pred.total <- data.table()

for(i in 1:nrow(test.auc)) {
  seed <- test.auc[i, seed]
  auc <- test.auc[i, auc]
  model <- xgb.load(paste0(fn.base, "/", seed, ".model"))
  pred <- predict(model, dtrain)
  c.pred <- copy(test.pred)
  c.pred[, pred.prob := pred]
  c.pred[, weight := auc]
  c.pred[, seed := seed]
  pred.total <- rbindlist(list(pred.total, c.pred), use.names = T)
  c.pred <- merge(c.pred, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
  fwrite(c.pred, paste0(fn.base, "/", seed, ".train.pred.tsv"), col.names = T, row.names = F, quote = T, sep = "\t", na = "NA")
}

pred.total.weighted <- pred.total[, .(pred.weighted = weighted.mean(pred.prob, weight)), by = protein_id]
pred.total.weighted <- merge(pred.total.weighted, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)

fwrite(pred.total.weighted[order(-pred.weighted)], file = paste0(fn.base, ".train.ext.pred.weighted.tsv"), sep = "\t", row.names = F, col.names = T, quote = T)