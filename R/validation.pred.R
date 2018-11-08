#!/usr/bin/env Rscript

library(data.table)
library(xgboost)

fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)
fn.base <- sub("input", "output", fn.base, fixed = T)

dt <- readRDS(fn)
dt <- dt[subset == "train"]

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

out.fn <- paste0(fn.base, "/", "test.pred.tsv")

if(file.exists(out.fn)) {
  file.remove(out.fn)
}

auc.test <- fread(paste0(fn.base, "/", "test.auc.tsv"), header = T, sep = "\t", quote = "\"")
for(seed in auc.test[, seed]) {
  set.seed(seed)
  train.pos.idx <- sample(dt[Y == "pos", id1], size = 0.8*dt[Y == "pos", .N])
  train.neg.idx <- sample(dt[Y == "neg", id1], size = 0.8*dt[Y == "neg", .N])
  X <- data.matrix(dt[!id1 %in% c(train.pos.idx, train.neg.idx), -c("Y", "id1", "subset")])
  y <- ifelse(dt[!id1 %in% c(train.pos.idx, train.neg.idx), Y] == "pos", 1, 0)
  test.meta <- dt[!id1 %in% c(train.pos.idx, train.neg.idx), .(Y, id1, subset)]
  dtest <- xgb.DMatrix(X, label = y)
  rm(X, y)
  model <- xgb.load(paste0(fn.base, "/", seed, ".model"))
  pred <- predict(model, dtest)
  test.meta[, pred.prob := pred]
  test.meta[, seed := seed]
  test.meta <- merge(test.meta, protein, by.x = "id1", by.y = "protein_id", all.x = T, sort = F)
  fwrite(test.meta, file = out.fn, sep = "\t", quote = T, na = "NA", col.names = !file.exists(out.fn), append = file.exists(out.fn))
  cat("Dataset with seed", seed, "processed\n")
}
