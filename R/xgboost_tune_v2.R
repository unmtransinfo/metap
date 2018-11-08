#!/usr/bin/env Rscript

library(data.table)
library(xgboost)
library(Matrix)
library(pROC)

fn <- "data/input/nature_schizophrenia.rds"

#fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

dt <- readRDS(fn)
dt <- dt[subset == "train"]

tune.grid.xgb <- expand.grid(
  max_depth = c(5, 7, 10),
  eta = c(0.05, 0.1, 0.15, 0.20),
  gamma = c(0.01, 0.1, 1),
  min_child_weight = c(0, 1, 2),
  subsample = seq(0.8, 1, 0.1),
  colsample_bytree = seq(0.5, 1, 0.1))
setDT(tune.grid.xgb)
test.auc <- 0.0
for(k in 1:10) {
  set.seed(1000 + k)
  train.pos.idx <- sample(dt[Y == "pos", id1], size = 0.8*dt[Y == "pos", .N])
  train.neg.idx <- sample(dt[Y == "neg", id1], size = 0.8*dt[Y == "neg", .N])
  X <- data.matrix(dt[id1 %in% c(train.pos.idx, train.neg.idx), -c("Y", "id1", "subset")])
  y <- ifelse(dt[id1 %in% c(train.pos.idx, train.neg.idx), Y] == "pos", 1, 0)
  dtrain <- xgb.DMatrix(X, label = y)
  sumpos <- sum(y == 1)
  sumneg <- sum(y == 0)
  cv.auc <- 0.0
  cv.pred <- data.table(actual = y)
  rm(X, y)
  X <- data.matrix(dt[!id1 %in% c(train.pos.idx, train.neg.idx), -c("Y", "id1", "subset")])
  y <- ifelse(dt[!id1 %in% c(train.pos.idx, train.neg.idx), Y] == "pos", 1, 0)
  test.pred <- dt[!id1 %in% c(train.pos.idx, train.neg.idx), .(id1, Y)]
  test.pred[, actual := y]
  test.pred <- merge(test.pred, protein, by.x = "id1", by.y = "protein_id", all.x = T, sort = F)
  dtest <- xgb.DMatrix(X, label = y)
  rm(X, y)
  tune.dt <- data.table()
  
  for(i in 1:nrow(tune.grid.xgb)) {
    model.param <- list(max_depth = tune.grid.xgb[i, max_depth], eta = tune.grid.xgb[i, eta], gamma = tune.grid.xgb[i, gamma], min_child_weight = tune.grid.xgb[i, min_child_weight], subsample = tune.grid.xgb[i, subsample], colsample_bytree = tune.grid.xgb[i, colsample_bytree], objective = "binary:logistic", nthread = 80, scale_pos_weight = sumneg/sumpos)
    model <- xgb.cv(params = model.param, data = dtrain, nrounds = 1000, nfold = 5, early_stopping_rounds = 5, metrics = "auc", prediction = T)
    nrounds <- model$best_ntreelimit
    auc <- model$evaluation_log[model$best_iteration, test_auc_mean]
    record <- data.table(seed = (1000+k), max_depth = tune.grid.xgb[i, max_depth], eta = tune.grid.xgb[i, eta], gamma = tune.grid.xgb[i, gamma], min_child_weight = tune.grid.xgb[i, min_child_weight], subsample = tune.grid.xgb[i, subsample], colsample_bytree = tune.grid.xgb[i, colsample_bytree], nrounds = nrounds, auc = auc)
    tune.dt <- rbindlist(list(tune.dt, record))
    if(auc > cv.auc) {
      cv.auc <- auc
      cv.pred[, pred.prob := model$pred]
    }
  }
  tune.dt <- tune.dt[order(-auc)]
  bestTune <- tune.dt[1, ]
  p <- list(eta = bestTune$eta, max_depth = bestTune$max_depth, gamma = bestTune$gamma, min_child_weight = bestTune$min_child_weight, subsample = bestTune$subsample, colsample_bytree = bestTune$colsample_bytree, objective = "binary:logistic", scale_pos_weight = sumneg/sumpos)
  model <- xgb.train(params = p, data = dtrain, nrounds = bestTune$nrounds, verbose = 1, nthread = 80)
  pred <- predict(model, dtest)
  test.pred[, pred.prob := pred]
  auc <- auc(test.pred$actual, test.pred$pred.prob)
  cat("best CV AUC = ", cv.auc, ", current test AUC = ", auc, ", best test AUC = ", test.auc, "\n")
  if(auc > test.auc) {
    test.auc <- auc
    xgb.save(model, paste0(fn.base, ".xgb.weighted.model"))
    fwrite(cv.pred[order(-pred.prob)], paste0(fn.base, ".weighted.cv.pred.tsv"), col.names = T, row.names = F, sep = "\t", quote = T, na = "NA")
    fwrite(test.pred[order(-pred.prob)], paste0(fn.base, ".weighted.test.pred.tsv"), col.names = T, row.names = F, sep = "\t", quote = T, na = "NA")
    fwrite(tune.dt, paste0(fn.base, ".xgb.tune.weighted.tsv"), col.names = T, row.names = F, sep = "\t", quote = T, na = "NA")
  }
}
