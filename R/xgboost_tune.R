#!/usr/bin/env Rscript

library(data.table)
library(xgboost)
library(Matrix)
library(caret)
library(DMwR)


set.seed(1234)

fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)

dt <- readRDS(fn)
dt <- dt[subset == "train"]
dt <- SMOTE(Y ~ ., data = dt)
# meta info before delete
dt.meta <- dt[, .(id1, Y)]
dt.meta[, id1 := floor(id1)]
dt.meta[, Y := ifelse(Y == "pos", 1, 0)]


X <- data.matrix(dt[, -c("Y")])
y <- ifelse(dt[, Y] == "pos", 1, 0)
dtrain <- xgb.DMatrix(X, label = y)
best.auc <- 0.0
best.pred <- data.table(actual = y)
rm(X, y)

tune.grid.xgb <- expand.grid(
  max_depth = c(5, 7, 10),
  eta = c(0.05, 0.1, 0.15, 0.20),
  gamma = c(0.01, 0.1, 1),
  min_child_weight = c(0, 1, 2),
  subsample = seq(0.6, 1, 0.1),
  colsample_bytree = seq(0.5, 1, 0.1))

setDT(tune.grid.xgb)

tune.dt <- data.table()

for(i in 1:nrow(tune.grid.xgb)) {
  model.param <- list(max_depth = tune.grid.xgb[i, max_depth], eta = tune.grid.xgb[i, eta], gamma = tune.grid.xgb[i, gamma], min_child_weight = tune.grid.xgb[i, min_child_weight], subsample = tune.grid.xgb[i, subsample], colsample_bytree = tune.grid.xgb[i, colsample_bytree], objective = "binary:logistic", nthread = 80)
  model <- xgb.cv(params = model.param, data = dtrain, nrounds = 1000, nfold = 5, early_stopping_rounds = 5, metrics = "auc", prediction = T)
  nrounds <- model$best_ntreelimit
  auc <- model$evaluation_log[model$best_iteration, test_auc_mean]
  if(auc > best.auc) {
    best.auc <- auc
    best.pred[, predicted := model$pred]
  }
  record <- data.table(max_depth = tune.grid.xgb[i, max_depth], eta = tune.grid.xgb[i, eta], gamma = tune.grid.xgb[i, gamma], min_child_weight = tune.grid.xgb[i, min_child_weight], subsample = tune.grid.xgb[i, subsample], colsample_bytree = tune.grid.xgb[i, colsample_bytree], nrounds = nrounds, auc = auc)
  tune.dt <- rbindlist(list(tune.dt, record))
}

fwrite(tune.dt[order(-auc)], file = paste0(fn.base, ".xgb.tune.balanced.sample.tsv"), sep = "\t", quote = T, row.names = F, col.names = T)
saveRDS(best.pred, file = paste0(fn.base, ".xgb.auc.balanced.sample.rds"))


dt <- readRDS(fn)
dt <- dt[subset == "train"]

# meta info before delete
dt.meta <- dt[, .(id1, Y)]
dt.meta[, Y := ifelse(Y == "pos", 1, 0)]


X <- data.matrix(dt[, -c("Y")])
y <- ifelse(dt[, Y] == "pos", 1, 0)
dtrain <- xgb.DMatrix(X, label = y)
sumpos <- sum(y == 1)
sumneg <- sum(y == 0)
best.auc <- 0.0
best.pred <- data.table(actual = y)
rm(X, y)

tune.dt <- data.table()

for(i in 1:nrow(tune.grid.xgb)) {
  model.param <- list(max_depth = tune.grid.xgb[i, max_depth], eta = tune.grid.xgb[i, eta], gamma = tune.grid.xgb[i, gamma], min_child_weight = tune.grid.xgb[i, min_child_weight], subsample = tune.grid.xgb[i, subsample], colsample_bytree = tune.grid.xgb[i, colsample_bytree], objective = "binary:logistic", nthread = 80, scale_pos_weight = sumneg/sumpos)
  model <- xgb.cv(params = model.param, data = dtrain, nrounds = 1000, nfold = 5, early_stopping_rounds = 5, metrics = "auc", prediction = T)
  nrounds <- model$best_ntreelimit
  auc <- model$evaluation_log[model$best_iteration, test_auc_mean]
  if(auc > best.auc) {
    best.auc <- auc
    best.pred[, predicted := model$pred]
  }
  record <- data.table(max_depth = tune.grid.xgb[i, max_depth], eta = tune.grid.xgb[i, eta], gamma = tune.grid.xgb[i, gamma], min_child_weight = tune.grid.xgb[i, min_child_weight], subsample = tune.grid.xgb[i, subsample], colsample_bytree = tune.grid.xgb[i, colsample_bytree], nrounds = nrounds, auc = auc)
  tune.dt <- rbindlist(list(tune.dt, record))
}

fwrite(tune.dt[order(-auc)], file = paste0(fn.base, ".xgb.tune.weighted.sample.tsv"), sep = "\t", quote = T, row.names = F, col.names = T)
saveRDS(best.pred, file = paste0(fn.base, ".xgb.auc.weighted.sample.rds"))
