#!/usr/bin/env Rscript

library(data.table)
library(xgboost)

fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)
fn.base <- sub("input", "output", fn.base, fixed = T)

auc.test <- fread(paste0(fn.base, "/", "test.auc.tsv"), header = T, sep = "\t", quote = "\"")
setorder(auc.test, -auc)
seed <- auc.test[1, seed]

dt <- readRDS(fn)
dt <- dt[subset == "train"]
set.seed(seed)
train.pos.idx <- sample(dt[Y == "pos", id1], size = 0.8*dt[Y == "pos", .N])
train.neg.idx <- sample(dt[Y == "neg", id1], size = 0.8*dt[Y == "neg", .N])
X <- data.matrix(dt[id1 %in% c(train.pos.idx, train.neg.idx), -c("Y", "id1", "subset")])
y <- ifelse(dt[id1 %in% c(train.pos.idx, train.neg.idx), Y] == "pos", 1, 0)
train.meta <- dt[id1 %in% c(train.pos.idx, train.neg.idx), .(Y, id1, subset)]
dtrain <- xgb.DMatrix(X, label = y)
sumpos <- sum(y == 1)
sumneg <- sum(y == 0)
rm(X, y)

tune.dt <- fread(paste0(fn.base, "/tune.tsv"), header = T, sep = "\t", quote = "\"", na.strings = "NA")
tune.dt <- tune.dt[order(-auc)]
bestTune <- tune.dt[1, ]
p <- list(eta = bestTune$eta, max_depth = bestTune$max_depth, gamma = bestTune$gamma, min_child_weight = bestTune$min_child_weight, subsample = bestTune$subsample, colsample_bytree = bestTune$colsample_bytree, objective = "binary:logistic", scale_pos_weight = sumneg/sumpos)
model <- xgb.cv(params = p, data = dtrain, nrounds = 1000, nfold = 5, early_stopping_rounds = 5, metrics = "auc", prediction = T)
train.meta[, pred.prob := model$pred]
setorder(train.meta, -pred.prob)
fwrite(train.meta, paste0(fn.base, "/", "cv.pred.tsv"), col.names = T, sep = "\t", row.names = F, quote = T, na = "NA")