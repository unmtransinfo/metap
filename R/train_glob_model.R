#!/usr/bin/env Rscript

library(data.table)
library(xgboost)

nthread <- 4

fn <- "data/input/nature_schizophrenia.rds"
#fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)

dt <- readRDS(fn)
dt <- dt[subset == "train"]

X <- data.matrix(dt[, -c("Y", "id1", "subset")])
y <- ifelse(dt[, Y] == "pos", 1, 0)
dtrain <- xgb.DMatrix(X, label = y)
rm(X, y)
sumpos <- dt[Y == "pos", .N]
sumneg <- dt[Y == "neg", .N]

tune.dt <- fread(paste0(fn.base, ".tune.tsv"), header = T, sep = "\t", quote = "\"", na.strings = "NA")
tune.dt <- tune.dt[order(-auc)]
bestTune <- tune.dt[1, ]
p <- list(eta = bestTune$eta, max_depth = bestTune$max_depth, gamma = bestTune$gamma, min_child_weight = bestTune$min_child_weight, subsample = bestTune$subsample, colsample_bytree = bestTune$colsample_bytree, objective = "binary:logistic", scale_pos_weight = sumneg/sumpos)
model <- xgb.train(params = p, data = dtrain, nrounds = bestTune$nrounds, verbose = 1, nthread = nthread)
xgb.save(model, paste0(fn.base, ".xgb.model"))
