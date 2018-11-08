#!/usr/bin/env Rscript

library(data.table)
library(pROC)
library(ggplot2)
library(plotROC)
library(caret)
library(xgboost)

#fn <- commandArgs(T)[1]
fn <- "data/input/nature_schizophrenia.rds"
fn.base <- sub(".rds", "", fn, fixed = T)
fn.base <- sub("/input/", "/output/", fn.base, fixed = T)

#cv.pred <- fread(paste0(fn.base, ".weighted.cv.pred.tsv"), header = T, sep = "\t", quote = "\"")
cv.auc <- fread(paste0(fn.base, "/", "test.auc.tsv"), header = T, sep = "\t")
setorder(cv.auc, -auc)
seed <- cv.auc[1, seed]
model <- xgb.load(paste0(fn.base, "/", seed, ".model"))

dt <- readRDS(fn)
dt <- dt[subset == "train"]
set.seed(seed)
train.pos.idx <- sample(dt[Y == "pos", id1], size = 0.8*dt[Y == "pos", .N])
train.neg.idx <- sample(dt[Y == "neg", id1], size = 0.8*dt[Y == "neg", .N])
X <- data.matrix(dt[!id1 %in% c(train.pos.idx, train.neg.idx), -c("Y", "id1", "subset")])
y <- ifelse(dt[!id1 %in% c(train.pos.idx, train.neg.idx), Y] == "pos", 1, 0)
test.pred <- dt[!id1 %in% c(train.pos.idx, train.neg.idx), .(id1, Y)]
test.pred[, actual := y]
#test.pred <- merge(test.pred, protein, by.x = "id1", by.y = "protein_id", all.x = T, sort = F)
dtest <- xgb.DMatrix(X, label = y)
rm(X, y)
pred <- predict(model, dtest)
test.pred[, pred.prob := pred]
auc <- auc(test.pred$actual, test.pred$pred.prob)


p <- ggplot(test.pred, aes(d = actual, m = pred.prob)) + geom_roc(n.cuts = 0, size = 0.5) + style_roc() + annotate("text", x = 0.28, y = 0.55, label = sprintf("AUC = %.2f", auc), color = "#E69F00", size = 5) + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_blank())
ggsave(paste0(fn.base, "/", "auc.plot.png"), plot = p, dpi = 300, width = 11, height = 8)

rocobj <- roc(test.pred$actual, test.pred$pred.prob)
cutoff <- coords(rocobj, "b", ret="t", best.method="closest.topleft")
predicted.label <- ifelse(test.pred$pred.prob >= cutoff, "pos", "neg")
predicted.label <- as.factor(predicted.label)
sink(paste0(fn.base, "/", "conf.matrix.tsv"))
confusionMatrix(predicted.label, test.pred$Y, positive = "pos")
sink()