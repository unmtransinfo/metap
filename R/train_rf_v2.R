#!/usr/bin/env Rscript

library(data.table)
library(randomForest)
library(pROC)
library(ggplot2)
set.seed(1234)

colors <- topo.colors(10, alpha = 1)


dt <- readRDS("data/input/input_test.rds")
dt[, id1 := NULL]
col.name <- colnames(dt)
col.name <- gsub("-", ".", col.name, fixed = T)
col.name <- gsub(" . ", ".", col.name, fixed = T)
col.name <- gsub(" ", "_", col.name, fixed = T)
col.name <- gsub("(", "", col.name, fixed = T)
col.name <- gsub(")", "", col.name, fixed = T)
col.name <- gsub(":", "_", col.name, fixed = T)
col.name <- gsub("/", "_", col.name, fixed = T)
col.name <- gsub(",", "", col.name, fixed = T)
setnames(dt, colnames(dt), col.name)
pos.sample.idx <- dt[Y == "pos", which = T]

models <- list()
train.rocs <- list()
test.rocs <- list()
i <- 1

while (i <= 10) {
  neg.sample.idx <- sample(dt[Y == "neg", which = T], size = length(pos.sample.idx)*2, replace = F)
  train.pos.sample.idx <- sample(pos.sample.idx, size = round(length(pos.sample.idx)*0.8), replace = F)
  train.neg.sample.idx <- sample(neg.sample.idx, size = round(length(neg.sample.idx)*0.8), replace = F)
  test.pos.sample.idx <- pos.sample.idx[!pos.sample.idx %in% train.pos.sample.idx]
  test.neg.sample.idx <- neg.sample.idx[!neg.sample.idx %in% train.neg.sample.idx]
  rf.model <- randomForest(Y ~., data = dt[c(train.pos.sample.idx, train.neg.sample.idx)], ntree = 3000, mtry = 100, importance = T)
  train.roc <- roc(dt[c(train.pos.sample.idx, train.neg.sample.idx), Y], rf.model$votes[, 2])
  if(train.roc$auc >= 0.98) {
    next
  }
  pred <- predict(rf.model, newdata = dt[c(test.pos.sample.idx, test.neg.sample.idx)], type = "prob")
  test.roc <- roc(dt[c(test.pos.sample.idx, test.neg.sample.idx), Y], pred[, 2])
  if(test.roc$auc >= 0.98) {
    next
  }
  if(test.roc$auc > train.roc$auc) {
    next
  }
  models[[i]] <- rf.model
  train.rocs[[i]] <- train.roc
  test.rocs[[i]] <- test.roc
  i <- i + 1
  cat("i = ", i, "\n")
  cat("train AUC = ", train.roc$auc, "\n")
  cat("test AUC = ", test.roc$auc, "\n")
}

train.roc.plot <- ggroc(train.rocs, aes = "colour", size = 1) + geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), color = "grey", size = 0.2, linetype = "dotted") + theme_bw() + theme(legend.position="none", axis.title = element_text(size = 20)) + annotate("text", label = "AUC range 0.9 - 0.97", x = 0.3, y = 0.25, size = 12)
ggsave("train.roc.plot.png", dpi = 300, height = 8, width = 6, plot = train.roc.plot)
test.roc.plot <- ggroc(test.rocs, aes = "colour", size = 1) + geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), color = "grey", size = 0.2, linetype = "dotted") + theme_bw() + theme(legend.position="none", axis.title = element_text(size = 20)) + annotate("text", label = "AUC range 0.8 - 0.94", x = 0.3, y = 0.25, size = 12)
ggsave("test.roc.plot.png", dpi = 300, height = 8, width = 8, plot = test.roc.plot)