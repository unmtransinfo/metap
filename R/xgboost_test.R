#!/usr/bin/env Rscript

library(data.table)
library(readxl)
library(plotROC)
library(pROC)
library(caret)

fn <- "data/input/104300.rds"
fn.base <- sub(".rds", "", fn, fixed = T)

dt <- readRDS(fn)

test <- read_xlsx("data/input/104300_text_mining.xlsx", sheet = 1)
setDT(test)
test <- test[include == T]

neg.idx <- readRDS("data/input/104300.test.neg.idx.rds")

neg <- dt[id1 %in% neg.idx]
pos <- dt[id1 %in% test$protein_id & subset != "train"]
pos[, Y := "pos"]

test <- rbindlist(list(pos, neg))
x <- test[, -c("id1", "Y", "subset")]
y <- as.numeric(test[, Y])
y[y == 1] <- 0
y[y == 2] <- 1
test.meta <- test[, .(id1, Y)]
test <- xgb.DMatrix(data.matrix(x), label = y)
rm(x)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

model.weighted <- xgb.load(paste0(fn.base, ".xgb.weighted.model"))
pred <- predict(model.weighted, test)
test.meta[, pred.prob := pred]
test.meta[, zscore := scale(pred.prob)]
setnames(test.meta, "id1", "protein_id")
test.meta <- merge(test.meta, protein, by.x = "protein_id", by.y = "protein_id")
y <- as.numeric(test.meta[, Y])
y[y == 1] <- 0
y[y == 2] <- 1

fwrite(test.meta[order(-zscore)], file = paste0(fn.base, ".xgb.pred.test.weighted.model.tsv"), sep = "\t", col.names = T, row.names = F, quote = T)

weighted.auc <- readRDS(paste0(fn.base, ".xgb.auc.weighted.sample.rds"))
weighted.auc[, Method := "Cross-validation"]
test.auc <- data.table(actual = y, predicted = test.meta$pred.prob, Method = "Test set")
weighted.auc <- rbindlist(list(weighted.auc, test.auc))
auc.cv <- auc(weighted.auc[Method == "Cross-validation", actual], weighted.auc[Method == "Cross-validation", predicted])
auc.test <- auc(weighted.auc[Method == "Test set", actual], weighted.auc[Method == "Test set", predicted])
p <- ggplot(weighted.auc, aes(d = actual, m = predicted, color = Method)) + geom_roc(n.cuts = 0, size = 0.5) + style_roc() + annotate("text", x = 0.25, y = 0.78, label = sprintf("AUC = %.2f", auc.cv), color = "#E69F00", size = 5) + annotate("text", x = 0.43, y = 0.55, label = sprintf("AUC = %.2f", auc.test), color = "#56B4E9", size = 5) + scale_color_manual(values=c("#E69F00", "#56B4E9")) + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_blank())
ggsave(paste0(fn.base, ".xgb.auc.weighted.sample.png"), plot = p, dpi = 300, width = 11, height = 8)
rocobj <- roc(test.meta$Y, test.meta$pred.prob)
cutoff <- coords(rocobj, "b", ret="t", best.method="closest.topleft")
predicted.label <- ifelse(test.meta$pred.prob >= cutoff, "pos", "neg")
sink(paste0(fn.base, ".xgb.conf.matrix.test.weighted.model.tsv"))
sink()

model.balanced <- xgb.load(paste0(fn.base, ".xgb.balanced.model"))
pred <- predict(model.balanced, test)
test.meta[, pred.prob := pred]
test.meta[, zscore := scale(pred.prob)]
fwrite(test.meta[order(-zscore)], file = paste0(fn.base, ".xgb.pred.test.balanced.model.tsv"), sep = "\t", col.names = T, row.names = F, quote = T)

balanced.auc <- readRDS(paste0(fn.base, ".xgb.auc.balanced.sample.rds"))
balanced.auc[, Method := "Cross-validation"]
test.auc <- data.table(actual = y, predicted = test.meta$pred.prob, Method = "Test set")
balanced.auc <- rbindlist(list(balanced.auc, test.auc))
auc.cv <- auc(balanced.auc[Method == "Cross-validation", actual], balanced.auc[Method == "Cross-validation", predicted])
auc.test <- auc(balanced.auc[Method == "Test set", actual], balanced.auc[Method == "Test set", predicted])
p <- ggplot(balanced.auc, aes(d = actual, m = predicted, color = Method)) + geom_roc(n.cuts = 0, size = 0.5) + style_roc() + annotate("text", x = 0.1, y = 0.8, label = sprintf("AUC = %.2f", auc.cv), color = "#E69F00", size = 5) + annotate("text", x = 0.12, y = 0.55, label = sprintf("AUC = %.2f", auc.test), color = "#56B4E9", size = 5) + scale_color_manual(values=c("#E69F00", "#56B4E9")) + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_blank())
ggsave(paste0(fn.base, ".xgb.auc.balanced.sample.png"), plot = p, dpi = 300, width = 11, height = 8)
rocobj <- roc(test.meta$Y, test.meta$pred.prob)
cutoff <- coords(rocobj, "b", ret="t", best.method="closest.topleft")
predicted.label <- ifelse(test.meta$pred.prob >= cutoff, "pos", "neg")
sink(paste0(fn.base, ".xgb.conf.matrix.test.balanced.model.tsv"))
confusionMatrix(predicted.label, test.meta$Y, positive = "pos")
sink()

