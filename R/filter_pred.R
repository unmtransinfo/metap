#!/usr/bin/env Rscript

library(data.table)
library(fst)

fn <- "data/input/nature_schizophrenia.rds"
#fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)
ref.protein <- "P14416"

test.auc <- fread(input = paste0(fn.base, ".test.auc.tsv"), header = T, sep = "\t", quote = "\"")

features <- 30
metrics.dt <- data.table()
for(seed in test.auc[, seed]) {
  train.dt <- read_fst(paste0(fn.base, "/", seed, ".train.pred.breakdown.fst"), as.data.table = T)
  train.dt <- train.dt[accession == ref.protein & feature != "intercept"][order(-abs(log.odds))]
  test.dt <- read_fst(paste0(fn.base, "/", seed, ".test.pred.breakdown.fst"), as.data.table = T)
  test.pred <- fread(paste0(fn.base, "/", seed, ".test.pred.tsv"), sep = "\t", header = T, quote = "\"", na.strings = "NA")
  test.pred <- test.pred[order(-pred.prob)]
  test.pred[, rank := 1:.N]
  
  for(a in test.dt[, unique(accession)]) {
    if(a == ref.protein) {
      next
    }
    dt <- test.dt[accession == a][order(-abs(log.odds))][1:features]
    x <- train.dt[1:features]
    x <- merge(x, dt[, .(feature.long, log.odds)], by.x = "feature.long", by.y = "feature.long", all.x = T, sort = F)
    complete.obs <- nrow(na.omit(x))
    if(complete.obs < 3) {
      next
    }
    tau <- cor(x[, log.odds.x], x[, log.odds.y], method = "kendall", use = "complete.obs")
    record <- data.table(ref.protein = ref.protein, target.protein = a, tau = tau, complete.obs = complete.obs, seed = seed, pred.prob = test.pred[accession == a, pred.prob], rank = test.pred[accession == a, rank], model.weight = test.pred[accession == a, weight])
    metrics.dt <- rbindlist(list(metrics.dt, record), use.names = T)
  }
}

metrics.dt.weighted <- metrics.dt[, .(rank.weighted = weighted.mean(rank, model.weight), tau.weighted = weighted.mean(tau, model.weight), pred.weighted = weighted.mean(pred.prob, model.weight), obs.weighted = weighted.mean(complete.obs, model.weight)), by = target.protein]
fwrite(metrics.dt.weighted, paste0(fn.base, "/", ref.protein, ".scored.tsv"), col.names = T, row.names = F, quote = T, sep = "\t", na = "NA")
