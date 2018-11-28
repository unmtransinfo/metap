#!/usr/bin/env Rscript

library(data.table)
library(waterfalls)
library(ggplot2)

fn <- "data/input/104300.rds"
fn.base <- sub(".rds", "", fn, fixed = T)

plot.sample.size <- 30

pred.breakdown <- fread(paste0(fn.base, ".xgb.weighted.pred.breakdown.tsv"), header = T, sep = "\t", quote = "\"", na.strings = "NA")

sample <- pred.breakdown[protein_id == 18706, .(log.odds, feature.long)]
sample[, log.odds := round(log.odds, 2)]
sample <- sample[order(-abs(log.odds))]
idx.intercept <- sample[feature.long == "intercept", which = T]
idx.keep <- unique(c(idx.intercept, 1:plot.sample.size))
sample <- sample[idx.keep]
plot <- waterfall(sample, calc_total = T, rect_text_size = 1.0) + theme_bw() + theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 20), axis.text.y = element_text(size = 20) , axis.title = element_text(size = 26)) + xlab("Feature") + ylab("Log odds")
ggsave("data/input/104300.SCGB3A1.feature.breackdown.png", dpi = 300, width = 18, height = 16)
