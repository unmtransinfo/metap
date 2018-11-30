#!/usr/bin/env Rscript
###
# Plot multiple ROC curves for a given model and set of random seeds (normally 100).
# ROC curve with median AUC is highlighted.
###
library(data.table)
library(plotROC)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0) {
  fn <- args[1]
} else {
  fn <- "data/input/167000.rds"
}
fn.base <- sub(".rds", "", fn, fixed = T)
fn.base <- sub("input", "output", fn.base, fixed = T)

# File produced by xgboost_tune.R
auc.dt <- fread(paste0(fn.base, "/", "test.auc.tsv"), header = T, sep = "\t", quote = "\"", na.strings = "NA")
max.auc <- auc.dt[, max(auc)]
min.auc <- auc.dt[, min(auc)]
median.auc <- auc.dt[, median(auc)]
max.seed <- auc.dt[auc == max.auc, seed]
min.seed <- auc.dt[auc == min.auc, seed]
median.seed <- NA
if(nrow(auc.dt[auc == median.auc]) == 0) {
  auc.dt[, median.diff := abs(median.auc - auc)]
  min.diff <- auc.dt[, min(median.diff)]
  median.auc <- auc.dt[median.diff == min.diff, max(auc)]
  median.seed <- auc.dt[auc == median.auc, seed]
} else {
  median.seed <- auc.dt[auc == median.auc, seed][1]
}
#
# File produced by xgboost_tune.R
dt <- fread(paste0(fn.base, "/", "test.pred.tsv"), header = T, sep = "\t", quote = "\"", na.strings = "NA")
dt[, D := ifelse(Y == "pos", 1, 0)]
dt.median <- dt[seed == median.seed]

# All ROC curves in datatable, by seed.
roc.dt <- data.table()
for(s in dt[, unique(seed)]) {
  d <- calculate_roc(dt[seed == s, pred.prob], dt[seed == s, D])
  setDT(d)
  d[, seed := s]
  roc.dt <- rbindlist(list(roc.dt, d))
}
# TPF/FPF = True/False positive fraction
fwrite(roc.dt[, .(FPF, TPF, seed)], paste0(fn.base, "/", "test.roc.coord.tsv"), col.names = T, row.names = F, sep = "\t", quote = T, na = "NA")
#
# Compute AUC mean, SD.
m1 <- roc.dt[, .SD[which.max(TPF)], by = FPF]
m2 <- roc.dt[, .SD[which.min(TPF)], by = FPF]
setnames(m1, "TPF", "ymax")
setnames(m2, "TPF", "ymin")
ribbon.dt <- merge(m1, m2, by.x = "FPF", by.y = "FPF", sort = F)

sem <- sd(auc.dt[, auc])/sqrt(auc.dt[, .N])
lab <- sprintf("mean AUC = %.2f\u00B1%.2f\nAUC range = %.2f - %.2f", 
               mean(auc.dt[, auc]), sem*1.96,
               min.auc, max.auc)

plot <- ggplot(dt.median, aes(d = D, m = pred.prob)) + geom_roc(n.cuts = 0, size = 0.8, show.legend = F, color = "blue") + style_roc() +
  geom_ribbon(data = ribbon.dt, aes(ymin = ymin, ymax = ymax, x = FPF), inherit.aes = F, alpha = 0.2) +
  annotate("text", x = 0.62, y = 0.4, label = lab)

ggsave(paste0(fn.base, "/", "roc.pdf"), plot, dpi = 300, width = 8, height = 8)
saveRDS(plot, file = paste0(fn.base, "/", "roc.rda"))
