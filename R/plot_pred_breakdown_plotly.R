#!/usr/bin/env Rscript
###
# For a given model and gene, generates plot of top features ranked by log odds.
# Cumulative and individual log.odds are shown, via position and height of bars, respectively.
###
library(data.table)
library(fst)
library(plotly)

fn <- "data/input/104300.rds"
#fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)
fn.base <- sub("input", "output", fn.base, fixed = T)

pred.breakdown <- read_fst(paste0(fn.base, "/test.pred.breakdown.fst"), as.data.table = T)

gene <- "LILRA3"

x <- pred.breakdown[symbol == gene][feature != "intercept"][order(-abs(log.odds))][1:20]
x[, text := sprintf("%.2f", log.odds)] #round for plot
#
# Separate features by log.odds sign (+/-).
pos <- x[log.odds > 0][order(-log.odds)]
neg <- x[log.odds < 0][order(log.odds)]
x <- rbindlist(list(pos,neg), use.names = T, fill = T)
#
# Compute cumulative sum of log.odds.
x[, base := cumsum(log.odds)]
x[, base := base - log.odds]
#
x[, pos := ifelse(log.odds > 0, log.odds, 0)]
x[, neg := ifelse(log.odds < 0, log.odds, 0)]
x[, `:=`(protein_id = NULL, feature = NULL)]
x[, neg := abs(neg)]
x[neg > 0, base := base - neg]
x[, feature.long := factor(feature.long, levels = feature.long, ordered = T)]
#
dt <- melt(x, id.vars = "feature.long", measure.vars = c("base", "pos", "neg"), 
           variable.name = "subset", value.name = "log.odds")
dt <- rbindlist(list(dt, data.table(feature.long = "Total", subset = "total", 
          log.odds=dt[subset=="pos", sum(log.odds)] - dt[subset=="neg", sum(log.odds)])), use.names=T, fill=T)
dt[, subset := factor(subset, levels=c("pos", "neg", "base", "total"), ordered = T)]

p <- ggplot(data = dt, aes(x = feature.long, y = log.odds, fill = subset, label = sprintf("%.2f", log.odds))) +
  geom_col() + theme_classic() + scale_fill_manual(values=c("#009E73", "#D55E00", "#FFFFFF", "#A9A9A9")) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 60, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), axis.title = element_text(size = 20), 
        plot.title = element_text(size = 18, hjust = 0.5)) +
  ggtitle(paste0("Top 20 important features for ", gene, " association with Alzheimer's Disease")) + 
  xlab("Feature") + ylab("Log odds") +
  geom_label(aes(label=ifelse(log.odds>0 & subset!="base", 
                  ifelse(subset=="neg", round(-1*log.odds, 2), round(log.odds, 2)), NA)), 
             position=position_stack(vjust=0.5))
#
ggsave("pred.breakdown.png", p, dpi = 300, width = 16, height = 12)
#

#
p <- plot_ly(data=x, type='bar', x=~feature.long, y=~base, marker=list(color='rgba(1,1,1, 0.0)'), hoverinfo='none') %>%
  add_trace(y = ~pos, marker = list(color = 'rgba(50, 171, 96, 0.7)',
                                        line = list(color = 'rgba(50, 171, 96, 1.0)',
                                                    width = 2))) %>%
  add_trace(y = ~neg, marker = list(color = 'rgba(219, 64, 82, 0.7)',
                                      line = list(color = 'rgba(219, 64, 82, 1.0)',
                                                  width = 2))) %>%
  layout(title = gene,
         xaxis = list(title = "Model features"),
         yaxis = list(title = "Log odds"),
         barmode = 'stack',
         paper_bgcolor = 'rgba(245, 246, 249, 1)',
         plot_bgcolor = 'rgba(245, 246, 249, 1)',
         showlegend = FALSE) %>%
  add_annotations(text = x$text,
                  x = x$feature.long,
                  y = x$base + 0.08,
                  xref = "x",
                  yref = "y",
                  font = list(family = 'Arial',
                              size = 14,
                              color = 'rgba(245, 246, 249, 1)'),
                  showarrow = FALSE)
