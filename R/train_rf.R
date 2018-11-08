#!/usr/bin/env Rscript

library(data.table)
library(randomForestSRC)
library(ggRandomForests)
library(pROC)
set.seed(1234)

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
neg.sample.idx <- sample(dt[Y == "neg", which = T], size = length(pos.sample.idx)*2, replace = F)
train.pos.sample.idx <- sample(pos.sample.idx, size = round(length(pos.sample.idx)*0.8), replace = F)
train.neg.sample.idx <- sample(neg.sample.idx, size = round(length(neg.sample.idx)*0.8), replace = F)
test.pos.sample.idx <- pos.sample.idx[!pos.sample.idx %in% train.pos.sample.idx]
test.neg.sample.idx <- neg.sample.idx[!neg.sample.idx %in% train.neg.sample.idx]

rf.model <- rfsrc(Y ~ ., data = dt[c(train.pos.sample.idx, train.neg.sample.idx)], ntree = 3000, mtry = 100, importance = T)
pred <- predict.rfsrc(rf.model, newdata = dt[c(test.pos.sample.idx, test.neg.sample.idx)])
var.imp <- gg_vimp(rf.model)
setDT(var.imp)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select 'pp.'||protein_id as protein_id,'PP:'||symbol as gene_sym from protein")
reactome <- dbGetQuery(conn, "select overlay(overlay(reactome_id PLACING '.' from 2 for 1) PLACING '.' FROM 6 for 1) as reactome_id,name as pathway from reactome where species = 'Homo sapiens'")
go <- dbGetQuery(conn, "select overlay(go_id PLACING '_' FROM 3 FOR 1) as go_id,name as go_term from go")
dbDisconnect(conn)
rm(conn)
setDT(protein)
setDT(reactome)
setDT(go)

var.imp <- merge(var.imp, protein, by.x = "vars", by.y = "protein_id", all.x = T)
var.imp <- merge(var.imp, reactome, by.x = "vars", by.y = "reactome_id", all.x = T)
var.imp <- merge(var.imp, go, by.x = "vars", by.y = "go_id", all.x = T)
var.imp[!is.na(gene_sym), vars2 := gene_sym]
var.imp[!is.na(pathway), vars2 := pathway]
var.imp[!is.na(go_term), vars2 := go_term]
var.imp[is.na(vars2), vars2 := vars]

vars.sel <- var.imp[set == "all"][order(-vimp)][, unique(vars)]
vars.sel <- vars.sel[1:20]
vars2.sel <- var.imp[set == "all"][order(-vimp)][, unique(vars2)]
vars2.sel <- vars2.sel[1:20]



var.imp <- gg_vimp(rf.model)
var.imp <- filter(var.imp, vars %in% vars.sel)
vars.lbls.map <- data.table(vars = vars.sel, vars2 = vars2.sel)
var.lbls <- data.table(vars = var.imp$vars, idx = 1:length(var.imp$vars))
var.lbls <- merge(var.lbls, vars.lbls.map, by.x = "vars", by.y = "vars", all.x = T)
var.lbls <- var.lbls[order(idx)]
use.lbls <- var.lbls$vars2
names(use.lbls) <- var.lbls$vars

class(var.imp) <- c("gg_vimp", "data.frame")
p <- plot(var.imp, lbls = use.lbls) + theme_bw() + theme(legend.position = "none", axis.text.y = element_text(size = 16))
ggsave(plot = p, filename = "vimp.png", dpi = 300, width = 12, height = 9)
