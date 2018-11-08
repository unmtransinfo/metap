#!/usr/bin/env Rscript

library(data.table)
library(xgboost)
library(xgboostExplainer)
library(fst)
library(RPostgreSQL)

nthread <- 86

readable.feature.names <- function(features) {
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
  protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
  drug_name <- dbGetQuery(conn, "select drug_id,drug_name from drug_name")
  kegg.pathways <- dbGetQuery(conn, "select kegg_pathway_id, kegg_pathway_name from kegg_pathway")
  reactome.pathways <- dbGetQuery(conn, "select reactome_id,name from reactome where species = 'Homo sapiens'")
  interpro <- dbGetQuery(conn, "select entry_ac,entry_name from interpro")
  ccle <- dbGetQuery(conn, "select distinct cell_id,tissue from ccle")
  dbDisconnect(conn)
  rm(conn)
  setDT(protein)
  setDT(drug_name)
  setDT(kegg.pathways)
  setDT(interpro)
  setDT(ccle)
  kegg.pathways <- unique(kegg.pathways, by = "kegg_pathway_id")
  setDT(reactome.pathways)
  features.names <- data.table(fid = features)
  features.names[fid %like% "pp\\.\\d+", c("protein_id") := tstrsplit(fid, ".", fixed = T, keep = 2L, type.convert = T)]
  features.names <- merge(features.names, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
  features.names[!is.na(symbol), fname := sprintf("PPI:%s", symbol)]
  features.names[, `:=`(protein_id = NULL, accession = NULL, symbol = NULL)]
  features.names[fid %like% "\\d+:[A-Z]", c("drug_id", "cell_id") := tstrsplit(fid, ":", fixed = T, keep = c(1L, 2L))]
  features.names[fid %like% "\\d+\\.P:[A-Z]", c("drug_id", "cell_id") := tstrsplit(fid, ":", fixed = T, keep = c(1L, 2L))]
  features.names <- merge(features.names, drug_name, by.x = "drug_id", by.y = "drug_id", all.x = T, sort = F)
  features.names[!is.na(drug_name), fname := sprintf("%s:%s signature", drug_name, cell_id)]
  features.names[, `:=`(drug_id = NULL, cell_id = NULL, drug_name = NULL)]
  features.names <- merge(features.names, kegg.pathways, by.x = "fid", by.y = "kegg_pathway_id", all.x = T, sort = F)
  features.names[!is.na(kegg_pathway_name), fname := kegg_pathway_name]
  features.names[, kegg_pathway_name := NULL]
  features.names <- merge(features.names, reactome.pathways, by.x = "fid", by.y = "reactome_id", all.x = T, sort = F)
  features.names[!is.na(name), fname := name]
  features.names[, name := NULL]
  features.names <- merge(features.names, interpro, by.x = "fid", by.y = "entry_ac", all.x = T, sort = F)
  features.names[!is.na(entry_name), fname := entry_name]
  features.names[, entry_name := NULL]
  ccle[is.na(tissue), `:=`(fid = cell_id, name = sprintf("expression in %s", cell_id))]
  ccle[!is.na(tissue), `:=`(fid = sprintf("%s_%s", cell_id, tissue), name = sprintf("%s (%s)", gsub("_", " ", tissue, fixed = T), cell_id))]
  features.names <- merge(features.names, ccle[, .(fid, name)], by.x = "fid", by.y = "fid", all.x = T, sort = F)
  features.names[!is.na(name), fname := name]
  features.names[, name := NULL]
  features.names[is.na(fname), fname := fid]
  return(features.names[, fname])
}

fn <- "data/input/nature_schizophrenia.rds"
fn.base <- sub(".rds", "", fn, fixed = T)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

dt <- readRDS(fn)
train <- dt[subset == "train"]
test <- dt[subset == "test"]
rm(dt)

X <- data.matrix(test[, -c("Y", "id1", "subset")])
y <- ifelse(test[, Y] == "pos", 1, 0)
dtest <- xgb.DMatrix(X, label = y)
rm(X, y)
meta.test <- test[, .(id1, Y, subset)]
setnames(meta.test, "id1", "protein_id")
meta.test <- merge(meta.test, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)

sumpos <- train[Y == "pos", .N]
sumneg <- train[Y == "neg", .N]

tune.dt <- fread(paste0(fn.base, ".tune.tsv"), header = T, sep = "\t", quote = "\"", na.strings = "NA")
tune.dt <- tune.dt[order(-auc)]
bestTune <- tune.dt[1, ]
p <- list(eta = bestTune$eta, max_depth = bestTune$max_depth, gamma = bestTune$gamma, min_child_weight = bestTune$min_child_weight, subsample = bestTune$subsample, colsample_bytree = bestTune$colsample_bytree, objective = "binary:logistic", scale_pos_weight = sumneg/sumpos)

for(k in 1:100) {
  seed <- 1000+k
  set.seed(seed)
  train.pos.idx <- sample(train[Y == "pos", id1], size = 0.8*train[Y == "pos", .N])
  train.neg.idx <- sample(train[Y == "neg", id1], size = 0.8*train[Y == "neg", .N])
  X <- data.matrix(train[id1 %in% c(train.pos.idx, train.neg.idx), -c("Y", "id1", "subset")])
  y <- ifelse(train[id1 %in% c(train.pos.idx, train.neg.idx), Y] == "pos", 1, 0)
  dtrain <- xgb.DMatrix(X, label = y)
  rm(X, y)
  
  meta.train <- train[id1 %in% c(train.pos.idx, train.neg.idx), .(id1, Y, subset)]
  setnames(meta.train, "id1", "protein_id")
  meta.train <- merge(meta.train, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
  
  
  model <- xgb.train(params = p, data = dtrain, nrounds = bestTune$nrounds, verbose = 1, nthread = nthread)
  explainer <- buildExplainer(xgb.model = model, type = "binary", base_score = 0.5, n_first_tree = model$niter - 1, trainingData = dtrain)
  pred.breakdown <- explainPredictions(xgb.model = model, explainer = explainer, data = dtest)
  pred.breakdown[, `:=`(protein_id = meta.test$protein_id, accession = meta.test$accession)]
  pred.breakdown <- melt(pred.breakdown, id.vars = c("protein_id", "accession"), variable.name = "feature", value.name = "log.odds", variable.factor = F, value.factor = F, verbose = T)
  pred.breakdown <- pred.breakdown[log.odds != 0.0]
  feature.name.long <- readable.feature.names(pred.breakdown[, feature])
  pred.breakdown[, feature.long := feature.name.long]
  pred.breakdown <- merge(pred.breakdown, protein[, .(protein_id, symbol)], by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
  write.fst(pred.breakdown, paste0(fn.base, "/", seed, ".test.pred.breakdown.fst"))
  cat("saved",paste0(fn.base, "/", seed, ".test.pred.breakdown.fst"),"\n")
  
  pred.breakdown <- explainPredictions(xgb.model = model, explainer = explainer, data = dtrain)
  pred.breakdown[, `:=`(protein_id = meta.train$protein_id, accession = meta.train$accession)]
  pred.breakdown <- melt(pred.breakdown, id.vars = c("protein_id", "accession"), variable.name = "feature", value.name = "log.odds", variable.factor = F, value.factor = F, verbose = T)
  pred.breakdown <- pred.breakdown[log.odds != 0.0]
  feature.name.long <- readable.feature.names(pred.breakdown[, feature])
  pred.breakdown[, feature.long := feature.name.long]
  pred.breakdown <- merge(pred.breakdown, protein[, .(protein_id, symbol)], by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
  write.fst(pred.breakdown, paste0(fn.base, "/", seed, ".train.pred.breakdown.fst"))
  cat("saved",paste0(fn.base, "/", seed, ".train.pred.breakdown.fst"),"\n")
}