#!/usr/bin/env Rscript

library(data.table)
library(xgboost)
library(xgboostExplainer)
library(fst)
library(RPostgreSQL)

nthread <- 86
eps <- 0.9999999999999999

readable.feature.names <- function(features) {
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
  protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
  drug_name <- dbGetQuery(conn, "select drug_id,drug_name from drug_name")
  kegg.pathways <- dbGetQuery(conn, "select kegg_pathway_id, kegg_pathway_name from kegg_pathway")
  go <- dbGetQuery(conn, "select go_id, name from go")
  reactome.pathways <- dbGetQuery(conn, "select reactome_id,name from reactome where species = 'Homo sapiens'")
  interpro <- dbGetQuery(conn, "select entry_ac,entry_name from interpro")
  ccle <- dbGetQuery(conn, "select distinct cell_id,tissue from ccle")
  hpa <- dbGetQuery(conn, "select distinct tissue||'.'||cell_type||level as fid,tissue||' ('||cell_type||')' as name from hpa_norm_tissue")
  dbDisconnect(conn)
  rm(conn)
  setDT(protein)
  setDT(drug_name)
  setDT(kegg.pathways)
  setDT(interpro)
  setDT(ccle)
  setDT(go)
  setDT(hpa)
  hpa[, fid := gsub(" ", "_", fid, fixed = T)]
  hpa[, fid := gsub("-", "_", fid, fixed = T)]
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
  features.names <- merge(features.names, go, by.x = "fid", by.y = "go_id", all.x = T, sort = F)
  features.names[!is.na(name), fname := name]
  features.names[, name := NULL]
  features.names <- merge(features.names, hpa, by.x = "fid", by.y = "fid", all.x = T, sort = F)
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

fn <- commandArgs(T)[1]
fn.base <- sub(".rds", "", fn, fixed = T)
fn.base <- sub("input", "output", fn.base, fixed = T)

conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol,name from protein where accession is not null and tax_id = 9606")
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

tune.dt <- fread(paste0(fn.base, "/tune.tsv"), header = T, sep = "\t", quote = "\"", na.strings = "NA")
tune.dt <- tune.dt[order(-auc)]
bestTune <- tune.dt[1, ]
p <- list(eta = bestTune$eta, max_depth = bestTune$max_depth, gamma = bestTune$gamma, min_child_weight = bestTune$min_child_weight, subsample = bestTune$subsample, colsample_bytree = bestTune$colsample_bytree, objective = "binary:logistic", scale_pos_weight = sumneg/sumpos)

X <- data.matrix(train[, -c("Y", "id1", "subset")])
y <- ifelse(train[, Y] == "pos", 1, 0)
dtrain <- xgb.DMatrix(X, label = y)
rm(X, y)

meta.train <- train[, .(id1, Y, subset)]
setnames(meta.train, "id1", "protein_id")
meta.train <- merge(meta.train, protein, by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)

set.seed(1234)
model <- xgb.train(params = p, data = dtrain, nrounds = bestTune$nrounds, verbose = 1, nthread = nthread)
xgb.save(model, paste0(fn.base, "/global.model"))
explainer <- buildExplainer(xgb.model = model, type = "binary", base_score = 0.5, n_first_tree = model$niter - 1, trainingData = dtrain)

pred.breakdown <- explainPredictions(xgb.model = model, explainer = explainer, data = dtest)
pred.breakdown[, `:=`(protein_id = meta.test$protein_id, accession = meta.test$accession)]
pred.breakdown <- melt(pred.breakdown, id.vars = c("protein_id", "accession"), variable.name = "feature", value.name = "log.odds", variable.factor = F, value.factor = F, verbose = T)
pred.breakdown <- pred.breakdown[log.odds != 0.0 & !is.na(log.odds)]
if(nrow(pred.breakdown) == 0) {
	stop("pred.breakdown empty!")
}
pred.intercept <- pred.breakdown[feature == "intercept"]
pred.breakdown <- pred.breakdown[feature != "intercept", head(.SD[order(-abs(log.odds))], 50), by = protein_id]
pred.breakdown <- rbindlist(list(pred.intercept, pred.breakdown), use.names = T, fill = T)
feature.name.long <- readable.feature.names(pred.breakdown[, feature])
pred.breakdown[, feature.long := feature.name.long]
pred.breakdown <- merge(pred.breakdown, protein[, .(protein_id, symbol)], by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
write.fst(pred.breakdown, paste0(fn.base, "/blind.pred.breakdown.fst"))
cat("saved",paste0(fn.base, "/blind.pred.breakdown.fst"),"\n")
pred <- predict(model, dtest)
meta.test[, pred.prob := pred]
ecdf.fn <- ecdf(meta.test[, pred.prob])
meta.test[, zscore := qnorm(ifelse(ecdf.fn(pred.prob) >= eps, eps, ecdf.fn(pred.prob)))]
fwrite(meta.test, file = paste0(fn.base, "/", "blind.pred.tsv"), quote = T, sep = "\t", na = "NA", col.names = T, row.names = F)


pred.breakdown <- explainPredictions(xgb.model = model, explainer = explainer, data = dtrain)
pred.breakdown[, `:=`(protein_id = meta.train$protein_id, accession = meta.train$accession)]
pred.breakdown <- melt(pred.breakdown, id.vars = c("protein_id", "accession"), variable.name = "feature", value.name = "log.odds", variable.factor = F, value.factor = F, verbose = T)
pred.breakdown <- pred.breakdown[log.odds != 0.0]
pred.intercept <- pred.breakdown[feature == "intercept"]
pred.breakdown <- pred.breakdown[feature != "intercept", head(.SD[order(-abs(log.odds))], 50), by = protein_id]
pred.breakdown <- rbindlist(list(pred.intercept, pred.breakdown), use.names = T, fill = T)
feature.name.long <- readable.feature.names(pred.breakdown[, feature])
pred.breakdown[, feature.long := feature.name.long]
pred.breakdown <- merge(pred.breakdown, protein[, .(protein_id, symbol)], by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
write.fst(pred.breakdown, paste0(fn.base, "/train.pred.breakdown.fst"))
cat("saved",paste0(fn.base, "/train.pred.breakdown.fst"),"\n")
pred <- predict(model, dtrain)
meta.train[, pred.prob := pred]
ecdf.fn <- ecdf(meta.train[, pred.prob])
meta.train[, zscore := qnorm(ifelse(ecdf.fn(pred.prob) >= eps, eps, ecdf.fn(pred.prob)))]
fwrite(meta.train, file = paste0(fn.base, "/train.pred.tsv"), quote = T, sep = "\t", na = "NA", col.names = T, row.names = F)
