#!/usr/bin/env Rscript

library(data.table)
library(xgboost)
library(DMwR)
library(xgboostExplainer)
library(RPostgreSQL)

readable.feature.names <- function(features) {
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
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
  ccle[!is.na(tissue), `:=`(fid = sprintf("%s_%s", cell_id, tissue), name = sprintf("expression in %s (%s)", cell_id, tissue))]
  features.names <- merge(features.names, ccle[, .(fid, name)], by.x = "fid", by.y = "fid", all.x = T, sort = F)
  features.names[!is.na(name), fname := name]
  features.names[, name := NULL]
  features.names[is.na(fname), fname := fid]
  return(features.names[, fname])
}

eps <- 0.9999999999999999

#fn <- commandArgs(T)[1]
fn <- "data/input/nature_schizophrenia.rds"
fn.base <- sub(".rds", "", fn, fixed = T)

model.param <- fread(paste0(fn.base, ".xgb.tune.weighted.tsv"), header = T, sep = "\t", quote = "\"")
bestTune <- model.param[1, ]
rm(model.param)


dt <- readRDS(fn)

X <- data.matrix(dt[subset == "train", -c("Y", "id1", "subset")])
y <- ifelse(dt[subset == "train", Y] == "pos", 1, 0)
train <- xgb.DMatrix(X, label = y)
sumpos <- sum(y == 1)
sumneg <- sum(y == 0)
rm(X)
rm(y)

x <- dt[subset == "test", -c("id1", "Y", "subset")]
y <- ifelse(dt[subset == "test", Y] == "pos", 1, 0)
test.meta <- dt[subset == "test", .(id1, Y)]
test <- xgb.DMatrix(data.matrix(x), label = y)
rm(x)
rm(y)

p <- list(eta = bestTune$eta, max_depth = bestTune$max_depth, gamma = bestTune$gamma, min_child_weight = bestTune$min_child_weight, subsample = bestTune$subsample, colsample_bytree = bestTune$colsample_bytree, objective = "binary:logistic", scale_pos_weight = sumneg/sumpos)
model <- xgb.train(params = p, data = train, nrounds = bestTune$nrounds, verbose = 1, nthread = 4)
pred <- predict(model, test)
test.meta[, pred.prob := pred]
ecdf.fn <- ecdf(test.meta[, pred.prob])
test.meta[, zscore := qnorm(ifelse(ecdf.fn(pred.prob) >= eps, eps, ecdf.fn(pred.prob)))]

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

test.meta <- merge(test.meta, protein, by.x = "id1", by.y = "protein_id", all.x = T, sort = F)
setnames(test.meta, "id1", "protein_id")
fwrite(test.meta[order(-zscore)], file = paste0(fn.base, ".xgb.pred.weighted.model.tsv"), sep = "\t", col.names = T, row.names = F, quote = T)

explainer <- buildExplainer(xgb.model = model, type = "binary", base_score = 0.5, n_first_tree = model$niter - 1, trainingData = train)
pred.breakdown <- explainPredictions(xgb.model = model, explainer = explainer, data = test)
pred.breakdown[, `:=`(protein_id = test.meta$protein_id, accession = test.meta$accession)]
pred.breakdown <- melt(pred.breakdown, id.vars = c("protein_id", "accession"), variable.name = "feature", value.name = "log.odds", variable.factor = F, value.factor = F, verbose = T)
pred.breakdown <- pred.breakdown[log.odds != 0.0]
feature.name.long <- readable.feature.names(pred.breakdown[, feature])
pred.breakdown[, feature.long := feature.name.long]
pred.breakdown <- merge(pred.breakdown, protein[, .(protein_id, symbol)], by.x = "protein_id", by.y = "protein_id", all.x = T, sort = F)
fwrite(pred.breakdown, file = paste0(fn.base, ".xgb.balanced.pred.breakdown.tsv"), col.names = T, row.names = F, sep = "\t", quote = T, na = "NA")

