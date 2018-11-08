#!/usr/bin/env Rscript

library(xgboost)
library(data.table)
library(DMwR)
library(ggplot2)
library(RPostgreSQL)

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

dt <- readRDS(fn)

auc.test <- fread(paste0(fn.base, "/", "test.auc.tsv"), header = T, sep = "\t", quote = "\"")
setorder(auc.test, -"auc")
seed <- auc.test[1, seed]

model <- xgb.load(paste0(fn.base, "/", seed, ".model"))

features <- colnames(dt)[4:ncol(dt)]
features.names <- readable.feature.names(features)

imp <- xgb.importance(feature_names = features.names, model = model)
fwrite(imp, paste0(fn.base, "/", seed, ".varimp.tsv"), col.names = T, row.names = F, quote = T, sep	= "\t", na = "NA")
plot <- xgb.ggplot.importance(importance_matrix = imp, top_n = 20, n_clusters = 1) + theme_bw() + theme(legend.position = "none", axis.title = element_text(size = 26), axis.text = element_text(size = 14), plot.title = element_blank())
ggsave(paste0(fn.base, "/model.feature.importance.pdf"), plot = plot, dpi = 300, width = 9, height = 8)
