#!/usr/bin/env Rscript

library(data.table)
library(readxl)
library(RPostgreSQL)

conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol,name from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)


pheno.list <- read_xlsx("data/output/pheno_list.xlsx", sheet = 1)
setDT(pheno.list)

train.gene.list <- data.table()
predictions <- data.table()

for(i in 1:nrow(pheno.list)) {
  if(!file.exists(paste0("data/input/", pheno.list[i, phenotype_id], ".rds", collapse = ""))) {
    next
  }
  dt <- readRDS(paste0("data/input/", pheno.list[i, phenotype_id], ".rds", collapse = ""))
  dt <- dt[Y == "pos" & subset == "train", .(id1)]
  dt <- merge(dt, protein, by.x = "id1", by.y = "protein_id")
  dt[, phenotype_id := pheno.list[i, phenotype_id]]
  dt[, phenotype := pheno.list[i, phenotype_term]]
  train.gene.list <- rbindlist(list(train.gene.list, dt), use.names = T, fill = T)
  if(!file.exists(paste0("data/output/", pheno.list[i, phenotype_id], "/blind.pred.tsv", collapse = ""))) {
    next
  }
  blind <- fread(paste0("data/output/", pheno.list[i, phenotype_id], "/blind.pred.tsv", collapse = ""), header = T, sep = "\t", quote = "\"", na.strings = "NA")
  blind[, phenotype_id := pheno.list[i, phenotype_id]]
  blind[, phenotype := pheno.list[i, phenotype_term]]
  predictions <- rbindlist(list(predictions, blind), use.names = T, fill = T)
  cat(pheno.list[i, phenotype_term], "\n")
}

fwrite(train.gene.list, "phenotype_genes_train.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "NA")
fwrite(predictions, "phenotype_predictions.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "NA")