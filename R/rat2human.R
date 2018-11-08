#!/usr/bin/env Rscript

library(data.table)
library(tidyr)

ortho <- fread("ftp://ftp.rgd.mcw.edu/pub/data_release/RGD_ORTHOLOGS.txt", header = T, sep = "\t", quote = "", na.strings = "", skip = 52, select = c(1:7), colClasses = c("charater", "integer", "integer", "character", "character", "character", "character", "character", "character", "character", "character", "character", "character"))
ortho <- ortho[!is.na(HUMAN_ORTHOLOG_NCBI_GENE_ID)]
ortho <- separate_rows(ortho, HUMAN_ORTHOLOG_NCBI_GENE_ID, sep = "\\|", convert = T)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
rat.ncbi <- dbGetQuery(conn, "select ncbi.protein_id,gene_id from ncbi,protein where protein.protein_id=ncbi.protein_id and protein.tax_id = '10116'")
human.ncbi <- dbGetQuery(conn, "select ncbi.protein_id,gene_id from ncbi,protein where protein.protein_id=ncbi.protein_id and protein.tax_id = '9606'")
dbDisconnect(conn)
rm(conn)
setDT(rat.ncbi)
setDT(human.ncbi)

setnames(rat.ncbi, "protein_id", "rat_protein_id")
setnames(human.ncbi, "protein_id", "human_protein_id")

ortho <- merge(ortho, rat.ncbi, by.x = "RAT_GENE_NCBI_GENE_ID", by.y = "gene_id")
ortho <- merge(ortho, human.ncbi, by.x = "HUMAN_ORTHOLOG_NCBI_GENE_ID", by.y = "gene_id")

fwrite(unique(ortho[, .(rat_protein_id, human_protein_id)]), "data/rgd/rat2human.tsv", col.names = T, row.names = F, sep = "\t", na = "None")