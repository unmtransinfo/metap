#!/usr/bin/env Rscript

library(data.table)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein where accession is not null")
dbDisconnect(conn)
rm(conn)
setDT(protein)

ncbi <- fread("data/uniprot/human.ncbi.tsv", header = T, sep = "\t")
ncbi <- rbindlist(list(ncbi, fread("data/uniprot/mouse.ncbi.tsv", header = T, sep = "\t", na.strings = "", quote = "")), use.names = T)
ncbi <- rbindlist(list(ncbi, fread("data/uniprot/rat.ncbi.tsv", header = T, sep = "\t", na.strings = "", quote = "")), use.names = T)


ncbi <- unique(ncbi)
ncbi <- merge(ncbi, protein, by.x = "ACCESSION", by.y = "accession")

fwrite(ncbi[, .(protein_id, NCBI_GENEID)], file = "data/uniprot/ncbi.mapping.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")