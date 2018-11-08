#!/usr/bin/env Rscript

library(data.table)

homology <- fread("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data", header = F, sep = "\t", col.names = c("homologene_group_id", "tax_id", "ncbi_gene_id", "symbol", "protein_gi", "ref_seq"))

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
ncbi <- dbGetQuery(conn, "select gene_id,protein_id from ncbi")
dbDisconnect(conn)
rm(conn)
setDT(ncbi)

homology <- merge(homology, ncbi, by.x = "ncbi_gene_id", by.y = "gene_id", allow.cartesian = T)

fwrite(homology[, .(homologene_group_id, tax_id, protein_id)], file = "data/homology/homologene.tsv", sep = "\t", quote = T, col.names = T, row.names = F, na = "None")