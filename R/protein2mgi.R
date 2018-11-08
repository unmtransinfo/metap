#!/usr/bin/env Rscript

library(data.table)
library(RPostgreSQL)

uni2mgi <- fread("data/uniprot/uni2mgi.tsv", header = T, sep = "\t", quote = "")

conn <- dbConnect(PostgreSQL(), dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein where tax_id = 10090")
dbDisconnect(conn)
rm(conn)

setDT(protein)

uni2mgi <- merge(uni2mgi, protein, by.x = "ACCESSION", by.y = "accession")

fwrite(uni2mgi[, .(protein_id, MGI)], "data/uniprot/protein2mgi.tsv", col.names = T, row.names = F, quote = T, na = "None", sep = "\t")