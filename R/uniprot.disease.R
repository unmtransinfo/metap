#!/usr/bin/env Rscript

library(data.table)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

dt <- fread("data/uniprot/human.disease.tsv", header = T, sep = "\t", quote = "")
dt <- merge(dt, protein, by.x = "ACCESSION", by.y = "accession")

fwrite(dt[, .(protein_id, DISEASE_ID, DISEASE_TERM, DB_REF, REF_ID)], "data/uniprot/human.disease.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")