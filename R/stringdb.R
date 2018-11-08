#!/usr/bin/env Rscript

library(data.table)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein where accession is not null")
dbDisconnect(conn)
rm(conn)
setDT(protein)

stringdb <- fread("data/uniprot/uni2string.tsv", header = T, sep = "\t", quote = "")
stringdb <- merge(stringdb, protein, by.x = "ACCESSION", by.y = "accession")

fwrite(stringdb[, .(protein_id, STRINGDB_ID)], file = "data/uniprot/stringdb.mapping.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")