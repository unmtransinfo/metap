#!/usr/bin/env Rscript

library(data.table)

entry.list <- fread("ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list", header = T, sep = "\t", quote = "")
fwrite(entry.list, file = "data/interpro/interpro.tsv", col.names = T, row.names = F, quote = T, sep = "\t", na = "None")

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), dbname = "metap_dev")
protein <- dbGetQuery(conn, "select protein_id,accession from protein")
dbDisconnect(conn)
rm(conn)
setDT(protein)

download.file("ftp://ftp.ebi.ac.uk/pub/databases/interpro/current/protein2ipr.dat.gz", destfile = "data/interpro/protein2ipr.dat.gz")
annon <- fread('gunzip -c data/interpro/protein2ipr.dat.gz', header = F, select = c(1,2), sep = "\t")
setnames(annon, c("V1", "V2"), c("accession", "ENTRY_AC"))

annon <- merge(annon, protein, by.x = "accession", by.y = "accession")
fwrite(unique(annon[, .(protein_id, ENTRY_AC)]), file = "data/interpro/interpro_annon.tsv",  quote = T, sep = "\t", na = "None", col.names = T, row.names = F)
