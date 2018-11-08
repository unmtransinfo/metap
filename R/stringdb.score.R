#!/usr/bin/env Rscript

library(data.table)
library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,stringdb_id from stringdb")
dbDisconnect(conn)
rm(conn)
setDT(protein)

download.file("https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz", destfile = "data/stringdb/9606.protein.links.v10.5.txt.gz")
string <- fread("gunzip -c data/stringdb/9606.protein.links.v10.5.txt.gz", header = T, sep = " ", quote = "")

string[, key := sprintf("%s.%s", ifelse(protein1 >= protein2, protein1, protein2), ifelse(protein2 < protein1, protein2, protein1))]
string <- unique(string, by = "key")
string[, key := NULL]
string <- merge(string, protein, by.x = "protein1", by.y = "stringdb_id")
setnames(string, "protein_id", "protein_id1")
string <- merge(string, protein, by.x = "protein2", by.y = "stringdb_id")
setnames(string, "protein_id", "protein_id2")
fwrite(string[, .(protein_id1, protein_id2, combined_score)], "data/stringdb/stringdb.score.tsv", quote = T, sep = "\t", col.names = T, row.names = F)
if(file.exists("data/stringdb/stringdb.score.tsv.gz")) {
  file.remove("data/stringdb/stringdb.score.tsv.gz")
}
system("gzip -9v data/stringdb/stringdb.score.tsv")