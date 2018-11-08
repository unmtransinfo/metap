#!/usr/bin/env Rscript

library(data.table)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,stringdb_id from stringdb")
dbDisconnect(conn)
rm(conn)
setDT(protein)
protein <- protein[, c("ensp") := tstrsplit(stringdb_id, ".", fixed = T, keep = 2)]
protein[, stringdb_id := NULL]

disease <- fread("http://download.jensenlab.org/human_disease_textmining_filtered.tsv", header = T, colClasses = c("character","character","character", "character", "numeric","numeric", "character"), col.names = c("id","symbol","doid","do_name","z_score","stars","url"), quote = "")
disease <- merge(disease, protein, by.x = "id", by.y = "ensp")

fwrite(disease[, .(protein_id, doid, do_name, z_score, stars, url)], "data/jensen/text_mining.tsv", quote = T, sep = "\t", row.names = F, col.names = T)

if(file.exists("data/jensen/text_mining.tsv.gz")) {
  file.remove("data/jensen/text_mining.tsv.gz")
}
system("gzip -9v data/jensen/text_mining.tsv")