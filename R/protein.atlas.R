#!/usr/bin/env Rscript

library(data.table)

download.file("http://www.proteinatlas.org/download/normal_tissue.tsv.zip", destfile = "data/hpa/normal_tissue.tsv.zip")

norm.tissue <- fread("7z x -so data/hpa/normal_tissue.tsv.zip", header = T, sep = "\t", quote = "", na.strings = "")

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,ensembl_gene_id from ensembl")
dbDisconnect(conn)
rm(conn)
setDT(protein)

norm.tissue <- merge(norm.tissue, protein, by.x = "Gene", by.y = "ensembl_gene_id")
norm.tissue[, `:=`(Level = tolower(Level), Reliability = tolower(Reliability))]
norm.tissue[, Level := factor(x = Level, levels = c("not detected", "low", "medium", "high"), ordered = T)]
norm.tissue[, Reliability := factor(x = Reliability, levels = c("uncertain", "supported", "approved"), ordered = T)]
norm.tissue[, Tissue := sub("\\s\\d+$","", Tissue)]
norm.tissue <- norm.tissue[, head(.SD[order(-Reliability, -Level)], 1), by = .(protein_id, Tissue, `Cell type`)]

fwrite(norm.tissue[, .(protein_id, Tissue, `Cell type`, Level, Reliability)], file = "data/hpa/normal_tissue.tsv", sep = "\t", quote = T, col.names = T, row.names = F, na = "None")

if(file.exists("data/hpa/normal_tissue.tsv.gz")) {
  file.remove("data/hpa/normal_tissue.tsv.gz")
}
system("gzip -9v data/hpa/normal_tissue.tsv")