#!/usr/bin/env Rscript

library(data.table)

pathways <- fread("http://reactome.org/download/current/ReactomePathways.txt", header = F, sep = "\t", col.names = c("reactome_id", "pathway_name", "organism"))
fwrite(pathways, file = "data/reactome/pathways.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")
annon <- fread("http://reactome.org/download/current/UniProt2Reactome_All_Levels.txt", header = F, sep = "\t", col.names = c("accession","reactome_id", "url", "name", "evidence_code", "species"))

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein")
dbDisconnect(conn)
rm(conn)
setDT(protein)
annon <- merge(annon, protein, by.x = "accession", by.y = "accession", all.x = F, all.y = F)
annon <- annon[reactome_id %chin% pathways$reactome_id]
fwrite(annon[, .(protein_id, reactome_id, evidence_code)], file = "data/reactome/reactome_annon.tsv", col.names = T, row.names = F, sep = "\t", na = "None", quote = T)
