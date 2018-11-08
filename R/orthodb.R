#!/usr/bin/env Rscript

library(data.table)

genes <- fread("gzcat data/orthodb/odb9v1_genes.tab.gz", header = F, sep = "\t", col.names = c("Ortho_DBid","tax_id","protein_sequence_id","Uniprot","ENSEMBL_gene_name","NCBI_gid","description"), quote = "", na.strings = "\\N")
genes <- genes[tax_id %in% c(9606,10116,10090)]
OG2genes <- fread("gzcat data/orthodb/odb9v1_OG2genes.tab.gz", header = F, sep = "\t", col.names = c("OG_id","Ortho_DBid"), na.strings = "\\N", quote = "")

genes <- merge(genes, OG2genes, by.x = "Ortho_DBid", by.y = "Ortho_DBid")

genes[, `:=`(Ortho_DBid = NULL, tax_id = NULL, protein_sequence_id = NULL, ENSEMBL_gene_name = NULL, NCBI_gid = NULL, description = NULL)]
genes <- genes[!is.na(Uniprot)]

conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein")
dbDisconnect(conn)
rm(conn)
setDT(protein)

genes <- merge(genes, protein, by.x = "Uniprot", by.y = "accession", all.x = T)
