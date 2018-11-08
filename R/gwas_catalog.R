#!/usr/bin/env Rscript

library(data.table)
library(tidyr)
library(RPostgreSQL)

download.file(ur = "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative", destfile = "data/gwas/ebi_gwas_catalog.tsv")

gwas <- fread("data/gwas/ebi_gwas_catalog.tsv", header = T, sep = "\t", quote = "", na.strings = "")

gwas <- gwas[, .(PUBMEDID, STUDY, `DISEASE/TRAIT`, MAPPED_GENE, SNPS, CONTEXT, INTERGENIC, `P-VALUE`, PVALUE_MLOG, `OR or BETA`, CNV, MAPPED_TRAIT, MAPPED_TRAIT_URI)]
gwas <- separate_rows(gwas, "MAPPED_GENE", convert = T)
gwas <- gwas[!is.na(MAPPED_GENE)]
gwas[, `P-VALUE` := as.numeric(`P-VALUE`)]

conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)
protein <- protein[!is.na(symbol)]

gwas <- merge(gwas, protein, by.x = "MAPPED_GENE", by.y = "symbol", sort = F)
fwrite(gwas[, .(protein_id, PUBMEDID, STUDY, `DISEASE/TRAIT`, CONTEXT, INTERGENIC, `P-VALUE`, PVALUE_MLOG, `OR or BETA`, CNV, MAPPED_TRAIT, MAPPED_TRAIT_URI)], "data/gwas/ebi_gwas_catalog_imp.tsv", col.names = T, row.names = F, sep = "\t", na = "None", quote = T)
