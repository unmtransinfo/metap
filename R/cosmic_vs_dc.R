#!/usr/bin/env Rscript

library(data.table)

dc.antineo <- fread("data/drugcentral/antineoplastic.tsv", header = F, sep = "\t", quote = "\"", col.names = c("struct_id", "drug_name", "target_id", "component_id", "accession", "swissprot", "symbol", "ncbi_gene_id"))
gene2drug <- dc.antineo[, .(drugs = paste0(drug_name, collapse = "; ")), by = ncbi_gene_id]
dc.antineo[, drugs := paste0(drug_name, collapse = "; "), by = ncbi_gene_id]
cosmic <- fread("data/cosmic/cancer_gene_census.csv", header = T, sep = ",", quote = "\"")
setnames(cosmic, "Entrez GeneId", "ncbi_gene_id")
cosmic <- merge(cosmic, unique(dc.antineo[, .(ncbi_gene_id, accession)]), by.x = "ncbi_gene_id", by.y = "ncbi_gene_id", all.x = T, sort = F)
cosmic <- merge(cosmic, gene2drug, by.x = "ncbi_gene_id", by.y = "ncbi_gene_id", all.x = T, sort = F)
fwrite(cosmic[order(-drugs)], "data/cosmic/cancer_gene_census_drugtargets.tsv", col.names = T, row.names = F, sep = "\t", quote = T)
dc.antineo[, drug_name := NULL]
dc.antineo[, struct_id := NULL]
dc.antineo <- unique(dc.antineo, by = "ncbi_gene_id")
fwrite(dc.antineo[!ncbi_gene_id %in% cosmic$ncbi_gene_id], "data/drugcentral/antineoplastic_not_cosmic.tsv", col.names = T, row.names = F, sep = "\t", quote = T)

dt <- readRDS("data/input/104300.rds")

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select ncbi.protein_id,protein.symbol,protein.accession,gene_id,protein.name from protein,ncbi where ncbi.protein_id = protein.protein_id and tax_id=9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

dc.antineo <- merge(dc.antineo, protein[, .(protein_id, gene_id)], by.x = "ncbi_gene_id", by.y = "gene_id", all.x = T, sort = F)
cosmic <- merge(cosmic, protein[, .(protein_id, gene_id)], by.x = "ncbi_gene_id", by.y = "gene_id", sort = F)

pos.idx <- dc.antineo[, unique(protein_id)]
neg.idx <- cosmic[!protein_id %in% pos.idx, unique(protein_id)]

dt[, subset := "test"]
dt[, Y := "neg"]
dt[id1 %in% c(pos.idx, neg.idx), subset := "train"]
dt[id1 %in% pos.idx, Y := "pos"]

attr(dt, "omim_title") <- NULL
attr(dt, "mim_id") <- NULL
attr(dt, "Date") <- Sys.Date()
attr(dt, which = "title") <- "COSMIC oncogenes vs antineoplastics Tclin"

saveRDS(dt, "data/input/cosmic_vs_dc.rds")
