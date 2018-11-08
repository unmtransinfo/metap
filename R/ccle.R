#!/usr/bin/env Rscript

library(data.table)

dt <- fread("gzcat data/CCLE/CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct.gz", header = T, skip = 2, sep = "\t", quote = "")
dt[, Description := NULL]

dt <- melt(dt, id.vars = "Name", variable.name = "cell_id", value.name = "expression", variable.factor = F, value.factor = F)
dt[, cell_id2 := tstrsplit(cell_id, "_", fixed = T, keep = 1)]
dt[cell_id2 != cell_id, tissue := substr(cell_id, nchar(cell_id2) + 2, nchar(cell_id))]
dt[, cell_id := NULL]
setnames(dt, c("cell_id2"), c("cell_id"))
dt[, ENSG := tstrsplit(Name, ".", fixed = T, keep = 1)]
dt[, Name := NULL]



library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein.protein_id,ensembl_gene_id from ensembl,protein where ensembl.protein_id = protein.protein_id and protein.tax_id=9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

dt <- merge(dt, protein, by.x = "ENSG", by.y = "ensembl_gene_id", all.x = T, allow.cartesian = T)
dt <- dt[!is.na(protein_id)]
dt[, ENSG := NULL]

setcolorder(dt, c("protein_id", "cell_id", "tissue", "expression"))
dt <- unique(dt)
fwrite(dt, file = "data/CCLE/ccle.tsv", quote = T, sep = "\t", col.names = T, row.names = F, na = "None")
if(file.exists("data/CCLE/ccle.tsv.gz")) {
  file.remove("data/CCLE/ccle.tsv.gz")
}
system("gzip -9v data/CCLE/ccle.tsv")
