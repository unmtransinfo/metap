#!/usr/bin/env Rscript

library(RPostgreSQL)
library(data.table)

#change this to run it on chromium and connect remotely to METAP db

lincs.conn <- dbConnect(PostgreSQL(), user = "oleg", password = "qas2wer", host = "chromium", dbname = "lincs")
lincs <- dbGetQuery(lincs.conn, "SELECT level5_lm.pr_gene_id, level5_lm.zscore, perturbagen.drugcentral_id, perturbagen.is_parent, signature.cell_id FROM level5_lm,perturbagen,signature where level5_lm.sig_id=signature.sig_id and signature.pert_id = perturbagen.pert_id and drugcentral_id is NOT NULL")
dbDisconnect(lincs.conn)
rm(lincs.conn)
setDT(lincs)

lincs[, pert_id := ifelse(is.na(is_parent), sprintf("%d", drugcentral_id), sprintf("%d.P", drugcentral_id) )]
lincs[, `:=`(drugcentral_id = NULL, is_parent = NULL)]
lincs <- lincs[, .SD[which.max(abs(zscore))], by = .(pr_gene_id, cell_id, pert_id)]

metap.conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(metap.conn, "select protein_id,gene_id from ncbi")
dbDisconnect(metap.conn)
rm(metap.conn)
setDT(protein)

lincs <- merge(lincs, protein, by.x = "pr_gene_id", by.y = "gene_id", allow.cartesian = T)
lincs[, pr_gene_id := NULL]
lincs <- lincs[, .SD[which.max(abs(zscore))], by = .(protein_id, cell_id, pert_id)]

fwrite(lincs, file = "data/lincs/level5_max.tsv", col.names = T, row.names = F, na = "None", quote = T, sep = "\t")
if(file.exists("data/lincs/level5_max.tsv.gz")) {
	file.remove("data/lincs/level5_max.tsv.gz")
}
system("gzip -9v data/lincs/level5_max.tsv")