#!/usr/bin/env Rscript

library(data.table)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein where accession is not null")
dbDisconnect(conn)
rm(conn)
setDT(protein)

ensembl <- fread("data/uniprot/human.ensembl.tsv", header = T, sep = "\t")
ensembl <- rbindlist(list(ensembl, fread("data/uniprot/mouse.ensembl.tsv", header = T, sep = "\t", na.strings = "", quote = "")), use.names = T)
ensembl <- rbindlist(list(ensembl, fread("data/uniprot/rat.ensembl.tsv", header = T, sep = "\t", na.strings = "", quote = "")), use.names = T)


download.file("ftp://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.92.uniprot.tsv.gz", destfile = "data/ensembl/Homo_sapiens.GRCh38.92.uniprot.tsv.gz")
dt <- fread("gzcat data/ensembl/Homo_sapiens.GRCh38.92.uniprot.tsv.gz", header = T, sep = "\t", na.strings = "-", quote = "")
dt <- dt[, .(xref, gene_stable_id)]
setnames(dt, c("xref", "gene_stable_id"), c("ACCESSION", "ENSEMBL_GENE_ID"))
ensembl <- rbindlist(list(ensembl, dt), use.names = T)

download.file("ftp://ftp.ensembl.org/pub/current_tsv/mus_musculus/Mus_musculus.GRCm38.92.uniprot.tsv.gz", destfile = "data/ensembl/Mus_musculus.GRCm38.92.uniprot.tsv.gz")
dt <- fread("gzcat data/ensembl/Mus_musculus.GRCm38.92.uniprot.tsv.gz", header = T, sep = "\t", na.strings = "-", quote = "")
dt <- dt[, .(xref, gene_stable_id)]
setnames(dt, c("xref", "gene_stable_id"), c("ACCESSION", "ENSEMBL_GENE_ID"))
ensembl <- rbindlist(list(ensembl, dt), use.names = T)

download.file("ftp://ftp.ensembl.org/pub/current_tsv/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.92.uniprot.tsv.gz", destfile = "data/ensembl/Rattus_norvegicus.Rnor_6.0.92.uniprot.tsv.gz")
dt <- fread("gzcat data/ensembl/Rattus_norvegicus.Rnor_6.0.92.uniprot.tsv.gz", header = T, sep = "\t", na.strings = "-", quote = "")
dt <- dt[, .(xref, gene_stable_id)]
setnames(dt, c("xref", "gene_stable_id"), c("ACCESSION", "ENSEMBL_GENE_ID"))
ensembl <- rbindlist(list(ensembl, dt), use.names = T)

ensembl <- unique(ensembl)
ensembl <- merge(ensembl, protein, by.x = "ACCESSION", by.y = "accession")

fwrite(ensembl[, .(protein_id, ENSEMBL_GENE_ID)], file = "data/ensembl/mapping.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")