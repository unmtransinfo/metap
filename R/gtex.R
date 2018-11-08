#!/usr/bin/env Rscript

library(data.table)

rnaseq <- fread('gunzip -c data/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz', header = T, skip = 2, sep = "\t", quote = "")
rnaseq[, c("ENSEMBL_GENEID") := tstrsplit(Name, ".", fixed = T, keep = 1)]
rnaseq[, `:=`(Name = NULL, Description = NULL)]
rnaseq <- melt(rnaseq, id.vars = "ENSEMBL_GENEID", variable.name = "SAMPID", value.name = "RPKM", na.rm = T, value.factor = F, verbose = T)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), dbname = "metap_dev")
protein <- dbGetQuery(conn, "select protein_id,ensembl_gene_id from ensembl")
dbDisconnect(conn)
rm(conn)
setDT(protein)

rnaseq <- merge(rnaseq, protein, by.x = "ENSEMBL_GENEID", by.y = "ensembl_gene_id")

sample.label <- fread("data/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt", header = T, sep = "\t", quote = "", na.strings = "")
sample.label <- sample.label[, .(SAMPID, SMATSSCR, SMTS, SMTSD)][!is.na(SMTSD)][, c("C1","C2") := tstrsplit(SAMPID, "-", fixed = T, keep = c(1,2))][, SUBJID := sprintf("%s-%s", C1, C2)][, `:=`(C1 = NULL, C2 = NULL)]
subject.label <- fread("data/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt", header = T, sep = "\t", quote = "", na.strings = "")
subject.label <- subject.label[DTHHRDY == 1 | DTHHRDY == 2] #subjects healthy @ death
sample.label <- merge(sample.label, subject.label, by.x = "SUBJID", by.y = "SUBJID")
sample.label[, GENDER := factor(SEX, levels = c(1,2), labels = c("male", "female"))]
sample.label[, SEX := NULL]
sample.label[is.na(SMTS) & SMTSD %like% "Skin - ", SMTS := "Skin"]
sample.label <- sample.label[SMATSSCR < 2 | is.na(SMATSSCR)] # remove samples with high degree of autolysis

rnaseq <- merge(rnaseq, sample.label, by.x = "SAMPID", by.y = "SAMPID")

all <- rnaseq[, .(MEDIAN_RPKM = median(RPKM, na.rm = T)), by = .(protein_id, SMTS, SMTSD)]
all[, GENDER := NA]
all <- all[!is.na(MEDIAN_RPKM)]
gender <- rnaseq[, .(MEDIAN_RPKM = median(RPKM, na.rm = T)), by = .(protein_id, SMTS, SMTSD, GENDER)]
gender <- gender[!is.na(MEDIAN_RPKM)]

fwrite(all, file = "data/gtex/all.tsv", col.names = T, row.names = F, sep = "\t", na = "None", quote = T)
if(file.exists("data/gtex/all.tsv.gz")) {
	file.remove("data/gtex/all.tsv.gz")
}
system("gzip -9v data/gtex/all.tsv")
fwrite(gender, file = "data/gtex/gender.tsv", col.names = T, row.names = F, sep = "\t", na = "None", quote = T)
if(file.exists("data/gtex/gender.tsv.gz")) {
  file.remove("data/gtex/gender.tsv.gz")
}
system("gzip -9v data/gtex/gender.tsv")