#!/usr/bin/env Rscript

library(data.table)
library(tidyr)

download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz", destfile = "data/clinvar/variant_summary.txt.gz")
vs <- fread("gzcat data/clinvar/variant_summary.txt.gz", header = T, sep = "\t", quote = "", na.strings = "-")
vs <- vs[!is.na(PhenotypeList)]
vs <- vs[!is.na(GeneID) & GeneID != -1]
vs <- separate_rows(vs, PhenotypeIDS, PhenotypeList, sep = ";")

disease <- unique(vs[, .(PhenotypeIDS, PhenotypeList)])
disease[PhenotypeIDS == "na", PhenotypeIDS := NA]
disease[PhenotypeIDS == "", PhenotypeIDS := NA]
disease <- separate_rows(disease, PhenotypeIDS, sep = ",")
disease <- disease[, c("Source", "ID", "X") := tstrsplit(PhenotypeIDS, ":", fixed = T, keep = c(1,2,3))]
disease[!is.na(X), ID := sprintf("%s:%s", ID, X)]
disease[, X:= NULL]
dis.seq <- disease[, .(PhenoType = unique(PhenotypeList))]
dis.seq[, SeqID := .I]
disease <- merge(disease, dis.seq, by.x = "PhenotypeList", by.y = "PhenoType")
fwrite(unique(disease[, .(SeqID, PhenotypeList)]), file = "data/clinvar/cv_disease.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")
fwrite(unique(disease[, .(SeqID, Source, ID)]), file = "data/clinvar/cv_disease_xref.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")

vs <- merge(vs, dis.seq, by.x = "PhenotypeList", by.y = "PhenoType")
vs[, `:=`(PhenotypeList = NULL, PhenotypeIDS = NULL, OtherIDs = NULL, Cytogenetic = NULL, AlternateAllele = NULL, ReferenceAllele = NULL, Assembly = NULL, ChromosomeAccession = NULL, Chromosome = NULL, Start = NULL, Stop = NULL, RCVaccession = NULL, HGNC_ID = NULL, Guidelines = NULL)]
setnames(vs, c("RS# (dbSNP)", "nsv/esv (dbVar)", "#AlleleID"), c("dbSNP_rs", "dbVar_id", "AlleleID"))
vs[dbSNP_rs == -1, dbSNP_rs := NA]
vs[, TestedInGTR := ifelse(TestedInGTR == "Y", T, F)]

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,gene_id from ncbi")
dbDisconnect(conn)
rm(conn)
setDT(protein)

vs <- merge(vs, protein, by.x = "GeneID", by.y = "gene_id")

fwrite(vs[, .(AlleleID, Type, Name, ClinicalSignificance, ClinSigSimple, LastEvaluated, dbSNP_rs, dbVar_id, Origin, OriginSimple, ReviewStatus, NumberSubmitters, TestedInGTR, SubmitterCategories, SeqID, protein_id)], file = "data/clinvar/clinvar.tsv", col.names = T, row.names = F, quote = T, sep = "\t", na = "None")

if(file.exists("data/clinvar/clinvar.tsv.gz")) {
  file.remove("data/clinvar/clinvar.tsv.gz")
}
system("gzip -9v data/clinvar/clinvar.tsv")