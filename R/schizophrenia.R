#!/usr/bin/env Rscript

library(readxl)
library(data.table)
library(tidyr)
library(RPostgreSQL)

set.seed(1234)

nature <- read_xlsx("data/gwas/nature13595_SI3.xlsx", sheet = 1)
nature <- separate_rows(nature, "Protein coding genes", sep = " ", convert = T)
setDT(nature)
nature <- nature[!is.na(`Protein coding genes`)]

conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession,symbol from protein where accession is not null and tax_id = 9606")
gwas.schizophrenia <- dbGetQuery(conn, "select DISTINCT protein_id from gwas WHERE mapped_trait is not null and (mapped_trait ilike '%schizophrenia%'
                                                                           or mapped_trait ilike '%bipolar%'
                                                                           or mapped_trait ilike '%psychosis%'
                                                                           or mapped_trait ilike '%psychotic%'
                                                                            or mapped_trait ilike '%schizoaffective%'
                                                                            or mapped_trait ilike '%behavioral disorder%'
                                                                            or mapped_trait ilike '%anxiety%'
                                                                            or mapped_trait ilike '%unipolar%'
                                                                            or mapped_trait ilike '%depressive%'
                                                                            or mapped_trait ilike '%depression%'
                                                                            or mapped_trait ilike '%attention deficit%'
                                                                            or mapped_trait ilike '%hyperactivity disorder%'
                                                                            or mapped_trait ilike '%dementia%'
                                                                            or mapped_trait ilike '%tourette%'
                                                                            or mapped_trait ilike '%alzheimers%'
                                                                            or mapped_trait ilike '%autism%'
                                                                            or mapped_trait ilike '%obsessive%'
                                                                            or mapped_trait ilike '%compulsive%'
                                                                            or mapped_trait ilike '%parkinson%'
                                                                            or mapped_trait ilike '%personality disorder%')")
gwas.all <- dbGetQuery(conn, "select DISTINCT protein_id from gwas where mapped_trait is not null")
dbDisconnect(conn)
rm(conn)
setDT(protein)
protein <- protein[!is.na(symbol)]
setDT(gwas.schizophrenia)
setDT(gwas.all)

nature <- merge(nature, protein, by.x = "Protein coding genes", by.y = "symbol", sort = F)
nature[, `:=`(OMIM = NULL, "NHGRI GWAS catalog" = NULL, "KO phenotype" = NULL)]
nature <- unique(nature, by = "protein_id")
gwas.all <- gwas.all[!protein_id %in% gwas.schizophrenia$protein_id]
gwas.all <- gwas.all[!protein_id %in% nature$protein_id]

dt <- readRDS("data/input/104300.rds")

dt[, subset := "test"]
dt[id1 %in% nature$protein_id, subset := "train"]
dt[id1 %in% gwas.all$protein_id, subset := "train"]
dt[id1 %in% nature$protein_id, Y := "pos"]
dt[id1 %in% gwas.all$protein_id, Y := "neg"]

attr(dt, "omim_title") <- NULL
attr(dt, "mim_id") <- NULL
attr(dt, "Date") <- Sys.Date()
attr(dt, which = "title") <- "Schizophrenia"
attr(dt, which = "DOI") <- "10.1038/nature13595"

saveRDS(dt, "data/input/nature_schizophrenia.rds")
