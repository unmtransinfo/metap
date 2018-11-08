#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

targets <- read_delim("http://www.guidetopharmacology.org/DATA/targets_and_families.csv", col_names = T, delim = ",", guess_max = Inf)
setDT(targets)
targets <- targets[!is.na(`Human SwissProt`), .(`Human SwissProt`, Type, `Family name`, `Family id`)]
targets <- separate_rows(targets, "Human SwissProt", sep = "\\|", convert = T)

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein where protein.tax_id=9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)

targets <- merge(targets, protein, by.x = "Human SwissProt", by.y = "accession")

targets <- unique(targets, by = "protein_id")

fwrite(targets[, .(protein_id, Type, `Family name`, `Family id`)], file = "data/iuphar/target_class.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")