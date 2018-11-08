#!/usr/bin/env Rscript

library(readr)
library(data.table)
library(RPostgreSQL)

conn <- dbConnect(PostgreSQL(), dbname = "metap_dev")
protein <- dbGetQuery(conn, "select * from protein2mgi")
dbDisconnect(conn)
rm(conn)
setDT(protein)
cat(protein[, uniqueN(mgi_id)],"unique MGI_IDs extracted from METAP\n")

readUrl <- function(url) {
  out <- tryCatch(
    {
      read_csv(url, na = c("", "NA", "null", "None", "NULL", "na", "N/A"))
      #fread(url, header = T, sep = ",", quote = "\"") 
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", url))
      # Choose a return value in case of error
      return(NULL)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", url))
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}

write.dt <- function(dt, curr.id = NA, fn = "data/impc/genotype-phenotype.tsv") {
  tryCatch(
    {
      fwrite(dt, file = fn, append = file.exists(fn), col.names = !file.exists(fn), row.names = F, quote = T, sep = "\t", na = "None")
    },
    error = function(cond) {
      message("Error writing data for: ", curr.id)
    },
    warning = function(cond) {
      message("Warning writing data for: ", curr.id)
    }
  )
}

protein <- protein[order(protein_id)]
max.prot.id <- NA
if(file.exists("data/impc/curr.protid.txt")) {
  max.prot.id <- scan(file = "data/impc/curr.protid.txt", what = integer())
  cat("Loading saved progress from previous run, last processed protein_id:", max.prot.id,"\n")
}

for(i in 1:nrow(protein)) {
  curr.prot.id <- protein[i, protein_id]
  if(!is.na(max.prot.id) & max.prot.id <= curr.prot.id) {
    cat("skipping current protein_id:",curr.prot.id,"already processed in previous run\n")
    next
  }
  mgi.id <- protein[i, mgi_id]
  url <- sprintf("http://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=marker_accession_id:\"%s\"&wt=csv&rows=1000000", mgi.id)
  data <- readUrl(url)
  cat(curr.prot.id, file = "data/impc/curr.protid.txt")
  if(is.null(data)) {
    next
  }
  if(nrow(data) == 0) {
    next
  }
  setDT(data)
  data <- data[, .(mp_term_id, mp_term_name, marker_symbol, marker_accession_id, p_value, effect_size, procedure_name, parameter_name)]
  data <- unique(data, by = c("mp_term_id", "marker_accession_id"))
  data[, association := T]
  data <- data[!is.na(mp_term_id)]
  write.dt(data, curr.id = mgi.id)
  cat("processed MGI:",mgi.id,"\n")
}

geno2pheno <- fread("data/impc/genotype-phenotype.tsv", header = T, sep = "\t", na.strings = "None", quote = "\"")


mpterm.list <- geno2pheno[, unique(mp_term_id)]
mpterm.list <- sort(mpterm.list)
max.mpterm.id <- NA
if(file.exists("data/impc/curr.mptermid.txt")) {
  max.mpterm.id <- scan(file = "data/impc/curr.mptermid.txt", what = character())
  cat("Loading saved progress from previous run, last processed mp_term_id:", max.mpterm.id,"\n")
}
for(curr.mpterm.id in mpterm.list) {
  if(!is.na(max.mpterm.id) & max.mpterm.id <= curr.mpterm.id) {
    cat("skipping current mp_term:",curr.mpterm.id,"already processed in previous run\n")
    next
  }
  url <- sprintf("https://www.ebi.ac.uk/mi/impc/solr/statistical-result/select?q=mp_term_id:\"%s\"&wt=csv&rows=1000000", curr.mpterm.id)
  data <- readUrl(url)
  cat(curr.mpterm.id, file = "data/impc/curr.mptermid.txt")
  if(is.null(data)) {
    next
  }
  if(nrow(data) == 0) {
    next
  }
  setDT(data)
  mgi.list <- geno2pheno[mp_term_id == curr.mpterm.id, unique(marker_accession_id)]
  data <- data[, .(mp_term_id, mp_term_name, marker_symbol, marker_accession_id, p_value, effect_size, procedure_name, parameter_name)]
  data <- data[!marker_accession_id %chin% mgi.list]
  data <- unique(data, by = c("mp_term_id", "marker_accession_id"))
  data[, association := F]
  data <- data[!is.na(mp_term_id)]
  write.dt(data, curr.id = curr.mpterm.id)
  cat("processed MP term:",curr.mpterm.id,"\n")
}

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,mgi_id from protein2mgi")
mpo <- dbGetQuery(conn, "select mp_term_id from mp_onto")
dbDisconnect(conn)
rm(conn)
setDT(protein)
setDT(mpo)

geno2pheno <- fread("data/impc/genotype-phenotype.tsv", header = T, sep = "\t", na.strings = "None", quote = "\"")
geno2pheno <- merge(geno2pheno, protein, by.x = "marker_accession_id", by.y = "mgi_id", allow.cartesian = T)
geno2pheno <- geno2pheno[, .(protein_id, mp_term_id, p_value, effect_size, procedure_name, parameter_name, association)]
geno2pheno <- unique(geno2pheno, by = c("protein_id", "mp_term_id"))
geno2pheno[, mp_term_id := sub(":", "_", mp_term_id, fixed = T)]
geno2pheno <- geno2pheno[mp_term_id %chin% mpo$mp_term_id]
fwrite(geno2pheno, file = "data/impc/genotype-phenotype.tsv", col.names = T, row.names = F, sep = "\t", na = "None", quote = T)
