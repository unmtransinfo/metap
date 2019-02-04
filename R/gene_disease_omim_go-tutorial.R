#!/usr/bin/env Rscript
#############################################################################
### Simplified version of gene_disease_omim.R for tutorial purposes.
### Only OMIM phenotypic series considered for queries.
### Only ClinVar and UniProt disease-gene associations considered as known.
### Only GO metapaths considered for features.
### However, all key steps here are consistent with full method:
###	(1) Data from db pulled into R (memory) for performance.
#############################################################################
library(data.table)
library(RPostgreSQL)
library(Matrix)

###
#Calculates mode of a factor variable
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
###
#Replaces NA in a data.table with specified value
replace_na <- function(DT, ind, replacement) {
  for(j in ind){
    set(DT, i = which(is.na(DT[[j]])), j = j, value = replacement)
  }
}
###
# Scales numeric columns in a data.table using R base scale function
scale.data.table <- function(dt) {
  col.names <- colnames(dt)[2:ncol(dt)]
  dt[, (col.names) := lapply(.SD, scale), .SDcols=col.names]
}
#
conn <- dbConnect(PostgreSQL(), host = "localhost", dbname = "metap")
###
#Extract all OMIM phenotype series to OMIM mappings
omim.ps <- dbGetQuery(conn, "SELECT * FROM omim2phen_series")
#Extract all OMIM phenotype series terms
ps.name <- dbGetQuery(conn, "SELECT * FROM omim_phen_series");
###
# Extract all UniProt (are the same as OMIM) diseases
uniprot.disease <- dbGetQuery(conn, "SELECT * FROM uniprot_disease");
# Extract all disease/phenotype gene associations from ClinVar. This list is roughly half the clinical_significance
# values. Could there be a simpler rule?
clinvar <- dbGetQuery(conn, "SELECT DISTINCT clinvar.protein_id,clinvar.clinical_significance,clinvar_disease.cv_dis_id,phenotype,source,source_id FROM clinvar_disease,clinvar_disease_xref,clinvar WHERE clinvar_disease.cv_dis_id=clinvar_disease_xref.cv_dis_id AND source = 'OMIM' AND clinvar.cv_dis_id=clinvar_disease.cv_dis_id AND clinical_significance IN ('Pathogenic, Affects','Benign, protective, risk factor','Pathogenic/Likely pathogenic','Pathogenic/Likely pathogenic, other','Pathogenic, other','Affects','Pathogenic, other, protective','Conflicting interpretations of pathogenicity, Affects, association, other','Pathogenic/Likely pathogenic, drug response','Pathogenic, risk factor','risk factor','Pathogenic, association','Conflicting interpretations of pathogenicity, Affects, association, risk factor','Pathogenic/Likely pathogenic, risk factor','Affects, risk factor','Conflicting interpretations of pathogenicity, association, other, risk factor','Likely pathogenic, association','association, protective','Likely pathogenic, Affects','Pathogenic','Conflicting interpretations of pathogenicity, association','Pathogenic/Likely pathogenic, Affects, risk factor','Conflicting interpretations of pathogenicity, other, risk factor','association, risk factor','Benign, protective','Conflicting interpretations of pathogenicity, risk factor','Uncertain significance, protective','association','Uncertain significance, Affects','protective, risk factor','Pathogenic, association, protective','Pathogenic, protective','Likely pathogenic, other','Pathogenic, protective, risk factor','Benign, association, protective','Conflicting interpretations of pathogenicity, Affects','Benign/Likely benign, protective','protective')")

# Extract list of human proteins
protein <- dbGetQuery(conn, "SELECT * FROM protein WHERE tax_id = 9606");
dbDisconnect(conn)
rm(conn)

#Convert to data.table
setDT(omim.ps)
setDT(uniprot.disease)
setDT(ps.name)
setDT(clinvar)
setDT(protein)

# Map OMIM phenotypes in UniProt to phenotypic series id
# (by OMIM id, aka "phenotype MIM number")
uniprot.disease <- merge(uniprot.disease, omim.ps, by.x = "ref_id", by.y = "mim", all.x = T)

# Map OMIM phenotypic series id to names
uniprot.disease <- merge(uniprot.disease, ps.name, by.x = "omim_phen_series_id", by.y = "omim_phen_series_id", all.x = T)
# Change column names to match OMIM column names
setnames(uniprot.disease, c("ref_id", "term", "title"), c("mim", "mim_title", "pheno_series_title"))
uniprot.disease[, `:=`(disease_id = NULL, dbref = NULL)]

# Map OMIM phenotypes in ClinVar to phenotypic series id
#(by OMIM id, aka "phenotype MIM number")
clinvar$source_id <- as.integer(clinvar$source_id)
clinvar <- merge(clinvar, omim.ps, by.x = "source_id", by.y = "mim", all.x = T)
# Map OMIM phenotypic series id to names
clinvar <- merge(clinvar, ps.name, by.x = "omim_phen_series_id", by.y = "omim_phen_series_id", all.x = T)
# Change column names to match OMIM column names
setnames(clinvar, c("source_id","phenotype","title"), c("mim","mim_title","pheno_series_title"))
clinvar[, `:=`(clinical_significance = NULL, cv_dis_id = NULL, source = NULL)]
clinvar <- unique(clinvar)

# Concatenate gene - OMIM phenotype associations from ClinVar, UniProt, and RGD
g2d <- rbindlist(list(uniprot.disease, clinvar), use.names = T)
g2d <- unique(g2d)
writeLines(sprintf("Total genes associated with ANY PS: %d", length(unique(g2d$protein_id[!is.na(g2d$omim_phen_series_id)]))))
#
# Select OMIM phentypic series associated with sufficient genes (MIN_NGENES).
MIN_NGENES <- 50L
top.phentype.series <- g2d[!is.na(pheno_series_title), .(uq_prot = uniqueN(protein_id)), by = pheno_series_title][uq_prot >= MIN_NGENES]
#
if (nrow(top.phentype.series)==0) {
  message(sprintf("ERROR: No phenotypic series found with NGENES>=%d", MIN_NGENES))
  stop()
}
writeLines(sprintf("PS: %-44s; NGENES: %3d", top.phentype.series$pheno_series_title, top.phentype.series$uq_prot))

# Add OMIM phenotypic series ids
top.phentype.series <- merge(top.phentype.series, unique(g2d[, .(pheno_series_title, omim_phen_series_id)]), by.x = "pheno_series_title", by.y = "pheno_series_title")

#Extract data matrix for OMIM phenotypic series (random selection for this demo).
i <- sample(1:nrow(top.phentype.series), 1)
pheno_series <- top.phentype.series[i, pheno_series_title]
writeLines(sprintf("PS: %s (%s); NGENES: %3d", pheno_series,
	top.phentype.series$omim_phen_series_id[i], 
	top.phentype.series$uq_prot[i]))
###
# Separate genes (protein_ids), from those associated with ANY OMIM Phenotypic Series (PS) into 2 sets, based on given PS:
#	- Known ASSOCIATED with this PS
#	- Known NOT ASSOCIATED with this PS
# (based on available data, OMIM and more).
# Assign ASSOCIATED to RIGHT SIDE, NOT_ASSOCIATED to LEFT SIDE.
###
right.side <- g2d[pheno_series_title == pheno_series, unique(protein_id)]
left.side <- g2d[pheno_series_title != pheno_series & !protein_id %in% right.side, unique(protein_id)]
###
#All other human proteins are the external "test" set subject to ML model association predictions.
### 
test.idx <- protein[!protein_id %in% c(left.side, right.side), protein_id]
writeLines(sprintf("KNOWN ASSOCIATED: %d; KNOWN NOT ASSOCIATED: %d; NOT KNOWN: %d", length(right.side), length(left.side), length(test.idx)))
### 
#Initialize data matrix for ML with protein_ids, subset (train/test), labels (pos/neg).
### 
dt <- data.table(id1 = c(right.side, left.side, test.idx),
	Y = c(rep_len("pos", length.out = length(right.side)), rep_len("neg", length.out = length(left.side)), rep_len("neg", length.out = length(test.idx))),
	subset = c(rep_len("train", length.out = length(right.side)), rep_len("train", length.out = length(left.side)), rep_len("test", length.out = length(test.idx))))
dt$Y <- factor(dt$Y, levels = c("neg", "pos"))
dt$subset <- factor(dt$subset, levels = c("train", "test"))
###
#All proteins not associated with the disease/phenotype are on the left side
left.side <- c(left.side, test.idx)
###
# Extract Gene Ontology annotations: find all genes sharing at least one annotation with 
# our KNOWN ASSOCIATED set (right side).
conn <- dbConnect(PostgreSQL(), host = "localhost", dbname = "metap")
go <- dbGetQuery(conn, sprintf("SELECT DISTINCT protein_id,go_id FROM goa WHERE go_id IN (SELECT DISTINCT go_id FROM goa WHERE protein_id IN (%s))", paste(right.side, collapse = ",")))
dbDisconnect(conn)
rm(conn)
###
setDT(go)
###
# Move all protein ids in GO to right/left side based on right/left sets defined above.
# Each row is one metapath instance: gene->gene->phenotypic_series.
###
go <- go[protein_id %in% c(left.side, right.side)]
go.right.side <- go[protein_id %in% right.side]
go.left.side <- go[protein_id %in% c(right.side, left.side)]
setnames(go.left.side, "protein_id", "id1")
setnames(go.right.side, "protein_id", "id2")
go <- merge(go.left.side, go.right.side, by.x = "go_id", by.y = "go_id", allow.cartesian = T)
go <- go[id1 != id2]
###

###
# Compute path degree weighted product (PDP). Generally de-weight associations according
# to commonality. PDP used by ML as feature representing the strength of metapath instance.
###
disease.conn.count <- go[, uniqueN(id2)] # Unique proteins on the right side (associated)
go[, id2.count := .N, by = id2] #Count id1's (unassociated) for each id2 (associated)
go[, id1.count := .N, by = id1] #Count id2's (associated) for each id1 (unassociated)
go[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
###
# Transform from long to wide format; features as columns.
go <- dcast(go, id1 ~ go_id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
# Scale the data
scale.data.table(go)
# Merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, go, by.x = "id1", by.y = "id1", all.x = T)
rm(go)
###
# Replace NA in numeric columns with 0.
# (Oleg's suggestion: Future iteration use multiple imputations package mice.)
col.ids <- which(sapply(dt, is.numeric))
col.ids <- col.ids[2:length(col.ids)] #Skip id1
replace_na(dt, col.ids, 0)
#
#Store metadata as attributes.
setattr(dt, "pheno_series_title", pheno_series)
setattr(dt, "omim_phen_series_id", top.phentype.series$omim_phen_series_id[i])
setattr(dt, "Date", Sys.Date())
#
print(table(dt$subset, dt$Y))
#
fpath <- sprintf("data/input/%s.rds", top.phentype.series$omim_phen_series_id[i])
saveRDS(dt, fpath)
writeLines(sprintf("Saved ML-ready dataset for \"%s\" to: %s", pheno_series, fpath))
#
