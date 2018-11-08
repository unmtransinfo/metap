#!/usr/bin/env Rscript

library(data.table)
library(RPostgreSQL)
library(Matrix)
library(readxl)

#calculates mode of a factor variable
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#replaces NA in a data.table with specified value
replace_na <- function(DT, ind, replacement) {
  for(j in ind){
    set(DT, i = which(is.na(DT[[j]])), j = j, value = replacement)
  }
}

# scales numeric columns in a data.table using R base scale function
scale.data.table <- function(dt) {
  col.names <- colnames(dt)[2:ncol(dt)]
  dt[, (col.names) := lapply(.SD, scale), .SDcols=col.names]
}

args <- commandArgs(T)

#connects to Metapath database instance on localhost with "oleg" user, replace with actual username
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
protein <- dbGetQuery(conn, "select * from protein where tax_id = 9606");
uni.dis <- dbGetQuery(conn, "select distinct protein_id from uniprot_disease")
dbDisconnect(conn)
rm(conn)
setDT(protein)
setDT(uni.dis)

g2d <- read_xlsx(args[1], sheet = 1)
setDT(g2d)
phenotype <- read_xlsx(args[1], sheet = 2)
setDT(phenotype)
g2d <- unique(g2d, by = "Name")
g2d[, pheno_series_title := phenotype[1, disease_name]]
g2d[, omim_phen_series_id := phenotype[1, disease_id]]

g2d <- merge(g2d, protein[, .(protein_id, accession, symbol)], by.x = "Name", by.y = "symbol")

i <- 1
pheno_series <- g2d[i, pheno_series_title]
#separate proteins into 2 sets one known to be associated with OMIM phenotype based on available data and the other one not known
#assign to right side protein_ids known to be associated with OMIM phenotype
right.side <- g2d[pheno_series_title == pheno_series, unique(protein_id)]
#assign to left side protein_ids known to be associated with OMIM phenotype
left.side <- uni.dis[!protein_id %in% right.side, unique(protein_id)]
#all the other protein in human proteome are the external set and will apply XGBoost models to predict potential associations
test.idx <- protein[!protein_id %in% c(left.side, right.side), protein_id]

#initialize the data matrix for machine learning with protein_ids, subset annotations (train/test), flag positive/negative associations
dt <- data.table(id1 = c(right.side, left.side, test.idx), Y = c(rep_len("pos", length.out = length(right.side)), rep_len("neg", length.out = length(left.side)), rep_len("neg", length.out = length(test.idx))), subset = c(rep_len("train", length.out = length(right.side)), rep_len("train", length.out = length(left.side)), rep_len("test", length.out = length(test.idx))))
dt$Y <- factor(dt$Y, levels = c("neg", "pos"))
dt$subset <- factor(dt$subset, levels = c("train", "test"))

#all proteins not associated with the disease/phenotype are on the left side
left.side <- c(left.side, test.idx)

#extract all PPI scores from StringDB for proteins on the right side
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
stringdb <- dbGetQuery(conn, sprintf("select * from stringdb_score where protein_id1 in (%s) or protein_id2 in (%s)", paste(right.side, collapse = ","), paste(right.side, collapse = ",")))
dbDisconnect(conn)
rm(conn)
setDT(stringdb)

#move all protein ids in PPI pairs to right/left side based on right/left sets defined above
stringdb[protein_id2 %in% right.side, id2 := protein_id2]
stringdb[is.na(id2), `:=`(id2 = protein_id1, id1 = protein_id2)]
stringdb[is.na(id1), id1 := protein_id1]
stringdb[, `:=`(protein_id1 = NULL, protein_id2 = NULL)]
stringdb <- stringdb[id1 %in% c(right.side, left.side)]
stringdb[, id2.count := .N, by = id2]
stringdb[, id1.count := .N, by = id1]
#count unique proteins on the right side (associated with phenotype)
disease.conn.count <- stringdb[, uniqueN(id2)]
# compute path degree weighted product weighted by PPI score (combined_score)
stringdb[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)*combined_score]
# assign name to metapth
stringdb[, col.id := sprintf("pp.%d", id2)]
# transform from long to wide format
stringdb <- dcast(stringdb, id1 ~ col.id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
# scale the data
scale.data.table(stringdb)
# merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, stringdb, by.x = "id1", by.y = "id1", all.x = T)
rm(stringdb)

#extract Reactome pathways annotations
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
reactome <- dbGetQuery(conn, sprintf("select * from reactomea where reactome_id in (select distinct reactome_id from reactomea where protein_id in (%s))", paste(right.side, collapse = ",")))
dbDisconnect(conn)
rm(conn)
setDT(reactome)

#move all protein ids in Reactome pathways to right/left side based on right/left sets defined above
reactome[, `:=`(evidence = NULL)]
reactome <- unique(reactome)
reactome <- reactome[protein_id %in% c(left.side, right.side)]
reactome.right.side <- reactome[protein_id %in% right.side]
reactome.left.side <- reactome[protein_id %in% c(right.side, left.side)]
setnames(reactome.left.side, "protein_id", "id1")
setnames(reactome.right.side, "protein_id", "id2")
reactome <- merge(reactome.left.side, reactome.right.side, by.x = "reactome_id", by.y = "reactome_id", allow.cartesian = T)
reactome <- reactome[id1 != id2]
#count unique proteins on the right side (associated with phenotype)
disease.conn.count <- reactome[, uniqueN(id2)]
reactome[, id2.count := .N, by = id2]
reactome[, id1.count := .N, by = id1]
#compute path degree weighted product
reactome[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
# transform from long to wide format
reactome <- dcast(reactome, id1 ~ reactome_id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
# scale the data
scale.data.table(reactome)
# merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, reactome, by.x = "id1", by.y = "id1", all.x = T)
rm(reactome)

#extract KEGG pathways annotations
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
kegg <- dbGetQuery(conn, sprintf("select protein_id,kegg_pathway_id from kegg_pathway where kegg_pathway_id in (select distinct kegg_pathway_id from kegg_pathway where protein_id in (%s))", paste(right.side, collapse = ",")))
dbDisconnect(conn)
rm(conn)
setDT(kegg)

#move all protein ids in KEGG pathways to right/left side based on right/left sets defined above
kegg <- kegg[protein_id %in% c(left.side, right.side)]
kegg.right.side <- kegg[protein_id %in% right.side]
kegg.left.side <- kegg[protein_id %in% c(right.side, left.side)]
setnames(kegg.left.side, "protein_id", "id1")
setnames(kegg.right.side, "protein_id", "id2")
kegg <- merge(kegg.left.side, kegg.right.side, by.x = "kegg_pathway_id", by.y = "kegg_pathway_id", allow.cartesian = T)
kegg <- kegg[id1 != id2]
#count unique proteins on the right side (associated with phenotype)
disease.conn.count <- kegg[, uniqueN(id2)]
kegg[, id2.count := .N, by = id2]
kegg[, id1.count := .N, by = id1]
#compute path degree weighted product
kegg[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
# transform from long to wide format
kegg <- dcast(kegg, id1 ~ kegg_pathway_id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
# scale the data
scale.data.table(kegg)
# merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, kegg, by.x = "id1", by.y = "id1", all.x = T)
rm(kegg)

#Extract GTEx expression data
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
gtex <- dbGetQuery(conn, "select protein_id,median_tpm,tissue_type_detail from gtex")
dbDisconnect(conn)
rm(conn)
setDT(gtex)

# transform from long to wide format
gtex <- dcast(gtex, protein_id ~ tissue_type_detail, fun.aggregate = median, value.var = "median_tpm", drop = T, fill = 0)
#scale the data
scale.data.table(gtex)
# merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, gtex, by.x = "id1", by.y = "protein_id", all.x = T)
rm(gtex)

#Extract data from L1000 drug gene perturbation profiles
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
lincs <- dbGetQuery(conn, "select protein_id,pert_id||':'||cell_id as col_id,zscore from lincs")
dbDisconnect(conn)
rm(conn)
setDT(lincs)

# transform from long to wide format
lincs <- dcast(lincs, protein_id ~ col_id, fun.aggregate = median, value.var = "zscore", drop = T, fill = 0)
#scale the data
scale.data.table(lincs)
# merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, lincs, by.x = "id1", by.y = "protein_id", all.x = T)
rm(lincs)

#Extract Gene Ontology annotations
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
go <- dbGetQuery(conn, sprintf("select distinct protein_id,go_id from goa where go_id in (select distinct go_id from goa where protein_id in (%s))", paste(right.side, collapse = ",")))
dbDisconnect(conn)
rm(conn)
setDT(go)
#move all protein ids in GO to right/left side based on right/left sets defined above
go <- go[protein_id %in% c(left.side, right.side)]
go.right.side <- go[protein_id %in% right.side]
go.left.side <- go[protein_id %in% c(right.side, left.side)]
setnames(go.left.side, "protein_id", "id1")
setnames(go.right.side, "protein_id", "id2")
go <- merge(go.left.side, go.right.side, by.x = "go_id", by.y = "go_id", allow.cartesian = T)
go <- go[id1 != id2]
#count unique proteins on the right side (associated with phenotype)
disease.conn.count <- go[, uniqueN(id2)]
go[, id2.count := .N, by = id2]
go[, id1.count := .N, by = id1]
#compute path degree weighted product
go[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
# transform from long to wide format
go <- dcast(go, id1 ~ go_id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
# scale the data
scale.data.table(go)
# merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, go, by.x = "id1", by.y = "id1", all.x = T)
rm(go)

#Extract InterPro annotations
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
interpro <- dbGetQuery(conn, sprintf("select distinct protein_id,entry_ac from interproa where entry_ac in (select distinct entry_ac from interproa where protein_id in (%s))", paste(right.side, collapse = ",")))
dbDisconnect(conn)
rm(conn)
setDT(interpro)
#move all protein ids in InterPro to right/left side based on right/left sets defined above
interpro <- interpro[protein_id %in% c(left.side, right.side)]
interpro.right.side <- interpro[protein_id %in% right.side]
interpro.left.side <- interpro[protein_id %in% c(right.side, left.side)]
setnames(interpro.left.side, "protein_id", "id1")
setnames(interpro.right.side, "protein_id", "id2")
interpro <- merge(interpro.left.side, interpro.right.side, by.x = "entry_ac", by.y = "entry_ac", allow.cartesian = T)
interpro <- interpro[id1 != id2]
#count unique proteins on the right side (associated with phenotype)
disease.conn.count <- interpro[, uniqueN(id2)]
interpro[, id2.count := .N, by = id2]
interpro[, id1.count := .N, by = id1]
#compute path degree weighted product
interpro[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
# transform from long to wide format
interpro <- dcast(interpro, id1 ~ entry_ac, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
#scale the data
scale.data.table(interpro)
# merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, interpro, by.x = "id1", by.y = "id1", all.x = T)
rm(interpro)

#Extract Protein Atlas normal tissue expression data
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
hpa <- dbGetQuery(conn, "select protein_id, tissue||'.'||cell_type as col_id,level from hpa_norm_tissue where reliability in ('supported','approved')")
dbDisconnect(conn)
rm(conn)
setDT(hpa)

#convert to factor expression levels
hpa$level <- factor(hpa$level, levels = c("not detected", "low", "medium", "high"), ordered = F)
hpa <- unique(hpa)
#normalize cell line names
hpa[, col_id := gsub(" ", "_", col_id, fixed = T)]
hpa[, col_id := gsub("/", "_", col_id, fixed = T)]
hpa[, col_id := gsub(",", "", col_id, fixed = T)]
hpa[, col_id := gsub("-", "_", col_id, fixed = T)]
#convert from long to wide format
hpa <- dcast(hpa, protein_id ~ col_id, fun.aggregate = getmode, value.var = "level", drop = T, fill = "not detected")
replace_na(hpa, 2:ncol(hpa), "not detected")
#convert HPA qualitative expression data to sparce matrix object and do one-hot encoding of factor variables
hpa.sparse.matrix <- sparse.model.matrix(~.-1, data = hpa)
#convert back to data.table and merge with the rest of the data matrix
hpa <- as.data.table(as.matrix(hpa.sparse.matrix), keep.rownames = F)
dt <- merge(dt, hpa, by.x = "id1", by.y = "protein_id", all.x = T)
rm(hpa)

#Extract expression data from CCLE
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
ccle <- dbGetQuery(conn, "select protein_id,cell_id,tissue,expression from ccle")
dbDisconnect(conn)
rm(conn)
setDT(ccle)

ccle[is.na(tissue), col_id := cell_id]
ccle[!is.na(tissue), col_id := sprintf("%s_%s", cell_id,tissue)]
ccle[, `:=`(tissue = NULL, cell_id = NULL)]
#convert from long to wide format
ccle <- dcast(ccle, protein_id ~ col_id, fun.aggregate = median, value.var = "expression", drop = T, fill = 0)
#scale the data
scale.data.table(ccle)
# merge with the rest of the data matrix by database protein_id (id1 in data.table)
dt <- merge(dt, ccle, by.x = "id1", by.y = "protein_id", all.x = T)
rm(ccle)


# replace NA in numeric columns with 0, future iteration use multiple imputations package mice
col.ids <- which(sapply(dt, is.numeric))
col.ids <- col.ids[2:length(col.ids)]
replace_na(dt, col.ids, 0)
#set data.table attribute with OMIM Phenotype name, id, and date stamp
setattr(dt, "pheno_series_title", pheno_series)
setattr(dt, "omim_phen_series_id", g2d[i, omim_phen_series_id])
setattr(dt, "Date", Sys.Date())
saveRDS(dt, sprintf("data/input/%s.rds", g2d[i, omim_phen_series_id]))
cat("saved data set for", pheno_series,"\n")

