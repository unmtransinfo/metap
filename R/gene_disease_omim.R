#!/usr/bin/env Rscript

library(data.table)
library(RPostgreSQL)
library(Matrix)

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

#connects to Metapath database instance on localhost with "oleg" user, replace with actual username
conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
#extract all OMIM phenotype series to OMIM mappings
omim.ps <- dbGetQuery(conn, "select * from omim2phen_series")
#extract all UniProt (are the same as OMIM) diseases
uniprot.disease <- dbGetQuery(conn, "select * from uniprot_disease");
#extract all OMIM phenotype series terms
ps.name <- dbGetQuery(conn, "select * from omim_phen_series");
#extract all disease/phenotype gene associations from ClinVar
clinvar <- dbGetQuery(conn, "select distinct clinvar.protein_id,clinvar.clinical_significance,clinvar_disease.cv_dis_id,phenotype,source,source_id from clinvar_disease,clinvar_disease_xref,clinvar WHERE clinvar_disease.cv_dis_id=clinvar_disease_xref.cv_dis_id and source = 'OMIM' and clinvar.cv_dis_id=clinvar_disease.cv_dis_id and clinical_significance in ('Pathogenic, Affects','Benign, protective, risk factor','Pathogenic/Likely pathogenic','Pathogenic/Likely pathogenic, other','Pathogenic, other','Affects','Pathogenic, other, protective','Conflicting interpretations of pathogenicity, Affects, association, other','Pathogenic/Likely pathogenic, drug response','Pathogenic, risk factor','risk factor','Pathogenic, association','Conflicting interpretations of pathogenicity, Affects, association, risk factor','Pathogenic/Likely pathogenic, risk factor','Affects, risk factor','Conflicting interpretations of pathogenicity, association, other, risk factor','Likely pathogenic, association','association, protective','Likely pathogenic, Affects','Pathogenic','Conflicting interpretations of pathogenicity, association','Pathogenic/Likely pathogenic, Affects, risk factor','Conflicting interpretations of pathogenicity, other, risk factor','association, risk factor','Benign, protective','Conflicting interpretations of pathogenicity, risk factor','Uncertain significance, protective','association','Uncertain significance, Affects','protective, risk factor','Pathogenic, association, protective','Pathogenic, protective','Likely pathogenic, other','Pathogenic, protective, risk factor','Benign, association, protective','Conflicting interpretations of pathogenicity, Affects','Benign/Likely benign, protective','protective')")
#extract rat phenotype - gene associations and mappings from rat to human proteins based on RGD provided mapping
rat.terms <- dbGetQuery(conn, "select DISTINCT rat_term.rgd_id, protein_id as rat_protein_id, rat2human.human_protein_id, term_id, rat_term.term_name, rdo_xref.db,rdo_xref.xref from rat_term,protein2rgd,rdo_xref,rat2human where ontology = 'RGD Disease Ontology' and rat_term.rgd_id = protein2rgd.rgd_id and rat_term.term_id=rdo_xref.rdo_id and rdo_xref.db='OMIM' and protein2rgd.protein_id = rat2human.rat_protein_id and rat_term.qualifier not in ('no_association', 'NOT')")
#extract the list of human proteins
protein <- dbGetQuery(conn, "select * from protein where tax_id = 9606");
dbDisconnect(conn)
rm(conn)

#convert all query results to data.table
setDT(omim.ps)
setDT(uniprot.disease)
setDT(ps.name)
setDT(clinvar)
setDT(rat.terms)
setDT(protein)

#map OMIM phenotypes in UniProt to phenotypic series id
uniprot.disease <- merge(uniprot.disease, omim.ps, by.x = "ref_id", by.y = "mim", all.x = T)
#map OMIM phenotypic series id to names
uniprot.disease <- merge(uniprot.disease, ps.name, by.x = "omim_phen_series_id", by.y = "omim_phen_series_id", all.x = T)
#change column names to match OMIM column names
setnames(uniprot.disease, c("ref_id", "term", "title"), c("mim", "mim_title", "pheno_series_title"))
uniprot.disease[, `:=`(disease_id = NULL, dbref = NULL)]

#map OMIM phenotypes in ClinVar to phenotypic series id
clinvar$source_id <- as.integer(clinvar$source_id)
clinvar <- merge(clinvar, omim.ps, by.x = "source_id", by.y = "mim", all.x = T)
#map OMIM phenotypic series id to names
clinvar <- merge(clinvar, ps.name, by.x = "omim_phen_series_id", by.y = "omim_phen_series_id", all.x = T)
#change column names to match OMIM column names
setnames(clinvar, c("source_id","phenotype","title"), c("mim","mim_title","pheno_series_title"))
clinvar[, `:=`(clinical_significance = NULL, cv_dis_id = NULL, source = NULL)]
clinvar <- unique(clinvar)

#map OMIM phenotypes in RGD database to phenotypic series id
rat.terms$xref <- as.integer(rat.terms$xref)
rat.terms <- merge(rat.terms, omim.ps, by.x = "xref", by.y = "mim", all.x = T)
#map OMIM phenotypic series id to names
rat.terms <- merge(rat.terms, ps.name, by.x = "omim_phen_series_id", by.y = "omim_phen_series_id", all.x = T)
#change column names to match OMIM column names
setnames(rat.terms, c("xref","term_name","title", "human_protein_id"), c("mim","mim_title","pheno_series_title", "protein_id"))
rat.terms[, `:=`(rgd_id = NULL, rat_protein_id = NULL, term_id = NULL, db = NULL)]
rat.terms <- unique(rat.terms)

#concatenate gene - OMIM phenotype associations from ClinVar, UniProt, and RGD
g2d <- rbindlist(list(uniprot.disease, clinvar, rat.terms), use.names = T)
g2d <- unique(g2d)

#Select OMIM phentypic series associated with at least 50 genes
top.phentype.series <- g2d[!is.na(pheno_series_title), .(uq_prot = uniqueN(protein_id)), by = pheno_series_title][uq_prot >= 50]
#Add OMIM phentypic series ids
top.phentype.series <- merge(top.phentype.series, unique(g2d[, .(pheno_series_title, omim_phen_series_id)]), by.x = "pheno_series_title", by.y = "pheno_series_title")

#extract data matrix for each of the top (>= genes associated) OMIM phenotypic series
for(i in 1:nrow(top.phentype.series)) {
  pheno_series <- top.phentype.series[i, pheno_series_title]
  #separate proteins into 2 sets one known to be associated with OMIM phenotype based on available data and the other one not known
  #assign to right side protein_ids known to be associated with OMIM phenotype
  right.side <- g2d[pheno_series_title == pheno_series, unique(protein_id)]
  #assign to left side protein_ids known to be associated with OMIM phenotype
  left.side <- g2d[pheno_series_title != pheno_series & !protein_id %in% right.side, unique(protein_id)]
  #all the other protein in human proteome are the external set and will apply XGBoost models to predict potential associations
  test.idx <- protein[!protein_id %in% c(left.side, right.side), protein_id]
  
  #initialize the data matrix for machine learning with protein_ids, subset annotations (train/test), flag positive/negative associations
  dt <- data.table(id1 = c(right.side, left.side, test.idx), Y = c(rep_len("pos", length.out = length(right.side)), rep_len("neg", length.out = length(left.side)), rep_len("neg", length.out = length(test.idx))), subset = c(rep_len("train", length.out = length(right.side)), rep_len("train", length.out = length(left.side)), rep_len("test", length.out = length(test.idx))))
  dt$Y <- factor(dt$Y, levels = c("neg", "pos"))
  dt$subset <- factor(dt$subset, levels = c("train", "test"))
  
  #all proteins not associated with the disease/phenotype are on the left side
  left.side <- c(left.side, test.idx)

  #extract all PPI scores from StringDB for proteins on the right side
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
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
  setattr(dt, "omim_phen_series_id", top.phentype.series[i, omim_phen_series_id])
  setattr(dt, "Date", Sys.Date())
  saveRDS(dt, sprintf("data/input/%s.rds", top.phentype.series[i, omim_phen_series_id]))
  cat("saved data set for", pheno_series,"\n")
}

#the steps below are equvalent to the steps above except these are applied to individual OMIM Phenotypes 

top.phenotype <- g2d[!is.na(mim_title), .(uq_prot = uniqueN(protein_id)), by = mim_title][uq_prot >= 50]
top.phenotype <- merge(top.phenotype, unique(g2d[, .(mim, mim_title)]), by.x = "mim_title", by.y = "mim_title")
top.phenotype <- unique(top.phenotype, by = "mim_title")

for(i in 1:nrow(top.phenotype)) {
  phenotype <- top.phenotype[i, mim_title]
  right.side <- g2d[mim_title == phenotype, unique(protein_id)]
  left.side <- g2d[mim_title != phenotype & !protein_id %in% right.side, unique(protein_id)]
  test.idx <- protein[!protein_id %in% c(left.side, right.side), protein_id]
  
  dt <- data.table(id1 = c(right.side, left.side, test.idx), Y = c(rep_len("pos", length.out = length(right.side)), rep_len("neg", length.out = length(left.side)), rep_len("neg", length.out = length(test.idx))), subset = c(rep_len("train", length.out = length(right.side)), rep_len("train", length.out = length(left.side)), rep_len("test", length.out = length(test.idx))))
  dt$Y <- factor(dt$Y, levels = c("neg", "pos"))
  dt$subset <- factor(dt$subset, levels = c("train", "test"))
  
  left.side <- c(left.side, test.idx)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  stringdb <- dbGetQuery(conn, sprintf("select * from stringdb_score where protein_id1 in (%s) or protein_id2 in (%s)", paste(right.side, collapse = ","), paste(right.side, collapse = ",")))
  dbDisconnect(conn)
  rm(conn)
  setDT(stringdb)
  
  stringdb[protein_id2 %in% right.side, id2 := protein_id2]
  stringdb[is.na(id2), `:=`(id2 = protein_id1, id1 = protein_id2)]
  stringdb[is.na(id1), id1 := protein_id1]
  stringdb[, `:=`(protein_id1 = NULL, protein_id2 = NULL)]
  stringdb <- stringdb[id1 %in% c(right.side, left.side)]
  stringdb[, id2.count := .N, by = id2]
  stringdb[, id1.count := .N, by = id1]
  disease.conn.count <- stringdb[, uniqueN(id2)]
  stringdb[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)*combined_score]
  stringdb[, col.id := sprintf("pp.%d", id2)]
  stringdb <- dcast(stringdb, id1 ~ col.id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
  scale.data.table(stringdb)
  dt <- merge(dt, stringdb, by.x = "id1", by.y = "id1", all.x = T)
  rm(stringdb)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  reactome <- dbGetQuery(conn, sprintf("select * from reactomea where reactome_id in (select distinct reactome_id from reactomea where protein_id in (%s))", paste(right.side, collapse = ",")))
  dbDisconnect(conn)
  rm(conn)
  setDT(reactome)
  
  reactome[, `:=`(evidence = NULL)]
  reactome <- unique(reactome)
  reactome <- reactome[protein_id %in% c(left.side, right.side)]
  reactome.right.side <- reactome[protein_id %in% right.side]
  reactome.left.side <- reactome[protein_id %in% c(right.side, left.side)]
  setnames(reactome.left.side, "protein_id", "id1")
  setnames(reactome.right.side, "protein_id", "id2")
  reactome <- merge(reactome.left.side, reactome.right.side, by.x = "reactome_id", by.y = "reactome_id", allow.cartesian = T)
  reactome <- reactome[id1 != id2]
  disease.conn.count <- reactome[, uniqueN(id2)]
  reactome[, id2.count := .N, by = id2]
  reactome[, id1.count := .N, by = id1]
  reactome[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
  reactome <- dcast(reactome, id1 ~ reactome_id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
  scale.data.table(reactome)
  dt <- merge(dt, reactome, by.x = "id1", by.y = "id1", all.x = T)
  rm(reactome)
  
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  kegg <- dbGetQuery(conn, sprintf("select protein_id,kegg_pathway_id from kegg_pathway where kegg_pathway_id in (select distinct kegg_pathway_id from kegg_pathway where protein_id in (%s))", paste(right.side, collapse = ",")))
  dbDisconnect(conn)
  rm(conn)
  setDT(kegg)
  
  kegg <- kegg[protein_id %in% c(left.side, right.side)]
  kegg.right.side <- kegg[protein_id %in% right.side]
  kegg.left.side <- kegg[protein_id %in% c(right.side, left.side)]
  setnames(kegg.left.side, "protein_id", "id1")
  setnames(kegg.right.side, "protein_id", "id2")
  kegg <- merge(kegg.left.side, kegg.right.side, by.x = "kegg_pathway_id", by.y = "kegg_pathway_id", allow.cartesian = T)
  kegg <- kegg[id1 != id2]
  disease.conn.count <- kegg[, uniqueN(id2)]
  kegg[, id2.count := .N, by = id2]
  kegg[, id1.count := .N, by = id1]
  kegg[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
  kegg <- dcast(kegg, id1 ~ kegg_pathway_id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
  scale.data.table(kegg)
  dt <- merge(dt, kegg, by.x = "id1", by.y = "id1", all.x = T)
  rm(kegg)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  gtex <- dbGetQuery(conn, "select protein_id,median_tpm,tissue_type_detail from gtex")
  dbDisconnect(conn)
  rm(conn)
  setDT(gtex)
  
  gtex <- dcast(gtex, protein_id ~ tissue_type_detail, fun.aggregate = median, value.var = "median_tpm", drop = T, fill = 0)
  scale.data.table(gtex)
  dt <- merge(dt, gtex, by.x = "id1", by.y = "protein_id", all.x = T)
  rm(gtex)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  lincs <- dbGetQuery(conn, "select protein_id,pert_id||':'||cell_id as col_id,zscore from lincs")
  dbDisconnect(conn)
  rm(conn)
  setDT(lincs)
  
  lincs <- dcast(lincs, protein_id ~ col_id, fun.aggregate = median, value.var = "zscore", drop = T, fill = 0)
  scale.data.table(lincs)
  dt <- merge(dt, lincs, by.x = "id1", by.y = "protein_id", all.x = T)
  rm(lincs)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  go <- dbGetQuery(conn, sprintf("select distinct protein_id,go_id from goa where go_id in (select distinct go_id from goa where protein_id in (%s))", paste(right.side, collapse = ",")))
  dbDisconnect(conn)
  rm(conn)
  setDT(go)
  
  go <- go[protein_id %in% c(left.side, right.side)]
  go.right.side <- go[protein_id %in% right.side]
  go.left.side <- go[protein_id %in% c(right.side, left.side)]
  setnames(go.left.side, "protein_id", "id1")
  setnames(go.right.side, "protein_id", "id2")
  go <- merge(go.left.side, go.right.side, by.x = "go_id", by.y = "go_id", allow.cartesian = T)
  go <- go[id1 != id2]
  disease.conn.count <- go[, uniqueN(id2)]
  go[, id2.count := .N, by = id2]
  go[, id1.count := .N, by = id1]
  go[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
  go <- dcast(go, id1 ~ go_id, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
  scale.data.table(go)
  dt <- merge(dt, go, by.x = "id1", by.y = "id1", all.x = T)
  rm(go)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  interpro <- dbGetQuery(conn, sprintf("select distinct protein_id,entry_ac from interproa where entry_ac in (select distinct entry_ac from interproa where protein_id in (%s))", paste(right.side, collapse = ",")))
  dbDisconnect(conn)
  rm(conn)
  setDT(interpro)
  
  interpro <- interpro[protein_id %in% c(left.side, right.side)]
  interpro.right.side <- interpro[protein_id %in% right.side]
  interpro.left.side <- interpro[protein_id %in% c(right.side, left.side)]
  setnames(interpro.left.side, "protein_id", "id1")
  setnames(interpro.right.side, "protein_id", "id2")
  interpro <- merge(interpro.left.side, interpro.right.side, by.x = "entry_ac", by.y = "entry_ac", allow.cartesian = T)
  interpro <- interpro[id1 != id2]
  disease.conn.count <- interpro[, uniqueN(id2)]
  interpro[, id2.count := .N, by = id2]
  interpro[, id1.count := .N, by = id1]
  interpro[, pdp := (id1.count^-0.5)*(id2.count^-0.5)*(disease.conn.count^-0.5)]
  interpro <- dcast(interpro, id1 ~ entry_ac, fun.aggregate = mean, value.var = "pdp", drop = T, fill = 0)
  scale.data.table(interpro)
  dt <- merge(dt, interpro, by.x = "id1", by.y = "id1", all.x = T)
  rm(interpro)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  hpa <- dbGetQuery(conn, "select protein_id, tissue||'.'||cell_type as col_id,level from hpa_norm_tissue where reliability in ('supported','approved')")
  dbDisconnect(conn)
  rm(conn)
  setDT(hpa)
  
  hpa$level <- factor(hpa$level, levels = c("not detected", "low", "medium", "high"), ordered = F)
  hpa <- unique(hpa)
  hpa[, col_id := gsub(" ", "_", col_id, fixed = T)]
  hpa[, col_id := gsub("/", "_", col_id, fixed = T)]
  hpa[, col_id := gsub(",", "", col_id, fixed = T)]
  hpa[, col_id := gsub("-", "_", col_id, fixed = T)]
  hpa <- dcast(hpa, protein_id ~ col_id, fun.aggregate = getmode, value.var = "level", drop = T, fill = "not detected")
  replace_na(hpa, 2:ncol(hpa), "not detected")
  hpa.sparse.matrix <- sparse.model.matrix(~.-1, data = hpa)
  hpa <- as.data.table(as.matrix(hpa.sparse.matrix), keep.rownames = F)
  dt <- merge(dt, hpa, by.x = "id1", by.y = "protein_id", all.x = T)
  rm(hpa)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", host = "127.0.0.1", dbname = "metap")
  ccle <- dbGetQuery(conn, "select protein_id,cell_id,tissue,expression from ccle")
  dbDisconnect(conn)
  rm(conn)
  setDT(ccle)
  
  ccle[is.na(tissue), col_id := cell_id]
  ccle[!is.na(tissue), col_id := sprintf("%s_%s", cell_id,tissue)]
  ccle[, `:=`(tissue = NULL, cell_id = NULL)]
  ccle <- dcast(ccle, protein_id ~ col_id, fun.aggregate = median, value.var = "expression", drop = T, fill = 0)
  scale.data.table(ccle)
  dt <- merge(dt, ccle, by.x = "id1", by.y = "protein_id", all.x = T)
  rm(ccle)
  
  
  # replace NA in numeric columns with 0, future iteration use multiple imputations package mice
  col.ids <- which(sapply(dt, is.numeric))
  col.ids <- col.ids[2:length(col.ids)]
  replace_na(dt, col.ids, 0)
  setattr(dt, "omim_title", phenotype)
  setattr(dt, "mim_id", top.phenotype[i, mim])
  setattr(dt, "Date", Sys.Date())
  saveRDS(dt, sprintf("data/input/%s.rds", top.phenotype[i, mim]))
  cat("saved data set for", phenotype,"\n")
}