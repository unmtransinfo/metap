#!/usr/bin/env Rscript

library(data.table)
library(RPostgreSQL)
library(Matrix)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

replace_na <- function(DT, ind, replacement) {
  for(j in ind){
    set(DT, i = which(is.na(DT[[j]])), j = j, value = replacement)
  }
}

scale.data.table <- function(dt) {
  col.names <- colnames(dt)[2:ncol(dt)]
  dt[, (col.names) := lapply(.SD, scale), .SDcols=col.names]
}

onto.level.list <- function(dt) {
  lev <- 0
  onto <- copy(dt)
  onto[mp_term_id == "MP_0000001", level := lev]
  children.count <- onto[parent_id %chin% onto[level == lev, mp_term_id], .N]
  while (children.count > 0) {
    lev <- lev + 1
    onto[parent_id %chin% onto[level == lev - 1, mp_term_id], level := lev]
    children.count <- onto[parent_id %chin% onto[level == lev, mp_term_id], .N]
  }
  terms <- onto[level > 1, mp_term_id]
  extended.parents <- data.table()
  for (term in terms) {
    parent <- onto[mp_term_id == term, parent_id]
    indirect.parent <- onto[mp_term_id == parent, parent_id]
    parents <- indirect.parent
    while (length(indirect.parent) > 0) {
      indirect.parent <- onto[mp_term_id == indirect.parent & !is.na(parent_id), parent_id]
      parents <- c(parents, indirect.parent)
    }
    t <- data.table(mp_term_id = rep(term, length(parents)), parent_id = parents, level = rep(onto[mp_term_id == term, level], length(parents)), name = rep(onto[mp_term_id == term, name], length(parents)))
    extended.parents <- rbindlist(list(extended.parents, t), use.names = T)
  }
  onto <- rbindlist(list(onto, extended.parents), use.names = T)
  return(onto)
}

conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
mp.onto <- dbGetQuery(conn, "SELECT * from mp_onto")
mouse.phen <- dbGetQuery(conn, "select * from mousephenotype")
mouse2human <- dbGetQuery(conn, "select h.protein_id human_protein_id,m.protein_id mouse_protein_id FROM homology h INNER JOIN homology m on h.homologene_group_id = m.homologene_group_id and h.tax_id = 9606 and m.tax_id = 10090")
rat.term <- dbGetQuery(conn, "select * from rat_term WHERE ontology = 'Mammalian Phenotype'")
rat2human <- dbGetQuery(conn, "select * from rat2human")
rgd2protein <- dbGetQuery(conn, "select * from protein2rgd")
protein <- dbGetQuery(conn, "select * from protein where tax_id = 9606");
dbDisconnect(conn)
rm(conn)

setDT(mp.onto)
setDT(mouse.phen)
setDT(mouse2human)
setDT(rat.term)
setDT(rat2human)
setDT(rgd2protein)
setDT(protein)

mp.onto.long <- onto.level.list(mp.onto)

mouse.phen <- merge(mouse.phen, mouse2human, by.x = "protein_id", by.y = "mouse_protein_id", sort = F)
mouse.phen <- mouse.phen[, .(human_protein_id, mp_term_id, association)]
rat.term <- merge(rat.term, rgd2protein, by.x = "rgd_id", by.y = "rgd_id", sort = F)
rat.term <- merge(rat.term, rat2human, by.x = "protein_id", by.y = "rat_protein_id", sort = F)
rat.term <- rat.term[, .(human_protein_id, term_id, qualifier)]
rat.term[!is.na(qualifier), association := F]
rat.term[is.na(qualifier), association := T]
rat.term[, qualifier := NULL]
setnames(rat.term, "term_id", "mp_term_id")
rat.term[, mp_term_id := sub(":", "_", mp_term_id, fixed = T)]

g2d <- rbindlist(list(rat.term, mouse.phen), use.names = T)
g2d <- unique(g2d)
setnames(g2d, "human_protein_id", "protein_id")
g2d <- merge(g2d, mp.onto, by.x = "mp_term_id", by.y = "mp_term_id", sort = F)
g2d <- merge(g2d, unique(mp.onto.long, by = "mp_term_id")[, .(mp_term_id, level)], by.x = "mp_term_id", by.y = "mp_term_id", sort = F)
g2d <- unique(g2d)
g2d <- g2d[level > 1]
g2d.count <- g2d[, .(prot_count = uniqueN(protein_id)), by = .(association, parent_id)]
g2d.count <- dcast(g2d.count, parent_id ~ association, value.var = "prot_count", fun.aggregate = sum)
g2d.count <- g2d.count[`TRUE` >= 50 & `FALSE` >= 50]

for(term_id in g2d.count[, unique(parent_id)]) {
  right.side <- g2d[parent_id == term_id & association == T, unique(protein_id)]
  left.side <- g2d[parent_id == term_id & association == F, unique(protein_id)]
  test.idx <- protein[!protein_id %in% c(left.side, right.side), protein_id]
  
  dt <- data.table(id1 = c(right.side, left.side, test.idx), Y = c(rep_len("pos", length.out = length(right.side)), rep_len("neg", length.out = length(left.side)), rep_len("neg", length.out = length(test.idx))), subset = c(rep_len("train", length.out = length(right.side)), rep_len("train", length.out = length(left.side)), rep_len("test", length.out = length(test.idx))))
  dt$Y <- factor(dt$Y, levels = c("neg", "pos"))
  dt$subset <- factor(dt$subset, levels = c("train", "test"))
  
  left.side <- c(left.side, test.idx)

  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
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
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
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
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
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
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
  gtex <- dbGetQuery(conn, "select protein_id,median_tpm,tissue_type_detail from gtex")
  dbDisconnect(conn)
  rm(conn)
  setDT(gtex)
  
  gtex <- dcast(gtex, protein_id ~ tissue_type_detail, fun.aggregate = median, value.var = "median_tpm", drop = T, fill = 0)
  scale.data.table(gtex)
  dt <- merge(dt, gtex, by.x = "id1", by.y = "protein_id", all.x = T)
  rm(gtex)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
  lincs <- dbGetQuery(conn, "select protein_id,pert_id||':'||cell_id as col_id,zscore from lincs")
  dbDisconnect(conn)
  rm(conn)
  setDT(lincs)
  
  lincs <- dcast(lincs, protein_id ~ col_id, fun.aggregate = median, value.var = "zscore", drop = T, fill = 0)
  scale.data.table(lincs)
  dt <- merge(dt, lincs, by.x = "id1", by.y = "protein_id", all.x = T)
  rm(lincs)
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
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
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
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
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
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
  
  conn <- dbConnect(PostgreSQL(), user = "oleg", dbname = "metap")
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
  setattr(dt, "MP_TERM_NAME", mp.onto[mp_term_id == term_id, name])
  setattr(dt, "MP_TERM_ID", term_id)
  setattr(dt, "Date", Sys.Date())
  saveRDS(dt, sprintf("data/input/%s.rds", term_id))
  cat("saved data set for", term_id,"\n")
}