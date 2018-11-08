#!/usr/bin/env Rscript

library(data.table)
library(stringdist)
library(RPostgreSQL)
library(ADPclust)
library(cluster)

#conn <- dbConnect(PostgreSQL(), user = "oleg", host = "10.234.37.25", dbname = "metap")
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
rat.term <- dbGetQuery(conn, "select term_id,term_name from rat_term")
dbDisconnect(conn)
rm(conn)
setDT(rat.term)

cat("pulled",nrow(rat.term),"rat terms from database\n")

rat.term <- unique(rat.term, by = c("term_id", "term_name"))

rat.term[, term_name := gsub("[[:punct:]]\\s+", " ", term_name)]
rat.term[, term_name := gsub("\\s+[[:punct:]]", " ", term_name)]
rat.term[, term_name := sub("[[:punct:]]$", "", term_name)]
rat.term[, term_name := sub("^[[:punct:]]$", "", term_name)]
rat.term[, term_name := gsub("\\s[[:graph:]]{1,3}\\s", " ", term_name)]
rat.term[, term_name := gsub("\\s[[:graph:]]{1,3}\\s", " ", term_name)]
rat.term[, term_name := sub("\\s[[:graph:]]{1,3}$", "", term_name)]
rat.term[, term_name := sub("^[[:graph:]]{1,3}\\s+", "", term_name)]
rat.term[, term_name := gsub("\\s+", " ", term_name)]
rat.term[, term_name := gsub("\\s+$", "", term_name)]

dm <- stringdistmatrix(rat.term[, tolower(term_name)], nthread = 6, method = "cosine")

cat("computed disimilarity matrix\n")

dv <- diana(dm)

cat("computed clustering\n")

max.sil <- 0
optim.k <- 0

for(i in seq(1000, nrow(rat.term) - 100, 100)) {
  cluster <- cutree(as.hclust(dv), k = i)
  sil <- silhouette(cluster, dm)
  ave.sil <- mean(sil[,3])
  if(ave.sil > max.sil) {
    max.sil <- ave.sil
    optim.k <- i
  }
  cat("i=",i,"ave silhouette=",ave.sil,"\n")
}

cat("max silhouette=",max.sil,"optimal k=",optim.k,"\n")

dv.clusters <- cutree(as.hclust(dv), k = optim.k) # based on silhouette index
rat.term[, clust_id := dv.clusters]
rat.term[, sql := sprintf("update rat_term set clust_id = %d where term_id = '%s';", clust_id, term_id)]

fwrite(rat.term[, .(sql)], file = "rat_term.clust.sql", sep = "\t", col.names = F, row.names = F, na = "", quote = F)
cat("Done\n")
