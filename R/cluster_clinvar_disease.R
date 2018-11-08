#!/usr/bin/env Rscript

library(data.table)
library(stringdist)
library(RPostgreSQL)
library(ADPclust)
library(cluster)

conn <- dbConnect(PostgreSQL(), user = "oleg", host = "10.234.37.25", dbname = "metap")
#conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
cv.dis <- dbGetQuery(conn, "select * from clinvar_disease")
dbDisconnect(conn)
rm(conn)
setDT(cv.dis)

cat("pulled",nrow(cv.dis),"disease terms from database\n")

cv.dis[, phenotype := gsub("[[:punct:]]\\s+", " ", phenotype)]
cv.dis[, phenotype := gsub("\\s+[[:punct:]]", " ", phenotype)]
cv.dis[, phenotype := sub("[[:punct:]]$", "", phenotype)]
cv.dis[, phenotype := sub("^[[:punct:]]$", "", phenotype)]
cv.dis[, phenotype := gsub("\\s[[:graph:]]{1,3}\\s", " ", phenotype)]
cv.dis[, phenotype := gsub("\\s[[:graph:]]{1,3}\\s", " ", phenotype)]
cv.dis[, phenotype := sub("\\s[[:graph:]]{1,3}$", "", phenotype)]
cv.dis[, phenotype := sub("^[[:graph:]]{1,3}\\s+", "", phenotype)]
cv.dis[, phenotype := gsub("\\s+", " ", phenotype)]
cv.dis[, phenotype := gsub("\\s+$", "", phenotype)]


dm <- stringdistmatrix(cv.dis[, tolower(phenotype)], nthread = 6, method = "cosine")

cat("computed disimilarity matrix\n")

dv <- diana(dm)

cat("computed clustering\n")

max.sil <- 0
optim.k <- 0

for(i in seq(1000, nrow(cv.dis) - 100, 100)) {
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
cv.dis[, clust_id := dv.clusters]
cv.dis[, sql := sprintf("update clinvar_disease set clust_id = %d where cv_dis_id = %d;", clust_id, cv_dis_id)]

fwrite(cv.dis[, .(sql)], file = "clinvar.dis.clust.sql", sep = "\t", col.names = F, row.names = F, na = "", quote = F)
cat("Done\n")
