#!/usr/bin/env Rscript

library(data.table)
library(KEGGREST)

hsa.gene <- keggList("hsa")
hsa.gene.id <- names(hsa.gene)
data <- data.table()
for(i in seq(1, length(hsa.gene.id), 10)) {
  end <- i + 9
  hsa.entries <- keggGet(hsa.gene.id[i:end])
  for(entry in hsa.entries) {
    if(!"PATHWAY" %in% names(entry)) {
      next
    }
    pathways <- entry$PATHWAY
    if("DBLINKS" %in% names(entry)) {
      for(link in entry$DBLINKS) {
        idx <- which(startsWith(entry$DBLINKS, "NCBI-GeneID:"))
        if(length(idx) > 0) {
          ncbi_gene_id <- sub("NCBI-GeneID: ", "", entry$DBLINKS[idx[1]], fixed = T)
          for(p in 1:length(pathways)) {
            data <- rbindlist(list(data, data.table(ncbi_gene_id = ncbi_gene_id, kegg_pathway_id = names(pathways)[p], kegg_pathway_name = pathways[p])))
          }
        }
      }
    }
    cat("processed entry:", entry$NAME, "\n")
  }
}

data[, ncbi_gene_id := as.integer(ncbi_gene_id)]
data <- unique(data, by = c("ncbi_gene_id", "kegg_pathway_id"))
fwrite(data, file = "data/kegg/kegg.pathway.notmapped.tsv", sep = "\t", quote = T, row.names = F, col.names = T, na = "None")

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,gene_id from ncbi")
dbDisconnect(conn)
rm(conn)
setDT(protein)

data <- merge(data, protein, by.x = "ncbi_gene_id", by.y = "gene_id")
data <- unique(data, by = c("kegg_pathway_id", "protein_id"))

fwrite(data[, .(protein_id, kegg_pathway_id, kegg_pathway_name)], file = "data/kegg/pathway.tsv", sep = "\t", quote = T, row.names = F, col.names = T, na = "None")

if(file.exists("data/kegg/pathway.tsv.gz")) {
  file.remove("data/kegg/pathway.tsv.gz")
}
system("gzip -9v data/kegg/pathway.tsv")
