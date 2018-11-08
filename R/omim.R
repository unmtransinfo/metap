#!/usr/bin/env Rscript

library(data.table)
library(RPostgreSQL)

download.file("https://data.omim.org/downloads/FIHB_S7QSr6ELan1-W3G3g/mimTitles.txt", destfile = "data/omim/mimTitles.txt")
omim <- read.delim2("data/omim/mimTitles.txt", skip = 2, header = T, comment.char = "#", stringsAsFactors = F, na.strings = c("NULL", "", "NA"), quote = "\"", col.names = c("Prefix", "Mim_Number", "Preferred_Title_symbol","Alternative_Title_symbol","Included_Title_symbol"))
setDT(omim)
#omim <- fread("data/omim/mimTitles.txt", skip = 2, header = T, sep = "\t", quote = "", na.strings = "NULL", strip.white = F, col.names = c("Prefix", "Mim_Number", "Preferred_Title_symbol","Alternative_Title_symbol","Included_Title_symbol"), comment.char = "#")
omim[Prefix == "NULL", Prefix := NA]
omim <- omim[, .(Mim_Number, Preferred_Title_symbol)]
setnames(omim, c("Mim_Number", "Preferred_Title_symbol"), c("mim", "title"))
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
dbSendQuery(conn, "delete from omim")
dbWriteTable(conn, "omim", omim, append = T, row.names = F)

dbCommit(conn)
ps.all <- read.delim("https://www.omim.org/phenotypicSeriesTitle/all?format=tsv", header = T, sep = "\t", skip = 4, quote = "", stringsAsFactors = F, col.names = c("title", "omim_phen_series_id"))
setDT(ps.all)
dbSendQuery(conn, "delete from omim_phen_series")
dbWriteTable(conn, "omim_phen_series", ps.all, append = T, row.names = F)

dbSendQuery(conn, "delete from omim2phen_series")
for(i in 1:nrow(ps.all)) {
  ps <- read.delim(sprintf("https://www.omim.org/phenotypicSeries/%s?format=tsv", ps.all[i, omim_phen_series_id]), header = T, skip = 6, sep = "\t", quote = "", skipNul = T)
  setDT(ps)
  ps <- ps[!is.na(Phenotype.MIM.number), .(Phenotype.MIM.number)]
  ps[, omim_phen_series_id := ps.all[i, omim_phen_series_id]]
  setnames(ps, "Phenotype.MIM.number", "mim")
  ps <- unique(ps)
  dbWriteTable(conn, "omim2phen_series", ps, append = T, row.names = F)
}

dbDisconnect(conn)
rm(conn)
