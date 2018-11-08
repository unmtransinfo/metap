#!/usr/bin/env Rscript

library(RJDBC)
library(data.table)

drv <- JDBC("org.apache.derby.jdbc.EmbeddedDriver","/Users/oleg/.m2/repository/org/apache/derby/derby/10.13.1.1/derby-10.13.1.1.jar", identifier.quote="\"")
conn <- dbConnect(drv, "jdbc:derby:/Users/oleg/Documents/dbase/drugdb/.config/localdb/db")
mol.name <- dbGetQuery(conn, "select id,name from STRUCTURES")
parent.name <- dbGetQuery(conn, "select cd_id,name from PARENTMOL")
dbDisconnect(conn)
rm(conn, drv)

setDT(mol.name)
setDT(parent.name)

parent.name[, ID := sprintf("%d.P", CD_ID)]
parent.name[, CD_ID := NULL]

names <- rbindlist(list(mol.name, parent.name), use.names = T)

fwrite(names, file = "data/drugcentral/drug_names.tsv", col.names = T, sep = "\t", quote = T, row.names = F, na = "None")