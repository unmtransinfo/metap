#!/usr/bin/env Rscript

library(data.table)
library(ontologyIndex)

go <- get_ontology("http://purl.obolibrary.org/obo/go.obo", extract_tags = "everything")
go.dt <- data.table(id = names(go$name), name = go$name)
go.dt <- go.dt[id %like% "GO:"]
aspect <- unlist(go$namespace)
aspect.dt <- data.table(id = names(aspect), aspect = aspect)
go.dt <- merge(go.dt, aspect.dt, by.x = "id", by.y = "id", all.x = T, all.y = F)
def.dt <- data.table(id = names(go$def), def = go$def)
go.dt <- merge(go.dt, def.dt, by.x = "id", by.y = "id", all.x = T, all.y = F)
go.dt[aspect == "molecular_function", aspect.short := "F"]
go.dt[aspect == "biological_process", aspect.short := "B"]
go.dt[aspect == "cellular_component", aspect.short := "C"]
go.dt[, aspect := NULL]
go.dt[, def := gsub("\"", "", def, fixed = T)]
fwrite(go.dt, file = "data/go/go.tsv", sep = "\t", quote = T, row.names = F, col.names = T, na = "None")

download.file("http://geneontology.org/gene-associations/goa_human.gaf.gz", destfile = "data/go/goa_human.gaf.gz")
goa <- fread("gzcat data/go/goa_human.gaf.gz", header = F, sep = "\t", skip = "UniProtKB", quote = "", col.names = c("DB","DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB_Reference", "Evidence_Code", "With", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"))

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein where tax_id = 9606")
dbDisconnect(conn)
rm(conn)
setDT(protein)
goa <- merge(goa, protein, by.x = "DB_Object_ID", by.y = "accession", all.x = T, all.y = F)
goa <- goa[!is.na(protein_id) & GO_ID %chin% go.dt$id, .(protein_id, GO_ID, Evidence_Code, Assigned_By)]
fwrite(goa, file = "data/go/goa_human.tsv", col.names = T, row.names = F, sep = "\t", quote = T, na = "None")
