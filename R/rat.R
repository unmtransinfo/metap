#!/usr/bin/env Rscript

library(data.table)
library(tidyr)

gene <- read.delim2("ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt", header = T, comment.char = "#", stringsAsFactors = F, na.strings = "", quote = "\"")
setDT(gene)
gene <- gene[!is.na(UNIPROT_ID), .(GENE_RGD_ID, UNIPROT_ID)]
gene <- separate_rows(gene, UNIPROT_ID, sep = ";")

library(RPostgreSQL)
conn <- dbConnect(PostgreSQL(), user = "oleg", host = "localhost", dbname = "metap")
protein <- dbGetQuery(conn, "select protein_id,accession from protein where accession is not null")
dbDisconnect(conn)
rm(conn)
setDT(protein)

gene <- merge(gene, protein, by.x = "UNIPROT_ID", by.y = "accession")
fwrite(gene[, .(protein_id, GENE_RGD_ID)], file = "data/rgd/protein2rgd.tsv", sep = "\t", col.names = T, row.names = F, na = "None", quote = T)
if(file.exists("data/rgd/protein2rgd.tsv.gz")) {
  file.remove("data/rgd/protein2rgd.tsv.gz")
}
system("gzip -9v data/rgd/protein2rgd.tsv")


qtl <- fread("ftp://ftp.rgd.mcw.edu/pub/data_release/QTLS_RAT.txt", header = T, sep = "\t", na.strings = "", quote = "\"", skip = 70, verbose = T, col.names = c("QTL_RGD_ID","SPECIES","QTL_SYMBOL","QTL_NAME","CHROMOSOME_FROM_REF","LOD","P_VALUE","VARIANCE","FLANK_1_RGD_ID","FLANK_1_SYMBOL","FLANK_2_RGD_ID","FLANK_2_SYMBOL","PEAK_RGD_ID","PEAK_MARKER_SYMBOL","TRAIT_NAME","MEASUREMENT_TYPE","(UNUSED)","PHENOTYPES","ASSOCIATED_DISEASES","CURATED_REF_RGD_ID","CURATED_REF_PUBMED_ID","CANDIDATE_GENE_RGD_IDS","CANDIDATE_GENE_SYMBOLS","INHERITANCE_TYPE","RELATED_QTLS","UNUSED","5.0_MAP_POS_CHR","5.0_MAP_POS_START","5.0_MAP_POS_STOP","5.0_MAP_POS_METHOD","3.4_MAP_POS_CHR","3.4_MAP_POS_START","3.4_MAP_POS_STOP","3.4_MAP_POS_METHOD","CROSS_TYPE","CROSS_PAIR","STRAIN_RGD_ID1","STRAIN_RGD_ID2","STRAIN_RGD_SYMBOL1","STRAIN_RGD_SYMBOL2","6.0_MAP_POS_CHR","6.0_MAP_POS_START","6.0_MAP_POS_STOP","6.0_MAP_POS_METHOD","STRAIN_RGD_ID3","STRAIN_RGD_SYMBOL3","SSTRAIN"))
qtl <- qtl[!is.na(CANDIDATE_GENE_RGD_IDS)]
qtl <- separate_rows(qtl, CANDIDATE_GENE_RGD_IDS, CANDIDATE_GENE_SYMBOLS, sep = ";", convert = T)
qtl <- qtl[, .(QTL_RGD_ID, QTL_SYMBOL, QTL_NAME, LOD, P_VALUE, TRAIT_NAME, MEASUREMENT_TYPE, ASSOCIATED_DISEASES, CANDIDATE_GENE_RGD_IDS, PHENOTYPES)]
qtl <- separate_rows(qtl, PHENOTYPES, sep = ";")
qtl <- qtl[CANDIDATE_GENE_RGD_IDS %in% gene$GENE_RGD_ID]
fwrite(qtl, file = "data/rgd/rat_qtl.tsv", quote = T, sep = "\t", na = "None", col.names = T, row.names = F)
if(file.exists("data/rgd/rat_qtl.tsv.gz")) {
  file.remove("data/rgd/rat_qtl.tsv.gz")
}
system("gzip -9v data/rgd/rat_qtl.tsv")

if(file.exists("data/rgd/rat_terms.tsv")) {
  file.remove("data/rgd/rat_terms.tsv")
}

rat.do <- fread("ftp://ftp.rgd.mcw.edu/pub/data_release/with_terms/rattus_terms_do", sep = "\t", na.strings = "", skip = 27, verbose = T, quote = "")
rat.do <- rat.do[OBJECT_TYPE == "gene"]
rat.do <- rat.do[, .(RGD_ID, OBJECT_SYMBOL, TERM_ACC_ID, TERM_NAME, QUALIFIER, EVIDENCE)]
rat.do[, ONTOLOGY := "Disease Ontology"]
rat.do <- rat.do[RGD_ID %in% gene$GENE_RGD_ID]
rat.do <- unique(rat.do, by = c("RGD_ID", "TERM_ACC_ID"))
fwrite(rat.do, file = "data/rgd/rat_terms.tsv", append = file.exists("data/rgd/rat_terms.tsv"), col.names = !file.exists("data/rgd/rat_terms.tsv"), sep = "\t", row.names = F, quote = T, na = "None")

rat.mp <- fread("ftp://ftp.rgd.mcw.edu/pub/data_release/with_terms/rattus_terms_mp", sep = "\t", na.strings = "", skip = 27, verbose = T, quote = "")
rat.mp <- rat.mp[OBJECT_TYPE == "gene"]
rat.mp <- rat.mp[, .(RGD_ID, OBJECT_SYMBOL, TERM_ACC_ID, TERM_NAME, QUALIFIER, EVIDENCE)]
rat.mp[, ONTOLOGY := "Mammalian Phenotype"]
rat.mp <- rat.mp[RGD_ID %in% gene$GENE_RGD_ID]
rat.mp <- unique(rat.mp, by = c("RGD_ID", "TERM_ACC_ID"))
fwrite(rat.mp, file = "data/rgd/rat_terms.tsv", append = file.exists("data/rgd/rat_terms.tsv"), col.names = !file.exists("data/rgd/rat_terms.tsv"), sep = "\t", row.names = F, quote = T, na = "None")

rat.rdo <- fread("ftp://ftp.rgd.mcw.edu/pub/data_release/with_terms/rattus_terms_rdo", sep = "\t", na.strings = "", skip = 27, verbose = T, quote = "")
rat.rdo <- rat.rdo[OBJECT_TYPE == "gene"]
rat.rdo <- rat.rdo[, .(RGD_ID, OBJECT_SYMBOL, TERM_ACC_ID, TERM_NAME, QUALIFIER, EVIDENCE)]
rat.rdo[, ONTOLOGY := "RGD Disease Ontology"]
rat.rdo <- rat.rdo[RGD_ID %in% gene$GENE_RGD_ID]
rat.rdo <- unique(rat.rdo, by = c("RGD_ID", "TERM_ACC_ID"))
fwrite(rat.rdo, file = "data/rgd/rat_terms.tsv", append = file.exists("data/rgd/rat_terms.tsv"), col.names = !file.exists("data/rgd/rat_terms.tsv"), sep = "\t", row.names = F, quote = T, na = "None")

if(file.exists("data/rgd/rat_terms.tsv.gz")) {
  file.remove("data/rgd/rat_terms.tsv.gz")
}
system("gzip -9v data/rgd/rat_terms.tsv")
