# MetaPathML codebase

 * [Semantics](#semantics)
 * [Database](#database)
 * [Dependencies](#dependencies)

See also: [MetaPathML History and Progress Notes](doc/HISTORY.md)

---
<a name="semantics" />
### Metapath Semantics:

A.k.a. network topology of evidence chains.

From implementation `gene_disease_omim.R` (by Oleg Ursu)
(See also `gene_disease_mp.R` for MPO-based workflow.)

  * All paths FROM genes (protein-coding) (UniProt ID).
  * All paths mapped TO OMIM disease IDs.
  * All known direct gene-disease associations are from OMIM, ClinVar, and RGD.
  * Since negative labels unavailable, use random selection of genes NOT known associated with disease/phenotype.
  * Starting ("left-side") genes NOT known associated with disease/phenotype.


| Group | Source | Semantics |
| :--- | :--- | :--- |
| Disease Association | OMIM | Gene to OMIM disease/phenotype/phenotype-series |
| Disease Association | ClinVar | Gene to OMIM disease/phenotype  |
| Disease Association | RGD | Rat gene to OMIM disease |
| Disease Association | IMPC | Mouse gene to OMIM disease |
| PPI | STRINGDB | Gene to gene via PPI to OMIM disease, weighted by STRING score and connectivity counts.  |
| Pathways | Reactome | Gene to gene via Reactome pathways, to OMIM disease |
| Pathways | KEGG | Gene to gene via KEGG pathways |
| Expression | LINCS | Drug gene perturbation profiles associate similar genes |
| Expression | HumanProteinAtlas (HPA) | Gene to tissue, to disease via localization (partial, manual).  |
| Expression | GTEx | NOT METAPATHS: expression profile fused with vectorized metapaths.  |
| Expression | CCLE | NOT METAPATHS: expression profile fused with vectorized metapaths.  |
| Function | GeneOntology (GO) | Gene to gene via GO annotations |
| Function | InterPro | Gene to gene via InterPro annotations |
| Ontology | DO |  Disease Ontology | 
| Ontology | MP |  Mammalian Phenotypes | 
| IDs | UniProt |  Protein IDs | 

---
<a name="database" />
### Database

This is the PostgreSql db built from various selected sources and used at runtime
for model building. A knowledge graph is represented via relational db. A separate
`metap_dev` also exists.

`$ echo '\d+' |psql -d metap`
```
                               List of relations
 Schema |          Name          |   Type   | Owner |    Size    |  
--------+------------------------+----------+-------+------------+
 public | ccle                   | table    | oleg  | 1237 MB    |
 public | ccle_ccle_id_seq       | sequence | oleg  | 8192 bytes |
 public | clinvar                | table    | oleg  | 189 MB     |
 public | clinvar_cv_id_seq      | sequence | oleg  | 8192 bytes |
 public | clinvar_disease        | table    | oleg  | 728 kB     |
 public | clinvar_disease_xref   | table    | oleg  | 960 kB     |
 public | disease_onto           | table    | oleg  | 696 kB     |
 public | do_xref                | table    | oleg  | 2432 kB    |
 public | drug_name              | table    | oleg  | 256 kB     |
 public | ensembl                | table    | oleg  | 5288 kB    |
 public | gene_info              | table    | oleg  | 912 kB     |
 public | go                     | table    | oleg  | 13 MB      |
 public | goa                    | table    | oleg  | 25 MB      |
 public | gtex                   | table    | oleg  | 188 MB     |
 public | gwas                   | table    | oleg  | 17 MB      |
 public | homology               | table    | oleg  | 3096 kB    |
 public | homology_hid_seq       | sequence | oleg  | 8192 bytes |
 public | hpa_norm_tissue        | table    | oleg  | 65 MB      |
 public | interpro               | table    | oleg  | 2680 kB    |
 public | interproa              | table    | oleg  | 20 MB      |
 public | iuphar_class           | table    | oleg  | 240 kB     |
 public | jensen_disease         | table    | oleg  | 8560 kB    |
 public | kegg_pathway           | table    | oleg  | 1920 kB    |
 public | lincs                  | table    | oleg  | 13 GB      |
 public | lincs_lincs_id_seq     | sequence | oleg  | 8192 bytes |
 public | mouse2human            | view     | oleg  | 0 bytes    |
 public | mousephenotype         | table    | oleg  | 202 MB     |
 public | mp_onto                | table    | oleg  | 1048 kB    |
 public | ncbi                   | table    | oleg  | 2640 kB    |
 public | omim                   | table    | oleg  | 1952 kB    |
 public | omim2phen_series       | table    | oleg  | 168 kB     |
 public | omim_phen_series       | table    | oleg  | 56 kB      |
 public | protein                | table    | oleg  | 15 MB      |
 public | protein2mgi            | table    | oleg  | 3352 kB    |
 public | protein2rgd            | table    | oleg  | 1064 kB    |
 public | protein_protein_id_seq | sequence | oleg  | 8192 bytes |
 public | rat2human              | table    | oleg  | 672 kB     |
 public | rat_qtl                | table    | oleg  | 512 kB     |
 public | rat_term               | table    | oleg  | 12 MB      |
 public | rdo                    | table    | oleg  | 1264 kB    |
 public | rdo_xref               | table    | oleg  | 1480 kB    |
 public | reactome               | table    | oleg  | 2248 kB    |
 public | reactomea              | table    | oleg  | 15 MB      |
 public | stringdb               | table    | oleg  | 1056 kB    |
 public | stringdb_score         | table    | oleg  | 215 MB     |
 public | uniprot_disease        | table    | oleg  | 488 kB     |
(46 rows)
```
---
<a name="dependencies" />
### Dependencies:

  * **R:**
    * Version 3.4.3 (Nov 2017)
    * Other versions may be ok, depending on package requirements.
  * **R Packages:**
    * **ML/Modeling:**
      ADPclust, xgboost, xgboostExplainer, randomForest, randomForestSRC, ggRandomForests,
      caret, DMwR (Data Mining with R), cluster, stringdist, pROC (Display and Analyze ROC Curves)
    * **Plotting:**
      ggplot2, plotly, plotROC, waterfalls
    * **R-enhancement:**
      data.table, Matrix, tidyr, tidyverse, fst
    * **I/O/Db/Misc:**
      readr, readxl, RJDBC, RPostgreSQL, KEGGREST, ontologyIndex
    * See table of selected R package dependencies below.
  * **Other:**
    * R packages may also have dependencies, additional R packages,
      and/or libraries for compilation and/or execution.
  * **PostgreSQL:**
    * Version 10.5 (for improved parallelization)
    * Versions 9.* should work with reduced performance.
    * Other versions may be ok.

### R package dependencies (selected)

From `R -e 'write.csv(data.table(installed.packages()))'`

| Package | Version | Depends | Imports | NeedsCompilation |
| :--- | :---: | :--- | :--- | :---: |
| data.table | 1.11.4 | R (>= 3.1.0) | methods | yes |
| DMwR | 0.4.1 | R(>= 2.10), methods, graphics, lattice (>= 0.18-3), grid (>= 2.10.1) | xts (>= 0.6-7), quantmod (>= 0.3-8), zoo (>= 1.6-4), abind (>= 1.1-0), rpart (>= 3.1-46), class (>= 7.3-1), ROCR (>= 1.0) | no |
| fst | 0.8.8 | R (>= 3.0.0) | Rcpp | yes |
| knitr | 1.2 | R (>= 3.1.0) | evaluate (>= 0.10), highr, markdown, stringr (>= 0.6), yaml, methods, tools | no |
| parallelDist | 0.2.1 | R (>= 3.0.2) | Rcpp (>= 0.12.6), RcppParallel (>= 4.3.20) | yes |
| readr | 1.1.1 | R (>= 3.0.2) | Rcpp (>= 0.12.0.5), tibble, hms, R6 | yes |
| recipes | 0.1.3 | R (>= 3.1), dplyr, broom | tibble, stats, ipred, dimRed (>= 0.1.0), lubridate, timeDate, ddalpha, purrr (>= 0.2.3), rlang (>= 0.1.1), gower, RcppRoll, tidyselect (>= 0.1.1), magrittr, Matrix, tidyr, pls | no |
| reshape2 | 1.4.3 | R (>= 3.1) | plyr (>= 1.8.1), Rcpp, stringr | yes |
| stringr | 1.3.1 | R (>= 3.1) | glue (>= 1.2.0), magrittr, stringi (>= 1.1.7) | no |
| tidyr | 0.8.1 | R (>= 3.1) | dplyr (>= 0.7.0), glue, magrittr, purrr, Rcpp, rlang, stringi, tibble, tidyselect | yes |
| tidyselect | 0.2.4 | R (>= 3.1) | glue, purrr, rlang (>= 0.2.0), Rcpp (>= 0.12.0) | yes |
| xgboost | 0.6.4.1 | R (>= 3.3.0) | Matrix (>= 1.1-0), methods, data.table (>= 1.9.6), magrittr (>= 1.5), stringi (>= 0.5.2) | yes |
| xgboostExplainer | 0.1 | R (>= 3.0.0) | data.table, xgboost, waterfalls, scales, ggplot2 | NA |
