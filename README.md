# MetaPathML codebase

Copied from Oleg Ursu's directory `metap` on seaborgium, minus data.

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
  * **Other:**
    * R packages may also have dependencies, additional R packages,
      and/or libraries for compilation and/or execution.
  * See table of selected R package dependencies below.

---
### Metapath Semantics:

A.k.a. network topology of evidence chains.

From implementation `gene_disease_omim.R` (by Oleg Ursu)
(See also `gene_disease_mp.R` for MPO-based workflow.)

  * All paths FROM genes (protein-coding) (UniProt ID).
  * All paths mapped TO OMIM disease IDs.
  * All known direct gene-disease associations are from OMIM, ClinVar, and RGD.
  * Since negative labels unavailable, use random selection of genes NOT known associated with disease/phenotype.
  * Starting ("left-side") genes NOT known associated with disease/phenotype.


| Source | Semantics |
| :--- | :--- |
| OMIM | Gene to OMIM disease, via OMIM phenotype and phenotype series.  |
| ClinVar | Gene to OMIM disease, via OMIM phenotype and phenotype series.  |
| RGD | Rat gene to OMIM disease via RGD mappings to OMIM |
| STRINGDB | Gene to gene via PPI to OMIM disease, weighted by STRING score and connectivity counts.  |
| Reactome | Gene to gene via Reactome pathways, to OMIM disease via known gene to OMIM disease see sources above |
| KEGG | Gene to gene via KEGG pathways |
| LINCS1000 | Drug gene perturbation profiles associate similar genes |
| GeneOntology (GO) | Gene to gene via GO annotations |
| InterPro | Gene to gene via InterPro annotations |
| HumanProteinAtlas (HPA) | Gene to tissue, to disease via localization (partial, manual).  |
| CCLE | NOT METAPATHS: gene expression signature simple ML features combined with vectorized metapaths.  |
| GTEx | NOT METAPATHS: gene expression signature simple ML features combined with vectorized metapaths.  |

---
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
### Progress Summary

In 2018, Oleg Ursu at UNM Translational Informatics Division developed a supervised
machine learning workflow designed to predict and rank disease-gene associations,
given a disease query, using a metapath-based knowledge-graph embedding approach,
coupled with the powerful XGBoost package from Tianqi Chen et al. This project built
upon prior metapath methodology including by Chen, Fu, and Himmelstein, most closely
resembling the Himmelstein method, but with novel modifications based on (1)
biomedical domain expertise, (2) software engineering, and (3) interpretability
considerations.  In short, this project endeavors to combine cutting edge ML tools
with deep understanding of the knowledge domain, and produce both a robust
application, amenable to reuse and enterprise deployment, and a transparent
methodology and workflow, amenable to scientific interpretation and improvement
through a virtuous cycle of validation and evaluation against experimental or other
gold standard ground truth.

As of October 2018, twelve datasets have been integrated via reusable pipelines
coded in R, Python and Java. The R makes extensive use of the "data.table" package.
Each dataset requires careful mapping of genes, proteins, diseases, phenotypes, and
other terms and IDs. A custom PostgreSql database is used for integration and
performance. The XGBoost input matrices are generated in R, with one row per gene,
and one column for each instance of an association for each dataset (e.g. gene to GO
term). In this way the results and variable importance scores are interpretable with
provenance to source data. Training labels, a.k.a. known associations, are from
selected high-confidence sources (OMIM, ClinVar, RGD). XGBoost access is via
supported R package. Each model is defined by a disease or phenotype or phenotype
series for which a sufficient known gene associations exist (default cutoff = 50)
Five-way cross validation is employed on 80% (?) of the dataset to optimize
parameters for the final model training. For this project, a server "seaborgium" was
purchased, with RAM 0.5TB and 44 cores, necessary to process large matrices
efficiently. The training for a typical model requires 24h or compute time.  Thus
far approximately ten models have been generated and are currently under review.
Scientific colleagues have been engaged to assist in usability and validation
efforts.

### Next Steps and Work Remaining

From the ML perspective the current system is fully functioning. For each model,
output files include performance measures, ROC curves, and variable importance
scores. However, this is not sufficient for interpretation and contextualization by
domain scientists. An interactive user interface is needed, to visualize and
interrogate the model, and integrate with relevant data such as via enrichment and
biochemical pathway tools. High ranking gene hits should be accompanied by names,
IDs, and link-outs to standard references.  Other essential next steps relate to
external validation. The gold standard for validation would be prospective
experimental confirmation of generated hypotheses. Another form of validation could
be comparison with existing methods, such as the Himmelstein method. It should be
noted that from a learning perspective in-validation also provides great value. It
is very likely that some models will perform much better than others, and there may
be patterns to be learned, informing which diseases and phenotypes and pre-requisite
evidence will yield better models. Another important step relates to integration
with Pharos/TCRD. Evidence criteria and presentation content need development and
testing. Finally, the components need to be engineered for automated updates,
enterprise deployment, and convenient generation of models by scientist users. In
software engineering terms this can be described as developing the prototype system
into a professional product.

Also important is organizing and documenting the code for sharing, collaboration,
and reuse. As of Oct 2018 the development, deployment and all use were on
one computer at UNM, via one login.

---

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
