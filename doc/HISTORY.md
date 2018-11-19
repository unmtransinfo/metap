# MetaPathML history and progress notes

Copied from Oleg Ursu's directory `metap` on seaborgium, minus data.

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

