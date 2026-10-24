Estrogen receptor beta programs associated with sex and bortezomib response in multiple myeloma
R analysis and figure-generation scripts accompanying the manuscript.

Overview
This study examines biological sex, ESR2 (encoding estrogen receptor beta, ERβ), and tumor–immune interactions in multiple myeloma (MM). It reanalyzes bone-marrow single-cell RNA-sequencing data from 20 healthy donors and 18 treatment-naïve MM patients sampled before bortezomib–melphalan–prednisone (VMP) therapy, alongside pharmacogenomic analyses and ERβ-overexpression experiments.

The manuscript reports sex-dependent ESR2 expression, communication networks and co-expression programs. ESR2 correlated with bortezomib resistance in female-derived cell lines (r = 0.806), but not male-derived lines (r = 0.038). ERβ overexpression increased drug-response AUC by 42.7% in OCI-MY5 and 17.3% in ARP-1. These are reported study findings, not results recalculated for this README.

Scripts

Figure1.R	 Cell annotation, marker expression, cell composition and sample-level PCA
Figure2.R	 ESR2 expression, response-group comparisons and sex-stratified ROC curves
Figure3.R	 Malignant plasma-cell subclusters and ESR2-positive cell proportions
Figure4.R	 Saved CytoTRACE scores, Monocle 3 trajectories and ESR2 along pseudotime
Figure5.R	 Sex- and response-stratified CellChat communication networks
Figure6.R	 Saved hdWGCNA modules, TOM networks, KEGG enrichment, GSVA and GSEA
Figure7.R	 Panels A–D: workflow, OR-labelled cell proportions and CTRP drug-response correlations

#R version （R 4.4.3）
