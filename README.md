Estrogen receptor beta programs associated with sex and bortezomib response in multiple myeloma
R analysis and figure-generation scripts accompanying the manuscript.

Repository: sex-specific-bortezomib-based-treatment-response-in-multiple-myeloma

Study overview
This study examines how biological sex, expression of ESR2 (encoding estrogen receptor beta, ERβ), and tumor–immune interactions relate to heterogeneous bortezomib-based treatment responses in multiple myeloma (MM). The study combines a reanalysis of bone-marrow single-cell RNA-sequencing data from 20 healthy donors and 18 treatment-naïve MM patients, sampled before bortezomib–melphalan–prednisone (VMP) therapy, with pharmacogenomic analyses and ERβ-overexpression experiments in MM cell lines.

The seven R files cover cell annotation and composition, ESR2 expression, response-associated ROC curves, malignant plasma-cell subclusters, inferred trajectories, CellChat communication networks, hdWGCNA modules, pathway analyses, and archived cell-line drug-response correlations.

Manuscript abstract
The following manuscript abstract was provided for this README. Its numerical findings describe the overall study and were not independently re-estimated while preparing this documentation.

Treatment response heterogeneity in multiple myeloma (MM) remains incompletely understood regarding biological sex, ESR2 expression, and tumor-immune interactions. We reanalyzed bone-marrow single-cell RNA-sequencing data from 20 healthy donors and 18 treatment-naïve MM patients prior to bortezomib-melphalan-prednisone (VMP) therapy, integrating clinical annotations with sex-stratified trajectory, co-expression, and cell communication analyses. Pharmacogenomic correlations and ERβ overexpression were examined in two MM cell lines. ESR2 expression exhibited cell-type-, sex-, and response-dependent variation. Predictive AUCs for ESR2-only, immune-only, and combined models were 0.75, 0.42, and 0.42 in females versus 0.60, 0.80, and 0.70 in males. Trajectory and communication analyses revealed sex-divergent differentiation states and signaling networks. hdWGCNA identified female-specific complement-coagulation modules and male-specific proteasome/oxidative phosphorylation signatures; GSEA confirmed a sex-specific ESR2-regulated proteasome network. Pharmacogenomically, ESR2 positively correlated with bortezomib resistance in female-derived lines (r = 0.806, p = 0.00087) but not male-derived lines (r = 0.038, p = 0.929). ERβ overexpression reduced bortezomib sensitivity, increasing drug-response AUC by 42.7% (OCI-MY5) and 17.3% (ARP-1). ESR2 marks sex- and context-dependent transcriptional states and modulates bortezomib sensitivity in MM cells, linking sex as a biological variable to therapeutic heterogeneity.

Code coverage
The abstract describes the complete study. The current Figure7.R generates panels A–D only; the overexpression experiments and associated E/F panels are outside this script. The seven files also start from processed objects and tables: they do not provide an end-to-end workflow from raw sequencing reads or rebuild the original upstream hdWGCNA networks.

This documentation was checked on 14 October 2026, against commit ea4b254. All seven scripts passed R syntax parsing at that revision. Full analyses were not rerun for this documentation.

Script guide
Place the English and Chinese README files alongside Figure1.R–Figure7.R in the repository root.

Script	Main analyses	Outputs implemented in the script
Figure1.R	Annotated cell embedding, marker DotPlot, HD/MM and OR/SOR cell composition, sample-composition PCA, patient-level proportions of malignant plasma cells and CD4+ T cells	Explicitly saves HM38.anno.dotplot.markers01.tiff; other panels require interactive evaluation or explicit printing/saving; Figure 1E plot and statistical objects remain in the R session
Figure2.R	ESR2 UMAP/violin/box plots, expression and detection-rate summaries, malignant-cell OR/SOR comparisons, sex-stratified patient ROC curves	PDF/PNG for D/F/G/H/I; source tables, AUC CSV, ROC RDS, and sessionInfo.txt
Figure3.R	Nine malignant plasma-cell subclusters; ESR2 expression and positivity; patient-by-subcluster heatmap; pooled-cell proportion comparisons	Figure3_integrated and A–F panels; CSV source data, methods, session and export-device information
Figure4.R	Saved CytoTRACE scores, Monocle 3 trajectories on existing UMAP coordinates, ESR2 counts versus pseudotime	Figure4_integrated and A–M panels; trajectory cache, roots, source tables, methods and session information
Figure5.R	Sex-stratified CellChat networks; ESR2-stratified networks; malignant plasma-cell/immune interactions; within-OR/SOR sex comparisons	Figure5_rebuilt and A–G panels; network-summary caches, source tables and provenance
Figure6.R	Saved female/male hdWGCNA modules and eigengenes; TOM networks; gene UMAP; KEGG enrichment, GSVA and proteasome GSEA	Figure6_A-L, A–L panels, configuration RDS, source data, pathway caches and reconstruction notes
Figure7.R	Workflow diagram; OR-labelled malignant-cell proportions; archived CTRP ESR2–bortezomib AUC correlations	Figure7_AD, A–D panels and six D subpanels; counts, correlations, membership/exclusion records and reference checks
Figures 3–7 use configurable built-in graphics exports. Figure 3 defaults to PDF/PNG/SVG; Figures 4–7 default to PDF/PNG and also support SVG/TIFF. Figure 1 and Figure 2 do not automatically export every displayed panel.

Required data
The inspected repository contains the seven R scripts, but not the RDS, TXT, CSV, RDA, image or reference assets they read. Obtain the original processed inputs and update the paths before execution. Similar filenames do not guarantee compatible annotations or coordinates.

Single-cell objects and annotations
Input	Used by	Required content
scRNA_harmony_HM38.anno.rds	Figures 1, 2, 5	Annotated Seurat atlas with RNA expression; orig.ident, anno, group.a; additionally group for Figures 1/2 and gender for Figures 2/5. Figures 1/2 use saved embeddings; Figures 2/5 require log-normalized RNA data, including ESR2
scRNA_harmony.HM38.annopre.rds and annoMM.txt	Figure 1E	Original metadata with orig.ident, group.a, seurat_clusters; annotation table with unique ClusterID and celltype. The script requires 10 OR and 8 SOR patients
scRNA_harmony_MM18_MM.dim50res0.5.anno.rds	Figure 3D–F	Malignant-cell Seurat object with orig.ident, gender, group.a, anno3 and RNA counts containing ESR2
scRNA_harmony_MM18.MM.anno.rds	Figure 3A–C and Figure 4	Historical malignant-cell object with the corresponding metadata, saved umap, RNA counts and normalized expression. Figure 3 verifies matching cells, annotations and ESR2 counts across its two objects
CytoTRACE_plot_table_anno3.txt, CytoTRACE_plot_table_group.a.txt, CytoTRACE_plot_table_gender.txt	Figure 4	Cell IDs in the first column; CytoTRACE, Component1, Component2, Phenotype. Cell sets and phenotype assignments must match the Seurat metadata; CytoTRACE scores and Component1/2 coordinates must agree across the three tables
Metadata labels are case-sensitive:

orig.ident: sample/patient identifier, with consistent sex and response assignments.
group: HD or MM.
gender: Female or Male.
group.a: OR or SOR for response comparisons. These are the study's existing response categories; the scripts do not derive their clinical cutoffs.
anno: cell-type annotation, including the exact label Malignant plasma cells.
anno3: the nine labels MM, MM_ANXA1, MM_AZGP1, MM_CRIP1, MM_IGHGP, MM_IGLL5, MM_IL5RA, PBM, PBM_CCND1.
Additional inputs
Analysis	Required files or fields
Figure 2H/I, full ROC comparisons	Patient-level CSV with one row per patient and columns orig.ident, gender, group.a, ESR2_mean, mono_n, CD8_cyto; set figure2_roc_csv to its path
Figure 5E, default mode	Figure5E_original.jpg, with the MD5 expected by figure5_e_reference_md5. This reference image is not included in the inspected repository
Figure 6, saved networks	female.hdWGCNA_object.rds, Male.MM.hdWGCNA_object.rds, FMM.umap_df.csv, Male.MM.umap_df.csv, female.MM_TOM.rda, male.MM_TOM.rda
Figure 6, object structure	Metadata gender, group.a, orig.ident, ESR2_group; active hdWGCNA results, harmonized module eigengenes, cell UMAP and normalized RNA expression. Expected modules are FMM1–13 and Male.MM1–3
Figure 6, gene UMAP tables	First column used as row names; fields gene, module, color, hub, kME, UMAP1, UMAP2, consistent with saved modules/TOM
Figure 6, pathways	Local MSigDB KEGG LEGACY and KEGG caches, or permitted database downloads. Optional paths: msigdb_legacy_file, kegg_cache_file
Figure 7B	cell_metadata_esr2.csv.gz or uncompressed CSV with cell_id, sample, sex, disease, response_group, celltype, ESR2_count
Figure 7C/D	df_ctrp_unique_MM.csv, df_ctrp_unique_MM_female.csv, df_ctrp_unique_MM_male.csv, each with cell_line, ESR2, bortezomib
Figure 7 reference checks	Figure7B_pairwise_pvalues.csv with group1, group2, pstar, p_value_numeric; required to pass the default strict historical checks
Figure 6 TOM matrices must match the module gene order; unnamed historical matrices depend on the saved ordering. Figure 7 expects the archived cell-line memberships defined in the code, as well as concordant values across overall and sex-specific tables.

R environment and dependencies
The repository has no renv.lock or pinned package environment. Use an R installation compatible with the saved objects and selected package versions, and retain sessionInfo() for each run. Figures 3–7 require Cairo graphics support:

capabilities("cairo")
The following is a dependency-installation starting point, not a reconstruction of the original package versions. Some packages are loaded by the scripts even where no corresponding analysis is performed.

install.packages(c(
  "Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr",
  "tibble", "tidyverse", "patchwork", "ggpubr", "ggsignif", "ggrepel",
  "scales", "RColorBrewer", "circlize", "cowplot", "viridis",
  "plotly", "htmlwidgets", "rmarkdown", "future", "furrr",
  "pROC", "reshape2", "gridGraphics", "digest", "igraph", "jpeg",
  "binom", "gson", "aplot", "msigdbr", "BiocManager"
))

BiocManager::install(c(
  "EnhancedVolcano", "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "SingleCellExperiment", "SummarizedExperiment",
  "batchelor", "GSVA", "AnnotationDbi", "BiocParallel"
))
Use the upstream installation instructions for CellChat, Monocle 3, and hdWGCNA, including their system and package dependencies. Use a compatible CellChat implementation providing the CellChat namespace and CellChatDB.human used in Figure 5. Bioconductor setup is documented in the official installation guide.

Figure1.R calls unqualified melt() without attaching its provider. Load reshape2 before running it, or change those calls in your local copy to reshape2::melt().

Running the scripts
1. Download and prepare inputs
git clone https://github.com/ctf239808-collab/sex-specific-bortezomib-based-treatment-response-in-multiple-myeloma.git
cd sex-specific-bortezomib-based-treatment-response-in-multiple-myeloma
Obtain the required processed inputs separately. There is no mandatory Figure 1 → Figure 7 execution chain: each file reads its own inputs. Figures 4–6 may additionally reuse or create analysis/database caches.

2. Configure paths and modes
Replace author-specific paths beginning with /Users/tf/, including input, output, cache and reference-image paths. Figures 1/2 also call setwd(). Main settings are near the beginning of Figures 3–7; Figure 1E and the later Figure 2 block contain additional paths.

Script	Main settings to review
Figure 1	setwd, atlas input, figure1e_input, figure1e_annotation; figure1e_out is assigned but not used to save files
Figure 2	setwd, atlas input, figure2_out, figure2_roc_csv
Figure 3	figure3_analysis_file, figure3_abc_file, figure3_output_dir, export options
Figure 4	Input/table paths, output/cache paths, figure4_root_mode and root inputs
Figure 5	Input/output/cache paths, figure5_c_mode, figure5_e_mode, reference image, network/selection options
Figure 6	All paths and analysis options in figure6_config
Figure 7	Metadata/CTRP/reference paths, figure7_output_dir, figure7_strict_reference
Figures 1/2 clear the workspace with rm(list = ls()). Run them in a dedicated R session. For Figure 2 full ROC analyses, supplying the CSV is more reliable than preloading session objects that the opening command removes.

3. Parse, then run the selected figure
From the repository root, syntax can be checked without running analyses:

invisible(lapply(sprintf("Figure%d.R", 1:7), parse))
After configuring paths, use RStudio or an R session to execute the relevant file. For example:

library(reshape2)  # Figure 1 uses melt().
source("/absolute/path/to/repository/Figure1.R")
Absolute source paths avoid confusion when Figures 1/2 change the working directory. Bare plotting expressions may need to be run interactively or explicitly printed/saved; check the script guide for automatic exports.

Figures 4–7 allow function loading without automatic analysis by setting FIGURE4_RUN, FIGURE5_RUN, FIGURE6_RUN, or FIGURE7_RUN to "0". Figure 3 runs its builder when sourced and has no equivalent environment-variable switch. Example for Figure 7:

Sys.setenv(FIGURE7_RUN = "0")
source("/absolute/path/to/repository/Figure7.R")

figure7_result <- build_figure7(
  b_metadata_file = "/absolute/path/to/cell_metadata_esr2.csv.gz",
  ctrp_dir = "/absolute/path/to/archived_ctrp_tables",
  b_reference_file = "/absolute/path/to/Figure7B_pairwise_pvalues.csv",
  output_dir = "/absolute/path/to/results/Figure7",
  strict_reference = TRUE
)
For Figures 4–6, the corresponding entry points are build_figure4(...), figure5_main(), and build_figure6(figure6_config). After loading functions, update the relevant arguments/configuration or global settings before calling them. Set a FIGURE*_RUN variable back to "1" to restore automatic execution when sourcing.

Figure-specific methods and reproduction notes
Figures 1–3: expression, proportions and ROC curves
Figure 1E uses the original pre-annotation object and cluster map. Each observation is one patient; the proportion denominator includes all retained cells from that patient. It uses two-sided, unpaired Welch t-tests.

Figure 2H/I implement historical patient scores: ESR2 alone; microenvironment score z(mono_n) - z(CD8_cyto); and combined score z(ESR2_mean) + microenvironment score. Standardization is within sex after complete-case filtering; OR is coded as 1, SOR as 0, with automatic pROC direction. The script does not fit LASSO or perform an independent validation split. Original mono_n and CD8_cyto values are needed to reproduce all three curves. Without them, the script exports ESR2-only curves with _ESR2_only filenames and an incomplete feature template containing NA values.

Figure 3 defines ESR2 positivity as RNA raw count > 0. Panel D summarizes patient–subcluster percentages, with fewer than 20 cells greyed out by default. E/F use pooled-cell proportions, Wilson intervals, and Fisher tests with BH correction within each nine-subcluster comparison. Those tests, and the Figure 2F/G cell-level Welch tests, do not model dependence among cells from the same patient.

Figure 4: trajectory roots
CytoTRACE scores and coordinates are read from historical tables. The script does not recompute CytoTRACE. Monocle 3 uses the saved cell UMAP; original trajectory roots were not retained.

The current defaults specify root_mode = "interactive" and three approximate reference positions. These positions are newly selected roots, rather than recovered original IDs. To choose roots in an interactive R session, set figure4_root_positions <- NULL and leave the other root inputs NULL. Alternatively, supply root cells or principal-graph nodes, using only one root input. The optional "cytotrace" mode creates a new automatic ordering and requires all explicit root inputs to be NULL.

Record the selected roots and pseudotime coverage. H–M are reconstructed results; pseudotime represents inferred ordering of cell states, rather than elapsed treatment time.

Figure 5: pooled CellChat comparisons
The pipeline uses normalized RNA expression, CellChatDB.human, triMean, a default 100 bootstraps and support threshold 0.05. A–C/F and exploratory E pool atlas cells by sex, including HD/MM; G constructs sex-specific networks within OR and SOR. D splits cells using the global atlas ESR2 median.

C: default "legacy" preserves historical receiver grouping, with overlapping receiver bars that should not be read as total incoming strength; "corrected" groups by the actual communication partner. Record the selected mode.
E: default "original_reference" embeds the verified original JPG; this panel is not recomputed. "exploratory" computes a new pathway-probability summary that differs from the historical panel.
F/G: sex-associated log2 probability ratios are descriptive network comparisons. The displayed p values reflect within-network CellChat support, not a statistical test of a female–male difference.
Network summaries are cached against input hashes, cells, parameters and package versions. Changing inputs or settings may trigger recomputation.

Figure 6: saved hdWGCNA networks and pathways
Saved module assignments, harmonized eigengenes, TOM matrices and gene UMAP coordinates are required; the script does not infer the original networks anew. It preserves historical ESR2 High/Low labels and exports an audit against positive/zero normalized expression. C–F are pooled-cell descriptive plots; K/L use patient-level expression summaries.

J/K use MSigDB C2:CP:KEGG_LEGACY; the FoxO component of K and proteasome GSEA in L also use KEGG definitions. allow_pathway_download = TRUE is the actual current default, despite an older comment about explicit opt-in. For a frozen local run, supply caches and set it to FALSE. Database revisions can change results; preserve the database provenance and source tables. Proteasome GSEA stops if the requested pathway is not returned under the historical, hard-coded pvalueCutoff = 0.05.

Figure 7: response labels and archived cell-line groups
B reports the proportion of OR-labelled cells within each sex × ESR2 group. Its historical axis title, “Bortezomib Response Rate (%)”, must be interpreted in that cell-level context. Wilson intervals and six unadjusted Fisher tests do not turn these pooled-cell counts into patient response rates.

C/D use Pearson correlations and linear-regression confidence bands. The script retains archived group assignments and exclusions, including NCIH929's exclusion only from Female Low and MM1S's exclusion from fits because of missing AUC. These historical groups are not independently re-annotated by the script.

strict_reference = TRUE checks fixed historical B counts/tests and C/D statistics. Setting it to FALSE still computes and records these checks, but does not stop solely because a historical result differs. Hard-coded cell-line membership requirements remain enforced, and a supplied reference CSV must still have the required fields. A run with changed data or group definitions requires an explicit new-analysis adaptation.

Outputs, provenance and citation
Use a separate output directory for each figure. Figures 3–7 export individual panels and numerical source data; Figures 4–6 also use analysis or database caches. Preserve source tables, input manifests/hashes, method/provenance notes, configuration and session information alongside the final figures. Figure 1 requires additional explicit saves for most panels; Figure 2 saves D/F/G/H/I automatically.

The manuscript title above can be used to identify the associated study. Author list, journal, DOI and publication status were not supplied for this README; add the finalized bibliographic citation when available. Cite the underlying datasets and analysis packages according to their original documentation.

No license file was present in the inspected repository snapshot. This README does not assign a software or data license.
