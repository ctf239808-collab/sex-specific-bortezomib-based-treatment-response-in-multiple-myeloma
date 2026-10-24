
Sys.setenv(LANGUAGE = "en")
# #Do not coerce to factor
options(stringsAsFactors = FALSE)
# #Clear workspace
rm(list=ls())


# ==== load ====
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpubr)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)
library(rmarkdown)
library(tibble)
library(future)
library(furrr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

setwd("/Users/tf/e/ndmm/HM38")

scRNA_harmony=readRDS("scRNA_harmony_HM38.anno.rds")

celltype_colors_set3_muted_hd1 <- c(
  
  ## =========================
  # # B lineage ( → ，)
  ## =========================
  "B cells"              = "#8FAFC9", #
  "Pro-B/Pre-B cells"    = "#6F93B3", #
  "LymP"                 = "#B3AFCF", #
  
  ## =========================
  # # Plasma cells (，)
  ## =========================
  "Malignant plasma cells"  = "#B22222", # ，
  
  ## =========================
  # # T / NK lineage ( → )
  ## =========================
  "CD4+ T cells"         = "#B7D6C2", #
  "CD8+ T cells"         = "#86BFA3", #
  "NK cells"             = "#7FBFC6", # ( T )
  
  ## =========================
  # # Monocytes / Macrophage ( → → )
  ## =========================
  "CD14+ monocytes"          = "#E1B65C", #
  "CD16+ Monocytes"          = "#EAD28F", #
  "Intermediate monocytes"   = "#B3DE69", #
  "Macrophages"              = "#DFAF9A", # ()
  
  ## =========================
  # # DC lineage ( → )
  ## =========================
  "pDCs"                 = "#D3BFE2", # 、
  
  ## =========================
  # # Progenitors / Cycling ( vs )
  ## =========================
  "HSPCs"                = "#BFA65A", #
  "GMPs"                 = "#AFAFAF", #
  "ERPs"                 = "#E6E6E6", #
  
  ## =========================
  # # Other BM cells ()
  ## =========================
  "Megakaryocytes"       = "#F768A1", #
  "Erythrocytes"         = "#D77486" # ( plasma )
)

# === Fig.1B ===

DimPlot(scRNA_harmony,group.by = "anno",label=F,raster=FALSE, cols = celltype_colors_set3_muted_hd1)+
  ggtitle(NULL)+
  theme(legend.title = element_blank()) +  # title
  theme(
    # axis
    axis.line = element_line(linewidth = 1),  # axis
    # axislabel (xy)
    axis.title = element_text(face = "bold", size = 12),
    # axisticklabel (tick)
    axis.text = element_text(face = "bold", size = 12),
    # legendtitle
    legend.title = element_text(face = "bold", size = 12),
    # legendlabel
    legend.text = element_text(face = "bold", size = 12)
  )

# === Fig.1C ===
markers.anno = c("CD34","AVP",#HSPCs,2
                 "AZU1","MPO",#GMPs,2
                 "IGLL1","VPREB1","DNTT",#LymP,3
                 "MME","BCL7A","VPREB3",#"DNTT",Pro-B,3
                 "TCL1A","CD79B",
                 "IGHD","IL4R",#Navie B,2
                 "CD19","CD79A","MS4A1",# B("CD79B",),4
                 "CD83","NFKB1",#Activated NFkB,3
                 "CD27",#CD27+ memory,"LTB",2
                 "FCRL5",
                 "IGHG1","IGHG2","IGHG3","MZB1","SDC1","JCHAIN", # Plasma,6
                 "IGHA1","IGHA2",#Plasmablast,2
                 "S100A9","CD14",#CD14+mono,2
                 "IL7R","LTB", #T_CD4
                 "CD3D","CD3E","CD8A",##T_CD8
                 "GNLY","NKG7",#NK,2
                 "NCAM1",##CD56brightCD16lo
                 "FCGR3A","MS4A7",#CD16+mono,2
                 "G0S2","NAMPT",# CD14++CD16+intermediate monocytes
                 "LILRA4","PLD4",##pDCs,2
                 "FCER1A","CD1C","CLEC10A",#mDCs(myeloid dendritic cells),3
                 "PRSS57",#MPP,1
                 "RNASE2",#Granulocyte,1
                 "CD68", "CD163","C1QA","C1QB","FOLR2",# macrophages(Mφ),5
                 "PPBP","PF4","GNG11",#Megakaryocyte巨核,"ITGA2B",3
                 "ITGA2B","APOC1","KLF1","TFR2",#erythroid progenitor,7
                 "NFIA","HBG2","AHSP","GYPA","SLC4A1",#Erythrocytes,9
                 "RRM2","TRIP13","NEK2",#CIN,3
                 "TOP2A","MKI67")#megakaryocyte progenitor,2
Idents(scRNA_harmony)="anno"
p <-DotPlot(scRNA_harmony, features = markers.anno)+ 
  theme(axis.text.x.bottom = element_text(angle = 90,hjust = 0.5,vjust = 0.5))
dot.data=p[["data"]]
colnames(dot.data)
p01=ggplot(dot.data,aes(y=id, x = features.plot, color = avg.exp.scaled, size = pct.exp)) + 
  geom_point() +theme(axis.text.x.bottom = element_text(angle = 90,hjust = 1,vjust = 1))+
  cowplot::theme_cowplot()+
  scale_color_viridis_c(name = 'expression') + 
  # theme(legend.position="top")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +scale_color_gradientn(colours = viridis::plasma(20), limits = c(0,3), oob = scales::squish, name = 'expression')+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1) 
ggsave("HM38.anno.dotplot.markers01.tiff", dpi=300, plot = p01, width = 8, height = 4)


# === Fig.1D ===

sc.combined=scRNA_harmony
Idents(sc.combined)="group"
table(sc.combined$group)# group.asamplegrouping
table(sc.combined$anno)# annocell type
pB2_df <- table(sc.combined@meta.data$anno,sc.combined@meta.data$group) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
cluster=c("B cells","CD14+ monocytes","CD16+ Monocytes",
          "CD4+ T cells","CD8+ T cells","ERPs","Erythrocytes","GMPs",
          "HSPCs","Intermediate monocytes","LymP","Malignant plasma cells",
          "mDCs","Macrophages","Megakaryocytes", "Naïve B cells",
          "NK cells","pDCs","Plasma cells","Pro-B/Pre-B cells")
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = cluster)
pB2_df$Sample <- factor(pB2_df$Sample,levels = c("HD","MM"))

table(sc.combined@meta.data$group)

#Fig.1D(left)
ggplot(pB2_df, aes(y = Number, x = Sample, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = celltype_colors_set3_muted_hd1) +
  theme_classic() +
  labs(y = "Proportion", x = "") +
  guides(fill = guide_legend(title = "")) +
  theme(
    legend.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1,face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  )


#Fig.1D(right)
Idents(sc.combined)="group.a"
sc.combined=subset(sc.combined,ident=c("OR","SOR"))
table(sc.combined$group.a)# group.asamplegrouping
table(sc.combined$anno)# annocell type
pB2_df <- table(sc.combined@meta.data$anno,sc.combined@meta.data$group.a) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
cluster=c("B cells","CD14+ monocytes","CD16+ Monocytes",
          "CD4+ T cells","CD8+ T cells","ERPs","Erythrocytes","GMPs",
          "HSPCs","Intermediate monocytes","LymP","Malignant plasma cells",
          "Macrophages","Megakaryocytes", "NK cells","pDCs","Pro-B/Pre-B cells")
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = cluster)
pB2_df$Sample <- factor(pB2_df$Sample,levels = c("OR","SOR"))
table(sc.combined@meta.data$group.a)

ggplot(pB2_df, aes(x = Sample, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = celltype_colors_set3_muted_hd1) +
  theme_classic() +
  labs(y = "Proportion", x = "") +
  guides(fill = guide_legend(title = "")) +
  theme(
    legend.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1,face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  )

# === Fig.1F ===

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

meta <- scRNA_harmony@meta.data

cell_comp <- meta %>%
  dplyr::count(orig.ident, anno, name = "n_cells") %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(freq = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

cell_comp_mat <- cell_comp %>%
  dplyr::select(orig.ident, anno, freq) %>%
  tidyr::pivot_wider(
    names_from  = anno,
    values_from = freq,
    values_fill = 0
  )

pca_res <- prcomp(
  cell_comp_mat %>% dplyr::select(-orig.ident),
  scale. = TRUE
)

pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$orig.ident <- cell_comp_mat$orig.ident

sample_info <- meta %>%
  dplyr::select(orig.ident, group, group.a) %>%
  dplyr::distinct()

sample_info$group    <- trimws(as.character(sample_info$group))
sample_info$group.a <- trimws(as.character(sample_info$group.a))

pca_df <- dplyr::left_join(
  pca_df,
  sample_info,
  by = "orig.ident"
)

stopifnot(!all(is.na(pca_df$group)))

grp <- tolower(pca_df$group)
rsp <- tolower(pca_df$group.a)

pca_df$PCA_group <- dplyr::case_when(
  grp == "hd" ~ "HD",
  grp == "mm" & rsp == "or"  ~ "MM_OR",
  grp == "mm" & rsp == "sor" ~ "MM_SOR",
  grp == "mm" & is.na(rsp)   ~ "MM_Unknown",
  TRUE ~ NA_character_
)

table(pca_df$PCA_group)

var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_exp[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_exp[2] * 100, 1), "%)")


ggplot(pca_df, aes(PC1, PC2, color = PCA_group)) +
  geom_point(size = 3) +
  stat_ellipse(type = "norm", linewidth = 1) +
  scale_color_manual(values = c(
    "HD"     = "#4E79A7",
    "MM_OR"  = "#59A14F",
    "MM_SOR" = "#E15759"
  )) +
  labs(
    x = pc1_lab,
    y = pc2_lab,
    color = "Group",
    title = "PCA of cell-type composition across samples"
  ) +
  theme_classic(base_size = 14)



# === Fig.1E ===

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggsignif)
  library(patchwork)
})

figure1e_input <- "/Users/tf/e/ndmm/HM38/scRNA_harmony.HM38.annopre.rds"
figure1e_annotation <- "/Users/tf/e/ndmm/HM38/annoMM.txt"
figure1e_out <- paste0(
  "/Users/tf/e/ndmm/ESR2/submit/NC/",
  "Figure1_code/outputs"
)

# Optional: provide figure1e_meta explicitly to reuse the ORIGINAL object's
# metadata. Do not substitute metadata from the re-clustered newer object.
# This does not change or subset the scRNA_harmony object used by other panels.
if (!exists("figure1e_meta", inherits = FALSE)) {
  suppressPackageStartupMessages(library(Seurat))
  figure1e_original <- readRDS(figure1e_input)
  figure1e_meta <- figure1e_original@meta.data
  rm(figure1e_original)
  invisible(gc())
}

figure1e_required <- c("orig.ident", "group.a", "seurat_clusters")
stopifnot(all(figure1e_required %in% names(figure1e_meta)))
figure1e_mm <- figure1e_meta %>%
  dplyr::select(dplyr::all_of(figure1e_required)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ trimws(as.character(.x)))) %>%
  dplyr::filter(group.a %in% c("OR", "SOR"))
if (nrow(figure1e_mm) == 0L ||
    anyNA(figure1e_mm) || any(figure1e_mm$orig.ident == "")) {
  stop("Missing or invalid MM sample metadata.")
}

figure1e_map <- read.delim(
  figure1e_annotation, colClasses = "character", check.names = FALSE
)
stopifnot(all(c("ClusterID", "celltype") %in% names(figure1e_map)))
stopifnot(!anyDuplicated(figure1e_map$ClusterID))
if (anyNA(figure1e_map$celltype) || any(figure1e_map$celltype == "")) {
  stop("Missing cell-type labels in annoMM.txt.")
}
figure1e_mm$celltype <- figure1e_map$celltype[
  match(figure1e_mm$seurat_clusters, figure1e_map$ClusterID)
]
if (anyNA(figure1e_mm$celltype)) stop("Unmapped original MM clusters.")

figure1e_sample_info <- figure1e_mm %>% dplyr::distinct(orig.ident, group.a)
if (anyDuplicated(figure1e_sample_info$orig.ident)) {
  stop("A patient has more than one OR/SOR group assignment.")
}
figure1e_n <- table(factor(figure1e_sample_info$group.a,
                           levels = c("OR", "SOR")))
if (!identical(as.integer(figure1e_n), c(10L, 8L))) {
  stop("The original Figure 1E requires 10 OR and 8 SOR patients.")
}

# Denominator: ALL retained cells in the same patient, including plasmablasts
# and every other cell type. Filter the two target cell types only afterward.
figure1e_targets <- c("Malignant plasma cells", "CD4+ T cells")
figure1e_totals <- figure1e_mm %>% dplyr::count(orig.ident, name = "total_cells")
figure1e_counts <- figure1e_mm %>%
  dplyr::filter(celltype %in% figure1e_targets) %>%
  dplyr::count(orig.ident, celltype, name = "n_cells")
figure1e_data <- tidyr::expand_grid(
  orig.ident = figure1e_sample_info$orig.ident,
  celltype = figure1e_targets
) %>%
  dplyr::left_join(figure1e_counts, by = c("orig.ident", "celltype")) %>%
  dplyr::mutate(n_cells = dplyr::coalesce(n_cells, 0L)) %>%
  dplyr::left_join(figure1e_totals, by = "orig.ident") %>%
  dplyr::left_join(figure1e_sample_info, by = "orig.ident") %>%
  dplyr::mutate(
    proportion = 100 * n_cells / total_cells,
    group.a = factor(group.a, levels = c("OR", "SOR")),
    celltype = factor(celltype, levels = figure1e_targets)
  ) %>%
  dplyr::arrange(celltype, group.a, orig.ident)

# One observation per patient; two-sided, unpaired Welch t-tests.
figure1e_stats <- dplyr::bind_rows(lapply(figure1e_targets, function(ct) {
  z <- dplyr::filter(figure1e_data, celltype == ct)
  test <- stats::t.test(proportion ~ group.a, data = z, var.equal = FALSE,
                        alternative = "two.sided")
  data.frame(
    celltype = ct, n_OR = sum(z$group.a == "OR"),
    n_SOR = sum(z$group.a == "SOR"),
    mean_OR_percent = mean(z$proportion[z$group.a == "OR"]),
    mean_SOR_percent = mean(z$proportion[z$group.a == "SOR"]),
    t_statistic = unname(test$statistic), df = unname(test$parameter),
    p_value = test$p.value,
    p_label = formatC(test$p.value, format = "f", digits = 4),
    method = test$method
  )
}))

figure1e_plot <- function(ct) {
  z <- dplyr::filter(figure1e_data, celltype == ct)
  p_text <- figure1e_stats$p_label[figure1e_stats$celltype == ct]
  y_range <- diff(range(z$proportion))
  ggplot(z, aes(group.a, proportion, fill = group.a)) +
    geom_boxplot(outlier.shape = NA, width = 0.75, linewidth = 0.5) +
    geom_point(
      position = position_jitter(width = 0.1, height = 0, seed = 20261013),
      colour = "black", size = 1.5
    ) +
    scale_fill_manual(values = c("OR" = "#EF8A62", "SOR" = "#F7F7F7")) +
    ggsignif::geom_signif(
      comparisons = list(c("OR", "SOR")),
      annotations = p_text,
      y_position = max(z$proportion) + 0.05 * y_range,
      margin_top = 0, tip_length = 0.03, textsize = 4, size = 0.65
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.1))) +
    labs(title = ct, x = NULL, y = paste0("% ", ct, " in sample")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10, colour = "black", face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.line = element_line(colour = "black", linewidth = 0.6),
      legend.position = "none", plot.margin = margin(8, 8, 8, 8)
    )
}

p1E_malignant <- figure1e_plot("Malignant plasma cells")
p1E_CD4 <- figure1e_plot("CD4+ T cells")
p1E <- p1E_malignant / p1E_CD4
p1E
