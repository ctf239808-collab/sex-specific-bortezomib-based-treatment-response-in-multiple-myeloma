
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





