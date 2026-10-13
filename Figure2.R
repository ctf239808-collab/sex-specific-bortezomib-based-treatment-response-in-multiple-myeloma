
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)       
library(ComplexHeatmap)
library(circlize)        
library(ggpubr)          
library(scales)          
library(ggrepel)         
library(enrichplot)      
library(clusterProfiler) 
library(org.Hs.eg.db)    

setwd("/Users/tf/e/ndmm/HM38")
scRNA_harmony=readRDS("scRNA_harmony_HM38.anno.rds")

# ===========================
# setparameters
# ===========================
esr2_gene <- "ESR2"
celltypes <- unique(scRNA_harmony$anno)
colors_celltype <- RColorBrewer::brewer.pal(n = length(celltypes), name = "Set3")

# ===========================
# Figure 2A. ESR2 UMAP
# ===========================
FeaturePlot(
  scRNA_harmony,
  features = esr2_gene,
  reduction = "umap",
  cols = c("#C9C9C9", "#6E140C"),
  pt.size = 0.5,
  order = TRUE,
  raster=FALSE
) +
  scale_color_gradientn(
    colours = c(
      "#C9C9C9",  
      "#FFFF56",  
      "#EE9438",  
      "#E73726", 
      "#6E140C"   
    ),
    values = c(0, 0.25, 0.50, 0.73, 1),
    limits = c(0, 3),
    breaks = 0:3,
    oob = scales::squish,
    name = "Expression\nLevel"
  ) +
  ggtitle("ESR2") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# ===========================
# Figure 2B. HD vs MM ESR2
# ===========================

Idents(scRNA_harmony)="group"
p1=FeaturePlot(
  subset(scRNA_harmony,ident="HD"),
  features = esr2_gene,
  reduction = "umap",
  cols = c("#C9C9C9", "#6E140C"),
  pt.size = 0.5,
  order = TRUE,
  raster=FALSE
) +
  scale_color_gradientn(
    colours = c(
      "#C9C9C9",  
      "#FFFF56",  
      "#EE9438",  
      "#E73726",  
      "#6E140C"   
    ),
    values = c(0, 0.25, 0.50, 0.73, 1),
    limits = c(0, 3),
    breaks = 0:3,
    oob = scales::squish,
    name = "Expression\nLevel"
  ) +
  ggtitle("HD") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

p2=FeaturePlot(
  subset(scRNA_harmony,ident="MM"),
  features = esr2_gene,
  reduction = "umap",
  cols = c("#C9C9C9", "#6E140C"),
  pt.size = 0.5,
  order = TRUE,
  raster=FALSE
) +
  scale_color_gradientn(
    colours = c(
      "#C9C9C9", 
      "#FFFF56", 
      "#EE9438", 
      "#E73726",  
      "#6E140C"   
    ),
    values = c(0, 0.25, 0.50, 0.73, 1),
    limits = c(0, 3),
    breaks = 0:3,
    oob = scales::squish,
    name = "Expression\nLevel"
  ) +
  ggtitle("MM") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

p <- (p1 + p2 +
        plot_layout(
          ncol = 2,          
          guides = "collect" 
        )) &
  theme(legend.position = "right")

p

# ===========================
# Figure 2C violin plot
# ===========================
VlnPlot(
  scRNA_harmony,
  features = esr2_gene,
  group.by = "anno",
  split.by = "group",
  pt.size = 0.0001,
  combine = TRUE,
  raster=FALSE
) + 
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ===========================
# Figure 2E. sexstratifiedboxplot (cell type)
# ===========================
table(scRNA_harmony$anno)
key_cells <- c("Malignant plasma cells", "NK cells", "CD4+ T cells","CD8+ T cells")
df <- FetchData(scRNA_harmony, vars = c(esr2_gene, "gender", "group", "anno")) %>%
  filter(anno %in% key_cells)

ggplot(df, aes(x = gender, y = ESR2, fill = gender)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  facet_grid(group ~ anno, scales = "free_y") +
  stat_compare_means(aes(group = gender), label = "p.signif") +
  theme_classic() +
  ggtitle("")


# Figure.2D, F, G, H and I
# Can be sourced after the user's existing A/B/C/E code.
# Expression: RNA/data (log-normalized); zero-expression cells are retained.
# F/G: exploratory cell-level Welch t-tests, matching the recovered code.
# H/I: patient-level historical scores, NOT fitted LASSO predictions.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pROC)
})

figure2_out <- paste0(
  "/Users/tf/e/ndmm/ESR2/submit/Communications Biology/",
  "Figure2_code/outputs"
)
dir.create(figure2_out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261013)

# Set this to a patient-level CSV containing the ORIGINAL microenvironment
# values, or keep df_6J in the current R session before sourcing this file.
# Required CSV columns:
# orig.ident, gender, group.a, ESR2_mean, mono_n, CD8_cyto
figure2_roc_csv <- ""

if (!exists("scRNA_harmony", inherits = TRUE)) {
  scRNA_harmony <- readRDS(
    "/Users/tf/e/ndmm/HM38/scRNA_harmony_HM38.anno.rds"
  )
}

# ===========================
# Shared input: align expression and metadata by cell barcode
# ===========================
figure2_meta <- scRNA_harmony[[]]
figure2_required <- c("orig.ident", "anno", "gender", "group", "group.a")
stopifnot(all(figure2_required %in% names(figure2_meta)))
if (!"RNA" %in% Assays(scRNA_harmony)) stop("RNA assay is missing.")
figure2_rna <- if (packageVersion("SeuratObject") >= "5.0.0") {
  GetAssayData(scRNA_harmony, assay = "RNA", layer = "data")
} else {
  GetAssayData(scRNA_harmony, assay = "RNA", slot = "data")
}
stopifnot("ESR2" %in% rownames(figure2_rna))
stopifnot(all(rownames(figure2_meta) %in% colnames(figure2_rna)))
figure2_cells <- figure2_meta[, figure2_required, drop = FALSE]
figure2_cells$ESR2 <- as.numeric(
  figure2_rna["ESR2", rownames(figure2_meta)]
)
figure2_cells <- figure2_cells %>%
  mutate(across(all_of(figure2_required), as.character))
if (any(!is.finite(figure2_cells$ESR2))) stop("Non-finite ESR2 expression.")
if (any(!figure2_cells$gender %in% c("Female", "Male"))) {
  stop("gender must contain Female/Male without missing values.")
}
if (any(!figure2_cells$group %in% c("HD", "MM"))) {
  stop("group must contain HD/MM without missing values.")
}
if (any(is.na(figure2_cells$anno)) || any(is.na(figure2_cells$orig.ident))) {
  stop("Missing anno or orig.ident.")
}
rm(figure2_rna)

# ===========================
# Figure 2D: expression Z-score tiles + detection-rate dots
# Z-score is calculated across the observed sex/status groups WITHIN each
# cell type. Mean expression includes zeros; detection means ESR2 > 0.
# ===========================
figure2_conditions <- c("Female | HD", "Male | HD", "Female | MM", "Male | MM")
figure2_celltypes <- sort(unique(figure2_cells$anno))

figure2_row_z <- function(x) {
  s <- sd(x)
  if (length(x) < 2L || !is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

figure2_D_observed <- figure2_cells %>%
  group_by(anno, gender, group) %>%
  summarise(
    n_cells = n(), mean_expr = mean(ESR2),
    detect_rate = mean(ESR2 > 0), .groups = "drop"
  ) %>%
  mutate(condition = paste(gender, group, sep = " | ")) %>%
  group_by(anno) %>%
  mutate(mean_expr_z = figure2_row_z(mean_expr)) %>%
  ungroup()

figure2_D_data <- expand_grid(
  anno = figure2_celltypes, condition = figure2_conditions
) %>%
  left_join(figure2_D_observed, by = c("anno", "condition")) %>%
  mutate(
    anno = factor(anno, levels = figure2_celltypes),
    condition = factor(condition, levels = figure2_conditions)
  )

p2D <- ggplot(figure2_D_data, aes(condition, anno)) +
  # Absent cell-type/group combinations remain blank, not zero expression.
  geom_tile(
    data = filter(figure2_D_data, !is.na(n_cells)),
    aes(fill = mean_expr_z), colour = "grey85", linewidth = 0.3
  ) +
  geom_point(
    data = filter(figure2_D_data, !is.na(n_cells)),
    aes(size = detect_rate), shape = 21,
    fill = "black", colour = "black", stroke = 0.2
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C7362A",
    midpoint = 0, limits = c(-1.5, 1.5),
    breaks = seq(-1.5, 1.5, 0.5), oob = scales::squish,
    na.value = "grey50", name = "ESR2 expression\n(Z-score)"
  ) +
  scale_size_continuous(
    range = c(1, 6), limits = c(0, max(0.05, figure2_D_observed$detect_rate)),
    breaks = seq(0, 0.05, 0.01), labels = scales::label_number(accuracy = 0.01),
    name = "Detection rate"
  ) +
  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = "Cell type") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10), legend.position = "right"
  )

# ===========================
# Figure 2F/G: malignant plasma cells, OR vs SOR (including ESR2 == 0)
# Actual response field is group.a; group is HD/MM, and response contains
# clinical categories such as CR/VGPR/PR/SD/PD.
# ===========================
figure2_mpc <- figure2_cells %>%
  filter(group == "MM", anno == "Malignant plasma cells")
if (nrow(figure2_mpc) == 0L) stop("No MM Malignant plasma cells in anno.")
if (any(!figure2_mpc$group.a %in% c("OR", "SOR"))) {
  stop("MM MPC group.a must be OR/SOR.")
}
figure2_mpc <- figure2_mpc %>%
  mutate(
    group.a = factor(group.a, levels = c("OR", "SOR")),
    gender = factor(gender, levels = c("Female", "Male"))
  )

figure2_ttest <- function(d) {
  if (!all(c("OR", "SOR") %in% as.character(d$group.a))) {
    stop("Each F/G comparison requires both OR and SOR.")
  }
  # Welch's two-sided t-test, as in the recovered geom_signif(test=t.test).
  p <- t.test(ESR2 ~ group.a, data = d, var.equal = FALSE)$p.value
  data.frame(
    n_OR_cells = sum(d$group.a == "OR"),
    n_SOR_cells = sum(d$group.a == "SOR"), p_value = p,
    label = if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  )
}
figure2_F_stats <- figure2_ttest(figure2_mpc)
figure2_G_stats <- bind_rows(lapply(c("Female", "Male"), function(sex) {
  cbind(gender = sex, figure2_ttest(filter(figure2_mpc, gender == sex)))
}))
figure2_ymax <- max(3, ceiling(max(figure2_mpc$ESR2) * 1.1))

figure2_violin <- function(d, stats, facet = FALSE) {
  stats <- mutate(stats, xmin = 1, xmax = 2,
                  y_bracket = figure2_ymax * 0.95,
                  y_tip = figure2_ymax * 0.93)
  p <- ggplot(d, aes(group.a, ESR2)) +
    geom_violin(fill = "grey95", colour = "black", linewidth = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.35) +
    geom_point(
      position = position_jitter(width = 0.18, height = 0, seed = 20261013),
      size = 0.25, alpha = 0.55, colour = "black"
    ) +
    geom_segment(
      data = stats, aes(x = xmin, xend = xmax, y = y_bracket,
                        yend = y_bracket), inherit.aes = FALSE
    ) +
    geom_segment(
      data = stats, aes(x = xmin, xend = xmin, y = y_tip,
                        yend = y_bracket), inherit.aes = FALSE
    ) +
    geom_segment(
      data = stats, aes(x = xmax, xend = xmax, y = y_tip,
                        yend = y_bracket), inherit.aes = FALSE
    ) +
    geom_text(
      data = stats, aes(x = 1.5, y = figure2_ymax * 0.97, label = label),
      inherit.aes = FALSE, size = 4
    ) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = c(-0.1, figure2_ymax)) +
    labs(x = NULL, y = "ESR2 expression") +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      axis.title.y = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(), strip.text = element_text(size = 12)
    )
  if (facet) p <- p + facet_wrap(~gender, nrow = 1, drop = FALSE)
  p
}
p2F <- figure2_violin(figure2_mpc, figure2_F_stats)
p2G <- figure2_violin(figure2_mpc, figure2_G_stats, facet = TRUE)

# Each orig.ident must describe exactly one sex and response group.
figure2_patient_checks <- figure2_mpc %>%
  group_by(orig.ident) %>%
  summarise(n_sex = n_distinct(gender), n_response = n_distinct(group.a), .groups = "drop")
stopifnot(all(figure2_patient_checks$n_sex == 1L),
          all(figure2_patient_checks$n_response == 1L))
figure2_patients <- figure2_mpc %>%
  group_by(orig.ident, gender, group.a) %>%
  summarise(ESR2_mean = mean(ESR2), n_MPC = n(), .groups = "drop")

# ===========================
# Figure 2H/I: patient-level ROC helper
# Historical code: OR = 1; SOR = 0; pROC direction = auto.
# Standardization is performed separately WITHIN each sex after complete-case
# filtering. Do not substitute monocyte proportions for mono_n or invent a
# cytotoxicity signature: the saved source does not define these inputs.
# ===========================
figure2_z <- function(x, name) {
  if (length(x) < 2L || any(!is.finite(x)) || sd(x) == 0) {
    stop("Cannot standardize ", name, ": need finite values with nonzero variance.")
  }
  as.numeric(scale(x))
}

figure2_roc_panel <- function(patient_data, sex, include_microenv = TRUE,
                              direction = "auto") {
  required <- c("orig.ident", "gender", "group.a", "ESR2_mean")
  if (include_microenv) required <- c(required, "mono_n", "CD8_cyto")
  if (!all(required %in% names(patient_data))) {
    stop("ROC input is missing: ", paste(setdiff(required, names(patient_data)), collapse = ", "))
  }
  if (anyDuplicated(patient_data$orig.ident)) stop("ROC input must have one row per patient.")
  d <- patient_data[as.character(patient_data$gender) == sex, required, drop = FALSE]
  if (any(is.na(d$group.a)) || any(!d$group.a %in% c("OR", "SOR"))) {
    stop("ROC group.a must contain OR/SOR only.")
  }
  score_cols <- c("ESR2_mean", if (include_microenv) c("mono_n", "CD8_cyto"))
  if (any(!vapply(d[score_cols], is.numeric, logical(1)))) {
    stop("ROC scores must be numeric.")
  }
  keep <- complete.cases(d) & apply(as.matrix(d[score_cols]), 1, function(x) all(is.finite(x)))
  if (any(!keep)) warning(sex, ": excluding ", sum(!keep), " patients with missing/non-finite scores.")
  d <- d[keep, , drop = FALSE]
  if (!all(c("OR", "SOR") %in% as.character(d$group.a))) stop("ROC requires OR and SOR.")
  predictors <- list(ESR2 = d$ESR2_mean)
  if (include_microenv) {
    env <- figure2_z(d$mono_n, "mono_n") - figure2_z(d$CD8_cyto, "CD8_cyto")
    predictors$Microenv <- env
    predictors$Combined <- figure2_z(d$ESR2_mean, "ESR2_mean") + env
  }
  rocs <- lapply(predictors, function(score) {
    pROC::roc(as.integer(d$group.a == "OR"), score,
              levels = c(0, 1), direction = direction, quiet = TRUE)
  })
  aucs <- bind_rows(lapply(names(rocs), function(model) {
    data.frame(gender = sex, model = model, n = nrow(d),
               n_OR = sum(d$group.a == "OR"), n_SOR = sum(d$group.a == "SOR"),
               AUC = as.numeric(pROC::auc(rocs[[model]])),
               direction = rocs[[model]]$direction)
  }))
  curves <- bind_rows(lapply(names(rocs), function(model) {
    r <- rocs[[model]]
    data.frame(model = model, FPR = 1 - r$specificities, sensitivity = r$sensitivities)
  }))
  legend_labels <- setNames(
    paste0(aucs$model, " (AUC=", sprintf("%.2f", aucs$AUC), ")"), aucs$model
  )
  p <- ggplot(curves, aes(FPR, sensitivity, colour = model)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_path(linewidth = 0.65) +
    scale_colour_manual(
      values = c(ESR2 = "#4575B4", Microenv = "#FDAE61", Combined = "#D73027"),
      breaks = names(rocs), labels = legend_labels
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0.02)) +
    coord_equal() +
    labs(title = sex, x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.6),
      legend.position = "bottom", legend.text = element_text(size = 9)
    ) + guides(colour = guide_legend(ncol = 1))
  list(plot = p, auc = aucs, rocs = rocs, patients = d)
}

# Use original patient scores if available. Otherwise generate only the ESR2
# curve and a clearly incomplete input template for the two missing features.
figure2_original_roc <- NULL
if (nzchar(figure2_roc_csv)) {
  figure2_original_roc <- read.csv(figure2_roc_csv, check.names = FALSE)
} else if (exists("df_6J", inherits = TRUE)) {
  figure2_original_roc <- get("df_6J", inherits = TRUE)
} else if (all(vapply(c("esr2_malignant", "mono_prop", "cd8_score"),
                      exists, logical(1), inherits = TRUE))) {
  figure2_original_roc <- esr2_malignant %>%
    left_join(mono_prop, by = "orig.ident") %>%
    left_join(cd8_score, by = c("orig.ident", "group.a", "gender"))
}

if (!is.null(figure2_original_roc)) {
  figure2_H_result <- figure2_roc_panel(figure2_original_roc, "Female")
  figure2_I_result <- figure2_roc_panel(figure2_original_roc, "Male")
  figure2_roc_suffix <- ""
} else {
  figure2_H_result <- figure2_roc_panel(figure2_patients, "Female", include_microenv = FALSE)
  figure2_I_result <- figure2_roc_panel(figure2_patients, "Male", include_microenv = FALSE)
  figure2_roc_suffix <- "_ESR2_only"
  figure2_roc_template <- figure2_patients %>%
    mutate(mono_n = NA_real_, CD8_cyto = NA_real_)
  write.csv(figure2_roc_template,
            file.path(figure2_out, "Figure2HI_original_features_TEMPLATE.csv"), row.names = FALSE)
  message("H/I: original mono_n/CD8_cyto are absent. Exporting ESR2-only curves and an incomplete template.")
}
p2H <- figure2_H_result$plot
p2I <- figure2_I_result$plot

# ===========================
# Exports and numerical source data
# ===========================
figure2_save <- function(p, name, width, height) {
  ggsave(file.path(figure2_out, paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(figure2_out, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
}
figure2_save(p2D, "Figure2D", 7.3, 8.2)
figure2_save(p2F, "Figure2F", 2.8, 4.6)
figure2_save(p2G, "Figure2G", 5.2, 4.6)
figure2_save(p2H, paste0("Figure2H", figure2_roc_suffix), 4.2, 4.8)
figure2_save(p2I, paste0("Figure2I", figure2_roc_suffix), 4.2, 4.8)
write.csv(figure2_D_data, file.path(figure2_out, "Figure2D_source_data.csv"), row.names = FALSE)
write.csv(bind_rows(mutate(figure2_F_stats, comparison = "All"),
                    mutate(figure2_G_stats, comparison = gender)),
          file.path(figure2_out, "Figure2FG_cell_level_Welch_tests.csv"), row.names = FALSE)
write.csv(figure2_patients, file.path(figure2_out, "Figure2_patient_ESR2.csv"), row.names = FALSE)
write.csv(bind_rows(figure2_H_result$auc, figure2_I_result$auc),
          file.path(figure2_out, paste0("Figure2HI_AUC", figure2_roc_suffix, ".csv")), row.names = FALSE)
saveRDS(list(Female = figure2_H_result, Male = figure2_I_result),
        file.path(figure2_out, paste0("Figure2HI_ROC", figure2_roc_suffix, ".rds")))
capture.output(sessionInfo(), file = file.path(figure2_out, "sessionInfo.txt"))
print(bind_rows(figure2_H_result$auc, figure2_I_result$auc))
message("Saved panels and source data in: ", figure2_out)

# In RStudio, display plots individually with:
# p2D; p2F; p2G; p2H; p2I
# p2F | p2G | p2H | p2I

