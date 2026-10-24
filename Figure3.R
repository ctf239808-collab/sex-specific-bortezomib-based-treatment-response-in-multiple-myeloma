
#Figure 3

# ======================== 1. User-configurable parameters ========================
figure3_analysis_file <- "/Users/tf/e/ndmm/HM38/scRNA_harmony_MM18_MM.dim50res0.5.anno.rds"
# A-C may use a separate historical object with the original coordinates; NULL uses the D-F object.
figure3_abc_file <- "/Users/tf/e/ndmm/HM38/scRNA_harmony_MM18.MM.anno.rds"
figure3_output_dir <- "/Users/tf/e/ndmm/ESR2/submit/NC/Figure3_code/outputs"
figure3_formats <- c("pdf", "png", "svg")  # Add "tiff" if needed; all formats use built-in Cairo devices.
figure3_width_mm <- 183
figure3_height_mm <- 174
figure3_dpi <- 300

# ======================== 2. Main function and dependency checks ========================
build_figure3 <- function(
    analysis_file, abc_file = NULL, output_dir,
    formats = c("pdf", "png", "svg"), width_mm = 183, height_mm = 174,
    dpi = 300, minimum_cells_heatmap = 20L) {
  required <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "scales")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages in this R installation: ", paste(missing, collapse = ", "),
         ". Use an R installation with these packages, or run install.packages().")
  }
  if (!all(formats %in% c("pdf", "png", "svg", "tiff"))) {
    stop("Supported formats are pdf, png, svg, and tiff.")
  }
  if (!capabilities("cairo")) {
    stop("Graphics export requires an R installation with Cairo support; on this machine, use /Library/Frameworks/R.framework/Resources/bin/Rscript.")
  }
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
  })
  if (!file.exists(analysis_file)) stop("Analysis input file not found: ", analysis_file)
  if (is.null(abc_file)) abc_file <- analysis_file
  if (!file.exists(abc_file)) stop("A-C input file not found: ", abc_file)
  source_dir <- file.path(output_dir, "source_data")
  panel_dir <- file.path(output_dir, "panels")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read the specified file each time rather than reusing a potentially modified session object.
  message("Reading the annotated malignant plasma-cell object...")
  obj <- readRDS(analysis_file)
  if (!inherits(obj, "Seurat")) stop("analysis_file must contain a Seurat object.")
  required_meta <- c("orig.ident", "gender", "group.a", "anno3")
  meta <- obj[[]]
  if (!all(required_meta %in% names(meta))) stop("Required cell metadata columns are missing.")
  meta <- meta[, required_meta, drop = FALSE]
  meta[] <- lapply(meta, as.character)
  if (anyNA(meta) || any(as.matrix(meta) == "")) stop("Cell metadata contain missing values.")
  if (!all(meta$gender %in% c("Female", "Male")) ||
      !all(meta$group.a %in% c("OR", "SOR"))) stop("Sex or response group labels do not match the original figure.")
  anno_order <- c("MM", "MM_ANXA1", "MM_AZGP1", "MM_CRIP1", "MM_IGHGP",
                  "MM_IGLL5", "MM_IL5RA", "PBM", "PBM_CCND1")
  if (!setequal(unique(meta$anno3), anno_order)) stop("anno3 labels do not match the nine subclusters in the original figure.")
  get_rna <- function(x, layer) {
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      SeuratObject::GetAssayData(x, assay = "RNA", layer = layer)
    } else {
      SeuratObject::GetAssayData(x, assay = "RNA", slot = layer)
    }
  }
  counts <- get_rna(obj, "counts")
  if (!"ESR2" %in% rownames(counts)) stop("ESR2 was not found in RNA counts.")
  cell_index <- match(rownames(meta), colnames(counts))
  if (anyNA(cell_index)) stop("Cell names in RNA counts do not match the metadata.")
  meta$ESR2_count <- as.numeric(counts["ESR2", cell_index])
  meta$ESR2_positive <- meta$ESR2_count > 0
  if (any(!is.finite(meta$ESR2_count)) || any(meta$ESR2_count < 0)) {
    stop("ESR2 counts must be finite and non-negative.")
  }
  subcluster_order <- names(sort(table(meta$anno3), decreasing = TRUE))
  
  # ======================== 3. Shared styling and patient metadata ========================
  sex_colors <- c(Female = "#B74E7A", Male = "#2F6F9F")
  response_colors <- c(OR = "#2A9D8F", SOR = "#D77A2B")
  anno_colors <- stats::setNames(scales::hue_pal()(length(anno_order)), anno_order)
  theme_fig3 <- function() {
    theme_classic(base_size = 6.5, base_family = "sans") +
      theme(
        axis.text = element_text(size = 6, colour = "black"),
        axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks = element_line(linewidth = 0.25),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5.7),
        legend.key.size = grid::unit(2.6, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(1.3, "mm"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0, 1),
        plot.margin = margin(4, 4, 4, 4)
      )
  }
  patient_info <- unique(meta[, c("orig.ident", "gender", "group.a")])
  if (anyDuplicated(patient_info$orig.ident)) {
    stop("At least one patient has multiple sex or response labels.")
  }
  patient_info <- patient_info[order(
    match(patient_info$gender, c("Female", "Male")),
    match(patient_info$group.a, c("OR", "SOR")), patient_info$orig.ident), ]
  rownames(patient_info) <- NULL
  patient_info$plot_index <- seq_len(nrow(patient_info))
  patient_order <- patient_info$orig.ident
  
  # ======================== 4. Figure 3A-C ========================
  # Use the input object's existing RNA/data layer, UMAP coordinates, and anno3 annotations.
  abc <- if (identical(normalizePath(abc_file), normalizePath(analysis_file))) {
    obj
  } else {
    readRDS(abc_file)
  }
  if (!inherits(abc, "Seurat") || !"anno3" %in% names(abc[[]])) {
    stop("The A-C input must contain a Seurat object with anno3 annotations.")
  }
  if (!setequal(unique(as.character(abc$anno3)), anno_order)) {
    stop("anno3 labels in the A-C input do not match the nine subclusters in the original figure.")
  }
  if (!setequal(colnames(abc), rownames(meta))) {
    stop("The A-C and D-F inputs contain different cell sets and cannot be directly combined.")
  }
  abc_meta <- abc[[]][match(rownames(meta), rownames(abc[[]])), required_meta, drop = FALSE]
  if (!all(vapply(required_meta, function(nm) {
    identical(as.character(abc_meta[[nm]]), as.character(meta[[nm]]))
  }, logical(1)))) stop("Patient, sex, response, or subcluster annotations differ between the A-C and D-F inputs.")
  abc_counts <- get_rna(abc, "counts")
  if (!identical(as.numeric(abc_counts["ESR2", rownames(meta)]), meta$ESR2_count)) {
    stop("Raw ESR2 counts differ between the A-C and D-F inputs.")
  }
  if (!"umap" %in% SeuratObject::Reductions(abc)) stop("The A-C input has no saved UMAP coordinates.")
  SeuratObject::DefaultAssay(abc) <- "RNA"
  abc$anno3 <- factor(as.character(abc$anno3), levels = anno_order)
  abc_data <- get_rna(abc, "data")
  marker_genes <- c("MZB1", "SDC1", "JCHAIN", "IGHG1", "IGHG3", "COBLL1",
                    "IGHA1", "IGHA2", "PRDM1", "IRF4", "XBP1", "MKI67", "RRM2",
                    "ANXA1", "CFI", "IL5RA", "SCYL2", "CRIP1", "AZGP1",
                    "TMEM176B", "IGLL5", "IGLV2-23", "CCND1", "WFDC2", "IGHGP", "IFI27")
  missing_genes <- setdiff(c("ESR2", marker_genes), rownames(abc_data))
  if (length(missing_genes)) stop("Genes missing from the A-C input: ", paste(missing_genes, collapse = ", "))
  umap <- SeuratObject::Embeddings(abc, reduction = "umap")[, 1:2, drop = FALSE]
  umap_source <- data.frame(
    cell = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
    anno3 = as.character(abc$anno3[match(rownames(umap), colnames(abc))]),
    ESR2_normalized = as.numeric(abc_data["ESR2", rownames(umap)]),
    row.names = NULL
  )
  # A: UMAP colored by the nine subclusters, with subcluster labels.
  p3A <- Seurat::DimPlot(
    abc, reduction = "umap", group.by = "anno3", cols = anno_colors,
    label = TRUE, label.size = 2.0, repel = FALSE, pt.size = 0.10,
    raster = FALSE
  ) + labs(title = NULL, tag = "A", x = "UMAP_1", y = "UMAP_2", colour = NULL) +
    theme_fig3() + guides(colour = guide_legend(override.aes = list(size = 1.5)))
  
  # B: Seurat DotPlot. Color: Z-score of average expression; size: percentage of expressing cells.
  # Seurat applies expm1 to log-normalized data, averages within subclusters, then applies log1p and scales each gene.
  marker_top_to_bottom <- c("MM_IGHGP", "PBM_CCND1", "MM_IGLL5", "MM_AZGP1",
                            "MM_CRIP1", "MM_IL5RA", "MM_ANXA1", "PBM", "MM")
  p3B <- Seurat::DotPlot(
    abc, features = marker_genes, assay = "RNA", group.by = "anno3",
    scale = TRUE, col.min = -2.5, col.max = 2.5, dot.scale = 3.3
  ) +
    scale_y_discrete(limits = rev(marker_top_to_bottom)) +
    scale_colour_gradient2(
      low = "#4682B4", mid = "white", high = "#B22222", midpoint = 0,
      limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2),
      name = "Average\nexpression"
    ) +
    labs(x = "Features", y = "Identity", tag = "B") + theme_fig3() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5.8)) +
    guides(size = guide_legend(title = "Percent\nexpressed", order = 2),
           colour = guide_colourbar(order = 1, barheight = grid::unit(15, "mm")))
  
  # C: Use the full RNA/data expression range without quantile cutoffs.
  p3C <- Seurat::FeaturePlot(
    abc, features = "ESR2", reduction = "umap", cols = c("lightgrey", "firebrick"),
    order = FALSE, pt.size = 0.10, raster = FALSE
  ) + labs(title = "ESR2", tag = "C", x = "UMAP_1", y = "UMAP_2", colour = NULL) +
    theme_fig3() + guides(colour = guide_colourbar(barheight = grid::unit(18, "mm")))
  
  # ======================== 5. Figure 3D ========================
  heatmap_data <- expand.grid(
    orig.ident = patient_order, anno3 = subcluster_order,
    stringsAsFactors = FALSE
  )
  heatmap_data$n_cells <- 0L
  heatmap_data$ESR2_positive_cells <- 0L
  for (i in seq_len(nrow(heatmap_data))) {
    idx <- meta$orig.ident == heatmap_data$orig.ident[i] & meta$anno3 == heatmap_data$anno3[i]
    heatmap_data$n_cells[i] <- sum(idx)
    heatmap_data$ESR2_positive_cells[i] <- sum(meta$ESR2_positive[idx])
  }
  heatmap_data$ESR2_positive_percent <- with(
    heatmap_data, ifelse(n_cells > 0, 100 * ESR2_positive_cells / n_cells, NA_real_))
  heatmap_data$display_percent <- with(
    heatmap_data, ifelse(n_cells >= minimum_cells_heatmap, ESR2_positive_percent, NA_real_))
  heatmap_data$patient_index <- match(heatmap_data$orig.ident, patient_order)
  heatmap_data$y <- length(subcluster_order) + 1L - match(heatmap_data$anno3, subcluster_order)
  heatmap_data$gender <- patient_info$gender[match(heatmap_data$orig.ident, patient_order)]
  heatmap_data$group.a <- patient_info$group.a[match(heatmap_data$orig.ident, patient_order)]
  heat_max <- max(heatmap_data$display_percent, na.rm = TRUE)
  if (!is.finite(heat_max)) stop("No patient-subcluster combination meets the cell-count threshold for panel D.")
  annotation_colors <- c(
    "Response: OR" = response_colors[["OR"]], "Response: SOR" = response_colors[["SOR"]],
    "Sex: Female" = sex_colors[["Female"]], "Sex: Male" = sex_colors[["Male"]],
    "<20 cells" = "#D9D9D9"
  )
  names(annotation_colors)[5] <- paste0("<", minimum_cells_heatmap, " cells")
  annotation_keys <- data.frame(
    x = NA_real_, y = NA_real_, key = factor(names(annotation_colors), levels = names(annotation_colors)))
  n_sub <- length(subcluster_order)
  p3D <- ggplot(heatmap_data, aes(patient_index, y)) +
    geom_tile(aes(fill = display_percent), width = 0.98, height = 0.98,
              colour = "white", linewidth = 0.15) +
    annotate("rect", xmin = patient_info$plot_index - 0.49, xmax = patient_info$plot_index + 0.49,
             ymin = n_sub + 0.47, ymax = n_sub + 0.76,
             fill = unname(sex_colors[patient_info$gender])) +
    annotate("rect", xmin = patient_info$plot_index - 0.49, xmax = patient_info$plot_index + 0.49,
             ymin = n_sub + 0.83, ymax = n_sub + 1.12,
             fill = unname(response_colors[patient_info$group.a])) +
    annotate("text", x = 0.35, y = n_sub + c(0.61, 0.98),
             label = c("Sex", "Response"), hjust = 1, size = 2.0) +
    geom_point(data = annotation_keys, aes(x, y, colour = key), shape = 15,
               inherit.aes = FALSE, na.rm = TRUE, show.legend = TRUE) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594"),
      limits = c(0, heat_max), na.value = "#D9D9D9",
      breaks = c(0, heat_max), labels = sprintf("%.1f", c(0, heat_max)), name = "ESR2+ (%)"
    ) +
    scale_colour_manual(values = annotation_colors, drop = FALSE, name = NULL) +
    scale_x_continuous(breaks = seq_along(patient_order), labels = patient_order,
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_len(n_sub), labels = rev(subcluster_order),
                       limits = c(0.5, n_sub + 1.2), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, length(patient_order) + 0.5), clip = "off") +
    labs(x = "Patient ID", y = NULL, tag = "D") +
    theme_fig3() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.7),
          axis.text.y = element_text(size = 5.8)) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2)),
           fill = guide_colourbar(order = 2, barheight = grid::unit(14, "mm")))
  
  # ======================== 6. E/F proportions, confidence intervals, and tests ========================
  wilson_ci <- function(successes, total) {
    if (total <= 0L) return(c(lower = NA_real_, upper = NA_real_))
    z <- stats::qnorm(0.975)
    p <- successes / total
    den <- 1 + z^2 / total
    center <- (p + z^2 / (2 * total)) / den
    half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / den
    c(lower = max(0, center - half), upper = min(1, center + half))
  }
  specs <- list(
    gender = list(variable = "gender", groups = c("Female", "Male")),
    group_a = list(variable = "group.a", groups = c("OR", "SOR"))
  )
  proportion_rows <- list()
  test_rows <- list()
  for (comparison in names(specs)) {
    spec <- specs[[comparison]]
    for (subcluster in subcluster_order) {
      group_counts <- matrix(0L, nrow = 2L, ncol = 2L,
                             dimnames = list(spec$groups, c("positive", "negative")))
      for (j in seq_along(spec$groups)) {
        group_name <- spec$groups[j]
        idx <- meta$anno3 == subcluster & meta[[spec$variable]] == group_name
        n_cells <- sum(idx)
        n_pos <- sum(meta$ESR2_positive[idx])
        ci <- 100 * wilson_ci(n_pos, n_cells)
        group_counts[j, ] <- c(n_pos, n_cells - n_pos)
        proportion_rows[[length(proportion_rows) + 1L]] <- data.frame(
          comparison = comparison, group = group_name, anno3 = subcluster,
          patients_with_cells = length(unique(meta$orig.ident[idx])),
          cells = n_cells, ESR2_positive_cells = n_pos,
          ESR2_negative_cells = n_cells - n_pos,
          ESR2_positive_percent = if (n_cells > 0L) 100 * n_pos / n_cells else NA_real_,
          ESR2_positive_percent_CI_low = ci[["lower"]],
          ESR2_positive_percent_CI_high = ci[["upper"]], row.names = NULL
        )
      }
      # For Fisher's test, place the contrast (Male/SOR) first and the reference (Female/OR) second.
      # This preserves the historical odds-ratio direction. Pass zero row margins to fisher.test without adding cells.
      ft <- stats::fisher.test(group_counts[2:1, , drop = FALSE],
                               alternative = "two.sided", conf.level = 0.95)
      test_rows[[length(test_rows) + 1L]] <- data.frame(
        comparison = comparison, anno3 = subcluster,
        reference_group = spec$groups[1], contrast_group = spec$groups[2],
        reference_cells = sum(group_counts[1, ]), contrast_cells = sum(group_counts[2, ]),
        reference_ESR2_positive_cells = group_counts[1, 1],
        contrast_ESR2_positive_cells = group_counts[2, 1],
        odds_ratio_fisher = unname(ft$estimate), odds_ratio_CI_low = ft$conf.int[1],
        odds_ratio_CI_high = ft$conf.int[2], fisher_exact_p = ft$p.value,
        row.names = NULL
      )
    }
  }
  group_summary <- do.call(rbind, proportion_rows)
  fisher_results <- do.call(rbind, test_rows)
  fisher_results$fisher_exact_q_BH_within_comparison <- NA_real_
  for (comparison in names(specs)) {
    idx <- fisher_results$comparison == comparison
    fisher_results$fisher_exact_q_BH_within_comparison[idx] <- stats::p.adjust(
      fisher_results$fisher_exact_p[idx], method = "BH")
  }
  fisher_results$inference_unit <- "Pooled cells; within-patient dependence is not modeled"
  stable <- group_summary$cells >= 20L
  proportion_xmax <- max(5, ceiling(1.12 * max(
    group_summary$ESR2_positive_percent_CI_high[stable], na.rm = TRUE)))
  format_p <- function(p) {
    if (!is.finite(p)) return("NA")
    if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
    sprintf("%.3f", p)
  }
  format_q <- function(q) {
    if (!is.finite(q)) return("q=NA")
    if (q < 0.001) return("q<0.001")
    paste0("q=", sprintf("%.3f", q))
  }
  make_proportion_panel <- function(comparison, tag) {
    spec <- specs[[comparison]]
    dat <- group_summary[group_summary$comparison == comparison, ]
    stat <- fisher_results[fisher_results$comparison == comparison, ]
    dat$group <- factor(dat$group, levels = spec$groups)
    dat$y <- n_sub + 1 - match(dat$anno3, subcluster_order) +
      ifelse(dat$group == spec$groups[1], 0.12, -0.12)
    dat$shape_id <- ifelse(dat$cells < 20L,
                           ifelse(dat$group == spec$groups[1], "open_circle", "open_triangle"),
                           ifelse(dat$group == spec$groups[1], "circle", "triangle"))
    dat$clipped_high <- pmin(proportion_xmax, dat$ESR2_positive_percent_CI_high)
    stat$y <- n_sub + 1 - match(stat$anno3, subcluster_order)
    stat$stat_label <- paste0("p=", vapply(stat$fisher_exact_p, format_p, character(1)),
                              ", ", vapply(stat$fisher_exact_q_BH_within_comparison, format_q, character(1)))
    stat$text_colour <- ifelse(stat$fisher_exact_q_BH_within_comparison < 0.05,
                               "#D73027", "#3F3F3F")
    group_colors <- if (comparison == "gender") sex_colors else response_colors
    clipped <- dat[is.finite(dat$ESR2_positive_percent_CI_high) &
                     dat$ESR2_positive_percent_CI_high > proportion_xmax, ]
    p <- ggplot(dat, aes(ESR2_positive_percent, y, colour = group)) +
      geom_segment(aes(x = ESR2_positive_percent_CI_low, xend = clipped_high, yend = y),
                   linewidth = 0.4, na.rm = TRUE) +
      geom_segment(aes(x = ESR2_positive_percent_CI_low, xend = ESR2_positive_percent_CI_low,
                       y = y - 0.09, yend = y + 0.09), linewidth = 0.3, na.rm = TRUE) +
      geom_segment(aes(x = clipped_high, xend = clipped_high, y = y - 0.09, yend = y + 0.09),
                   linewidth = 0.3, na.rm = TRUE) +
      geom_point(aes(shape = shape_id), size = 1.5, na.rm = TRUE) +
      geom_text(data = clipped, aes(x = proportion_xmax, y = y, colour = group),
                label = ">", hjust = 1, size = 1.8, inherit.aes = FALSE, show.legend = FALSE) +
      geom_text(data = stat, aes(x = proportion_xmax + 0.7, y = y, label = stat_label),
                colour = stat$text_colour, hjust = 0, size = 1.95,
                inherit.aes = FALSE) +
      annotate("text", x = proportion_xmax + 0.7, y = n_sub + 0.8,
               label = "p value, q value", hjust = 0, size = 2.0) +
      scale_colour_manual(values = group_colors,
                          name = if (comparison == "gender") "Sex" else "Response") +
      scale_shape_manual(values = c(circle = 16, triangle = 17, open_circle = 1, open_triangle = 2),
                         guide = "none") +
      scale_x_continuous(breaks = seq(0, proportion_xmax, by = 5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq_len(n_sub), labels = rev(subcluster_order),
                         expand = c(0, 0)) +
      coord_cartesian(xlim = c(0, proportion_xmax), ylim = c(0.5, n_sub + 0.5), clip = "off") +
      labs(x = "ESR2-positive cells (%)", y = NULL, tag = tag) +
      theme_fig3() +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_text(size = 6.2),
            legend.position = "top", legend.direction = "horizontal",
            plot.margin = margin(7, 61, 5, 4)) +
      guides(colour = guide_legend(override.aes = list(shape = c(16, 17), size = 1.8)))
    p
  }
  p3E <- make_proportion_panel("gender", "E")
  p3F <- make_proportion_panel("group_a", "F")
  
  # ======================== 7. A-F assembly and export ========================
  # Allocate widths to complete panels so the right margins for E/F statistics do not compress B's gene labels.
  # The original ggplot objects remain available in panels for independent editing.
  panel_block <- function(p) patchwork::wrap_elements(full = ggplot2::ggplotGrob(p))
  top_row <- (panel_block(p3A) | panel_block(p3B)) + plot_layout(widths = c(0.43, 0.57))
  middle_row <- (panel_block(p3C) | panel_block(p3D)) + plot_layout(widths = c(0.36, 0.64))
  bottom_row <- (panel_block(p3E) | panel_block(p3F)) + plot_layout(widths = c(1, 1))
  figure <- (top_row / middle_row / bottom_row) +
    plot_layout(heights = c(1.05, 1.05, 1.12))
  panels <- list(A = p3A, B = p3B, C = p3C, D = p3D, E = p3E, F = p3F)
  write_source <- function(x, name) {
    utils::write.csv(x, file.path(source_dir, name), row.names = FALSE, na = "")
  }
  write_source(umap_source, "Figure3AC_UMAP_ESR2_source.csv")
  write_source(p3B$data, "Figure3B_marker_DotPlot_source.csv")
  write_source(patient_info, "Figure3_patient_metadata.csv")
  write_source(heatmap_data, "Figure3D_patient_subcluster_ESR2_source.csv")
  write_source(group_summary, "Figure3EF_pooled_proportions_Wilson95CI.csv")
  write_source(fisher_results, "Figure3EF_Fisher_BH_tests.csv")
  input_manifest <- data.frame(
    panels = c("A-C", "D-F"), input_file = normalizePath(c(abc_file, analysis_file)),
    input_md5 = unname(tools::md5sum(c(abc_file, analysis_file))),
    cells = c(ncol(abc), nrow(meta)),
    same_cell_set = setequal(colnames(abc), rownames(meta)),
    row.names = NULL
  )
  write_source(input_manifest, "Figure3_input_manifest.csv")
  save_plot <- function(plot, stem, w_mm, h_mm) {
    for (format in formats) {
      filename <- paste0(stem, ".", format)
      # Open devices directly, bypassing ggsave device selection and ragg/svglite/systemfonts.
      # This avoids locate_fonts_c / SET_VECTOR_ELT font errors in some interactive R sessions.
      # Cairo PDF retains selectable text; Cairo SVG commonly stores text as glyph paths.
      export_one <- function() {
        switch(format,
               pdf = grDevices::cairo_pdf(
                 filename = filename, width = w_mm / 25.4, height = h_mm / 25.4,
                 family = "sans", bg = "white", fallback_resolution = dpi),
               png = grDevices::png(
                 filename = filename, width = w_mm / 25.4, height = h_mm / 25.4,
                 units = "in", res = dpi, type = "cairo", family = "sans", bg = "white"),
               svg = grDevices::svg(
                 filename = filename, width = w_mm / 25.4, height = h_mm / 25.4,
                 family = "sans", bg = "white", onefile = FALSE),
               tiff = grDevices::tiff(
                 filename = filename, width = w_mm / 25.4, height = h_mm / 25.4,
                 units = "in", res = dpi, type = "cairo", compression = "lzw",
                 family = "sans", bg = "white")
        )
        device_id <- grDevices::dev.cur()
        on.exit({
          if (device_id %in% grDevices::dev.list()) grDevices::dev.off(device_id)
        }, add = TRUE)
        print(plot)
      }
      tryCatch(export_one(), error = function(e) {
        stop("Export failed [", filename, "]: ", conditionMessage(e), call. = FALSE)
      })
    }
  }
  message("Exporting the integrated Figure 3 and six individual panels...")
  save_plot(figure, file.path(output_dir, "Figure3_integrated"), width_mm, height_mm)
  panel_dims <- list(A = c(84, 65), B = c(108, 65), C = c(74, 65),
                     D = c(123, 65), E = c(94, 67), F = c(94, 67))
  for (tag in names(panels)) {
    dims <- panel_dims[[tag]]
    save_plot(panels[[tag]], file.path(panel_dir, paste0("Figure3", tag)), dims[1], dims[2])
  }
  legend <- c(
    "Figure 3 | Heterogeneous ESR2 expression across malignant plasma-cell subclusters.",
    "A: Saved UMAP coordinates and nine anno3 annotations; no reclustering or new embedding is performed.",
    "B: Canonical-marker DotPlot from the RNA/data layer. Dot size: percent expressing cells; colour: gene-wise scaled average expression (Seurat DotPlot, clipped to [-2.5, 2.5]).",
    "C: ESR2 log-normalized RNA expression on the same saved UMAP, without quantile cutoffs.",
    sprintf("D: Raw RNA count >0 defines ESR2-positive cells. Each tile is one patient-subcluster percentage; < %d cells or absent subclusters are grey. Patients are ordered by sex, response, then ID.", minimum_cells_heatmap),
    "E/F: ESR2-positive proportions in pooled cells by sex/response; error bars: Wilson 95% confidence intervals. Open symbols indicate fewer than 20 pooled cells; > indicates a confidence interval extending beyond the displayed axis.",
    "P: two-sided Fisher exact test of ESR2-positive/negative counts. q: Benjamini-Hochberg correction separately across nine subclusters within each comparison family.",
    "E/F treat cells as observations and ignore within-patient dependence. These are exploratory pooled-cell comparisons and do not establish independent patient-level associations.",
    sprintf("D-F input: %d cells, %d patients, %d ESR2-positive cells.", nrow(meta), nrow(patient_info), sum(meta$ESR2_positive)),
    "This integrated layout is reconstructed from the available processed objects and plotting calculations; it is not the recovered historical composite source."
  )
  writeLines(legend, file.path(output_dir, "Figure3_legend_and_methods.txt"))
  capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  writeLines(c("All formats use built-in R grDevices with Cairo.",
               "No ragg/svglite/systemfonts font lookup is used for exporting.",
               "PDF retains selectable text; Cairo SVG may store text as glyph paths.",
               paste("R:", R.version.string)),
             file.path(output_dir, "export_device_info.txt"))
  message(sprintf("Completed: %d cells, %d patients, %d ESR2-positive cells; outputs: %s",
                  nrow(meta), nrow(patient_info), sum(meta$ESR2_positive), output_dir))
  invisible(list(figure = figure, panels = panels, patient_metadata = patient_info,
                 heatmap_data = heatmap_data, proportions = group_summary,
                 fisher_results = fisher_results, input_manifest = input_manifest))
}

# ======================== 8. Run and retain editable ggplot objects ========================
figure3_result <- build_figure3(
  analysis_file = figure3_analysis_file, abc_file = figure3_abc_file,
  output_dir = figure3_output_dir, formats = figure3_formats,
  width_mm = figure3_width_mm, height_mm = figure3_height_mm, dpi = figure3_dpi
)
p3A <- figure3_result$panels$A
p3B <- figure3_result$panels$B
p3C <- figure3_result$panels$C
p3D <- figure3_result$panels$D
p3E <- figure3_result$panels$E
p3F <- figure3_result$panels$F
Figure3 <- figure3_result$figure
if (interactive()) print(Figure3)
