#!/usr/bin/env Rscript
# Figure5.new.R | Six CellChat panels after removing the former pathway panel E.
# Run in RStudio: source("/Users/tf/e/ndmm/ESR2/submit/Communications Biology/Figure5.new.R")
# Or use /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript.
# To load functions only: Sys.setenv(FIGURE5_RUN = "0"); source(...)
#
# Panel mapping from Figure5_rebuilt.pdf: A->A, B->B, C->C, D->D, F->E, G->F.
# The former E, its reference image and its exploratory pathway method are removed.
# This script is self-contained and does not source Figure.5.R or read old PDFs.
# Existing network caches are the default input; no CellChat inference is run.
# To rerun the full input pipeline, set figure5_use_cached_networks = FALSE.
# Set figure5_rebuild_networks = TRUE only to force fresh CellChat inference.
#
# Scientific definitions are retained from Figure.5.R:
# A/B use network counts; D uses weights and omits filterCommunication as before.
# C retains the historical receiving-axis/bar definition (legacy mode).
# E (formerly F) retains the documented reconstructed MPC/immune focus filter;
# the original exact focus-filter source was not recovered.
# F (formerly G) retains the historical within-OR/SOR comparisons and filters.
# E/F use their original distinct absent-sex and zero-p-value conventions.
# pval_min is WITHIN-NETWORK interaction support, not a between-sex test.
# A-C/E pool all atlas donors by sex; D also splits by global-median ESR2;
# F uses only OR/SOR cells. All sex comparisons remain descriptive.

# ======================== 1. Parameters ========================
figure5_project_dir <- "/Users/tf/e/ndmm/ESR2/submit/Communications Biology"
figure5_input_file <- "/Users/tf/e/ndmm/HM38/scRNA_harmony_HM38.anno.rds"
figure5_output_dir <- file.path(figure5_project_dir, "Figure5_code", "outputs_new")
figure5_figure_stem <- file.path(figure5_project_dir, "Figure5.new")
figure5_cache_dir <- file.path(figure5_project_dir, "Figure5_code", "cache")
figure5_formats <- c("pdf", "png")
figure5_export_panels <- FALSE
figure5_use_cached_networks <- TRUE
figure5_verify_cached_input_md5 <- TRUE
figure5_expected_input_md5 <- "83784df6f57f7d0a6ac493f19da750bc"
figure5_rebuild_networks <- FALSE
figure5_width_mm <- 270
figure5_height_mm <- 380
figure5_dpi <- 300

# Physical row sizes sum to the page height. D's four networks now occupy
# the full middle-left block, making each close in size to A/B.
figure5_layout_top_mm <- 95
figure5_layout_middle_mm <- 165
figure5_layout_bottom_mm <- 120
figure5_left_column_mm <- 110.4
figure5_ab_height_mm <- 79.8
figure5_ab_edge_width_max <- 0.9
figure5_d_edge_width_max <- 0.35
figure5_edge_alpha <- 0.35
figure5_arrow_width <- 0.30
figure5_arrow_size <- 0.06
figure5_network_label_cex <- 0.45
figure5_network_label_wrap <- 16L
figure5_network_label_y_extent <- 1.40
figure5_network_title_cex <- 0.70
figure5_network_leader_width <- 0.45
figure5_c_x_label_angle <- 65
figure5_c_standalone_width_mm <- 180
figure5_c_standalone_height_mm <- 105
figure5_f_y_title_margin_pt <- 16
figure5_f_standalone_width_mm <- 240
figure5_seed <- 1L
figure5_nboot <- 100L
figure5_min_cells <- 10L
figure5_p_threshold <- 0.05
figure5_c_mode <- "legacy"  # Or "corrected": separate source/target labels.
figure5_e_top_cellpairs <- 12L
figure5_e_min_abs_log2fc <- 0.25
figure5_e_scope <- "mpc_immune"  # Or "historical_all_pairs".
figure5_f_top_per_cellpair <- 5L
figure5_mpc <- "Malignant plasma cells"
figure5_immune_cells <- c("CD8+ T cells", "CD4+ T cells", "NK cells",
  "CD14+ monocytes", "CD16+ Monocytes", "Macrophages", "Intermediate monocytes")

# ======================== 2. Dependencies and input ========================
figure5_check_packages <- function() {
  required <- c("SeuratObject", "CellChat", "Matrix", "dplyr", "tidyr",
    "ggplot2", "patchwork", "gridGraphics", "future", "digest", "igraph")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
    ". On this machine use /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript.")
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
  })
  if (!capabilities("cairo")) stop("Cairo graphics support is required for export.")
}

figure5_get_expression <- function(obj) {
  # Use the existing log-normalized RNA data, not SCT/integrated values or counts.
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    assay <- obj[["RNA"]]
    if (inherits(assay, "Assay5")) {
      layers <- SeuratObject::Layers(assay, search = "^data($|\\.)")
      if (!length(layers)) stop("RNA assay has no normalized data layer.")
      if (length(layers) > 1L) {
        obj <- SeuratObject::JoinLayers(obj, assay = "RNA", layers = layers, new = "data")
        layer <- "data"
      } else layer <- layers[[1L]]
    } else layer <- "data"
    return(SeuratObject::GetAssayData(obj, assay = "RNA", layer = layer))
  }
  SeuratObject::GetAssayData(obj, assay = "RNA", slot = "data")
}

figure5_read_input <- function(path) {
  if (!file.exists(path)) stop("Input RDS not found: ", path)
  message("Reading annotated atlas: ", path)
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("Input must be an annotated Seurat object.")
  meta <- obj[[]]
  required <- c("gender", "anno", "group.a", "orig.ident")
  if (!all(required %in% names(meta))) stop("Missing metadata: ",
    paste(setdiff(required, names(meta)), collapse = ", "))
  expr <- figure5_get_expression(obj)
  if (!inherits(expr, "dgCMatrix")) expr <- methods::as(expr, "dgCMatrix")
  if (anyDuplicated(colnames(expr)) || anyDuplicated(rownames(expr)) ||
      anyDuplicated(rownames(meta))) stop("Duplicated gene or cell identifiers.")
  if (!all(colnames(expr) %in% rownames(meta))) stop("Expression/metadata cell IDs differ.")
  meta <- meta[colnames(expr), , drop = FALSE]
  if (!"ESR2" %in% rownames(expr)) stop("ESR2 is missing from the RNA assay.")
  if (any(!is.finite(expr@x)) || any(expr@x < 0)) stop("Invalid normalized RNA data.")
  meta$gender <- as.character(meta$gender)
  meta$anno <- as.character(meta$anno)
  meta$group.a <- as.character(meta$group.a)
  meta$ESR2 <- as.numeric(expr["ESR2", ])
  # Preserve the historical global median across ALL atlas cells before splitting.
  cutoff <- stats::median(meta$ESR2, na.rm = TRUE)
  meta$ESR2_group <- ifelse(meta$ESR2 > cutoff, "High", "Low")
  keep <- !is.na(meta$gender) & meta$gender %in% c("Female", "Male")
  if (anyNA(meta$anno[keep]) || any(!nzchar(meta$anno[keep]))) stop("Missing cell-type annotations.")
  if (!all(c("Female", "Male") %in% meta$gender[keep])) stop("Both sexes are required.")
  if (!figure5_mpc %in% meta$anno[keep]) stop("Malignant plasma cell label not found.")
  rm(obj)
  invisible(gc())
  list(expr = expr, meta = meta, esr2_cutoff = cutoff,
    input_md5 = unname(tools::md5sum(path)))
}

figure5_write_csv <- function(x, name) {
  utils::write.csv(x, file.path(figure5_output_dir, "source_data", paste0(name, ".csv")),
    row.names = FALSE, na = "")
}

# ======================== 3. Reusable CellChat pipeline ========================
figure5_run_network <- function(expr, meta, cells, key, input_md5,
    filter_small_groups = TRUE) {
  cells <- as.character(cells)
  if (!length(cells)) stop("No cells in network: ", key)
  if (anyDuplicated(cells) || !all(cells %in% colnames(expr))) stop("Invalid cells: ", key)
  cache_file <- file.path(figure5_cache_dir, paste0(key, ".rds"))
  settings <- list(schema = 1L, input_md5 = input_md5,
    cell_ids = digest::digest(cells, algo = "xxhash64"), n_cells = length(cells),
    group_by = "anno", database = "CellChatDB.human", type = "triMean",
    trim = 0.1, raw_use = TRUE, population_size = FALSE, nboot = figure5_nboot,
    seed = figure5_seed, min_cells = if (filter_small_groups) figure5_min_cells else NULL,
    p_threshold = figure5_p_threshold,
    package_versions = vapply(c("CellChat", "SeuratObject", "Matrix"),
      function(x) as.character(utils::packageVersion(x)), character(1)))
  if (!figure5_rebuild_networks && file.exists(cache_file)) {
    saved <- readRDS(cache_file)
    if (identical(saved$settings, settings)) {
      message("Reusing network: ", key)
      figure5_write_csv(saved$edges, paste0("network_", key, "_LR_edges"))
      return(saved)
    }
    message("Cache settings changed; rebuilding: ", key)
  }
  message("Running CellChat: ", key, " (", length(cells), " cells)")
  sub_meta <- meta[cells, , drop = FALSE]
  sub_meta$anno <- droplevels(factor(sub_meta$anno))
  cc <- CellChat::createCellChat(object = expr[, cells, drop = FALSE],
    meta = sub_meta, group.by = "anno")
  cc@DB <- CellChat::CellChatDB.human
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = "triMean", trim = 0.1,
    raw.use = TRUE, population.size = FALSE, nboot = figure5_nboot, seed.use = figure5_seed)
  # Historical D omitted filterCommunication; retain that choice for D only.
  if (filter_small_groups) cc <- CellChat::filterCommunication(cc, min.cells = figure5_min_cells)
  cc <- CellChat::computeCommunProbPathway(cc, thresh = figure5_p_threshold)
  cc <- CellChat::aggregateNet(cc, thresh = figure5_p_threshold)
  edges <- CellChat::subsetCommunication(cc, thresh = figure5_p_threshold)
  if (!nrow(edges)) stop("No supported communication in network: ", key)
  edges$source <- as.character(edges$source)
  edges$target <- as.character(edges$target)
  result <- list(settings = settings, key = key, edges = edges,
    net = cc@net, netP = cc@netP, group_sizes = table(cc@idents))
  saveRDS(result, cache_file)
  figure5_write_csv(edges, paste0("network_", key, "_LR_edges"))
  rm(cc)
  invisible(gc())
  result
}

# ======================== 4. Network panels A/B/D ========================
figure5_network_grob <- function(network, title, esr2_panel = FALSE) {
  # echoGrob draws later, after the sex/ESR2 loop has finished. Freeze labels
  # now so all four delayed plots do not inherit the final loop's title.
  force(title)
  force(esr2_panel)
  mat <- if (esr2_panel) network$net$weight else network$net$count
  weights <- as.numeric(network$group_sizes[rownames(mat)])
  if (anyNA(weights)) stop("Network node names do not match cell counts.")
  if (!length(mat) || any(!is.finite(mat)) || max(mat) <= 0) stop("Empty network: ", title)
  # Preserve CellChat's adjacency, colors, node-size formula and edge-width
  # formula; replace only the circular text layout. Full names remain visible.
  g <- igraph::graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE)
  n <- nrow(mat)
  angle <- 2 * pi * (seq_len(n) - 1L) / n
  coords <- cbind(cos(angle), sin(angle)) * 0.95
  colors <- CellChat::scPalette(n)
  starts <- igraph::ends(g, igraph::E(g), names = FALSE)[, 1]
  edge_weights <- igraph::edge_attr(g, "weight")
  edge_widths <- if (esr2_panel) 0.3 + figure5_d_edge_width_max * edge_weights else
    0.3 + figure5_ab_edge_width_max * edge_weights / max(edge_weights)
  vertex_max <- if (esr2_panel) 8 else if (length(unique(weights)) == 1L) 5 else 15
  vertex_size <- 5 + weights / max(weights) * vertex_max
  labels <- vapply(rownames(mat), function(x) paste(strwrap(x,
    width = figure5_network_label_wrap), collapse = "\n"), character(1))
  label_x <- ifelse(coords[, 1] >= 0, 1.15, -1.15)
  label_y <- numeric(n)
  for (side in c(-1, 1)) {
    idx <- which(sign(label_x) == side)
    idx <- idx[order(coords[idx, 2], decreasing = TRUE)]
    if (length(idx) == 1L) {
      label_y[idx] <- 0
    } else if (length(idx) > 1L) {
      lines <- vapply(strsplit(labels[idx], "\n"), length, integer(1))
      gaps <- (head(lines, -1) + tail(lines, -1)) * 0.10 + 0.05
      # Allocate more separation beside wrapped labels; keep all labels inside
      # the fixed outer band while preserving their vertical order.
      band <- 2 * figure5_network_label_y_extent
      if (sum(gaps) > band) gaps <- gaps * band / sum(gaps) else
        gaps <- gaps + (band - sum(gaps)) / length(gaps)
      label_y[idx] <- figure5_network_label_y_extent - c(0, cumsum(gaps))
    }
  }
  draw <- function() {
    graphics::par(mar = c(0.2, 0.2, 1.8, 0.2), xpd = NA)
    graphics::plot(g, layout = coords, rescale = FALSE,
      xlim = c(-2.3, 2.3), ylim = c(-1.6, 1.6), asp = 1,
      vertex.label = NA, vertex.shape = "circle", vertex.size = vertex_size,
      vertex.color = colors, vertex.frame.color = colors,
      edge.width = edge_widths, edge.color = grDevices::adjustcolor(colors[starts],
        alpha.f = figure5_edge_alpha), edge.curved = 0.2,
      edge.arrow.width = figure5_arrow_width, edge.arrow.size = figure5_arrow_size,
      edge.loop.angle = -atan2(coords[starts, 2], coords[starts, 1]))
    elbow_x <- sign(label_x) * 1.05
    graphics::segments(coords[, 1], coords[, 2], elbow_x, label_y,
      col = "grey65", lwd = figure5_network_leader_width)
    graphics::segments(elbow_x, label_y, label_x * 0.98, label_y,
      col = "grey65", lwd = figure5_network_leader_width)
    for (side in c(-1, 1)) {
      idx <- which(sign(label_x) == side)
      graphics::text(label_x[idx], label_y[idx], labels[idx],
        adj = c(if (side > 0) 0 else 1, 0.5),
        cex = figure5_network_label_cex, family = "sans", col = "black")
    }
    if (esr2_panel) {
      # D now has taller subpanels. Attach each title to its network coordinates
      # so it stays close to the labels instead of the distant device border.
      graphics::text(0, figure5_network_label_y_extent + 0.45, title,
        font = 2, cex = figure5_network_title_cex, family = "sans", col = "black")
    } else graphics::title(main = title, cex.main = figure5_network_title_cex, line = 0.5)
  }
  gridGraphics::echoGrob(draw)
}

figure5_tag <- function(plot, tag) {
  plot + ggplot2::labs(tag = tag) + ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 16, face = "bold"),
    plot.tag.position = "topleft", plot.margin = ggplot2::margin(8, 8, 8, 8))
}

figure5_network_plot <- function(network, title, tag = NULL, esr2_panel = FALSE) {
  p <- patchwork::wrap_elements(full = figure5_network_grob(network, title, esr2_panel),
    clip = FALSE)
  if (!is.null(tag)) p <- figure5_tag(p, tag)
  p
}

# ======================== 5. C: directional plasma-cell strengths ========================
figure5_strength_data <- function(networks, mode = figure5_c_mode) {
  if (!mode %in% c("legacy", "corrected")) stop("C mode must be legacy or corrected.")
  scores <- dplyr::bind_rows(lapply(names(networks), function(sex) {
    networks[[sex]]$edges %>% dplyr::group_by(source, target) %>%
      dplyr::summarise(strength = sum(prob), .groups = "drop") %>%
      dplyr::mutate(gender = sex)
  }))
  if (mode == "legacy") {
    # Exact historical axis choice: incoming edges all share the MPC target.
    # Their bars overlap; their visible height is NOT total incoming strength.
    dat <- scores %>% dplyr::filter(source == figure5_mpc | target == figure5_mpc) %>%
      dplyr::mutate(direction = ifelse(source == figure5_mpc, "sending", "receiving"),
        partner = target)
  } else {
    dat <- dplyr::bind_rows(
      scores %>% dplyr::filter(target == figure5_mpc) %>%
        dplyr::mutate(direction = "receiving", partner = source),
      scores %>% dplyr::filter(source == figure5_mpc) %>%
        dplyr::mutate(direction = "sending", partner = target))
    # Include zero-valued partner categories; MPC->MPC appears once in EACH view.
    types <- sort(unique(unlist(lapply(networks, function(x) names(x$group_sizes)))))
    dat <- tidyr::expand_grid(gender = c("Female", "Male"),
      direction = c("receiving", "sending"), partner = types) %>%
      dplyr::left_join(dplyr::select(dat, gender, direction, partner, strength),
        by = c("gender", "direction", "partner")) %>%
      dplyr::mutate(strength = tidyr::replace_na(strength, 0))
  }
  dat$gender <- factor(dat$gender, levels = c("Female", "Male"))
  dat$direction <- factor(dat$direction, levels = c("receiving", "sending"))
  dat
}

figure5_plot_C <- function(dat) {
  figure5_tag(ggplot2::ggplot(dat, ggplot2::aes(x = partner, y = strength, fill = gender)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::facet_wrap(~direction) +
    ggplot2::scale_fill_manual(values = c(Female = "#F8766D", Male = "#00BFC4"), drop = FALSE) +
    ggplot2::theme_classic(base_size = 9) + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = figure5_c_x_label_angle,
        hjust = 1, vjust = 1, size = 6, margin = ggplot2::margin(t = 4)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
      panel.spacing.x = grid::unit(7, "mm"),
      legend.position = "top", legend.direction = "horizontal",
      legend.key.size = grid::unit(3, "mm")) +
    ggplot2::labs(title = "Plasma cell communication strength",
      x = if (figure5_c_mode == "legacy") "target" else "Partner cell type",
      y = "Communication strength"), "C")
}

# ======================== 6. E (formerly F): pooled-sex LR comparison ========================
figure5_compare_sexes <- function(male, female, absent = c("zero", "drop")) {
  absent <- match.arg(absent)
  keys <- c("ligand", "receptor", "source", "target")
  prepare <- function(edges, suffix) {
    if (!all(c(keys, "prob", "pval") %in% names(edges))) stop("Invalid LR edge table.")
    x <- edges[, c(keys, "prob", "pval"), drop = FALSE]
    x[keys] <- lapply(x[keys], as.character)
    if (anyDuplicated(x[keys])) stop("Duplicate ligand/receptor/source/target keys.")
    names(x)[names(x) %in% c("prob", "pval")] <- paste0(c("prob", "pval"), suffix)
    x
  }
  out <- dplyr::full_join(prepare(male, "_M"), prepare(female, "_F"), by = keys)
  if (absent == "drop") out <- dplyr::filter(out, !is.na(prob_M), !is.na(prob_F))
  out %>% dplyr::mutate(prob_M = tidyr::replace_na(prob_M, 0),
    prob_F = tidyr::replace_na(prob_F, 0), pval_M = tidyr::replace_na(pval_M, 1),
    pval_F = tidyr::replace_na(pval_F, 1),
    log2FC = log2((prob_F + 1e-6) / (prob_M + 1e-6)), pval_min = pmin(pval_M, pval_F),
    LR = paste(ligand, receptor, sep = "\u2013"),
    cellpair = paste(source, target, sep = " \u2192 "))
}

figure5_E_data <- function(comparison) {
  # The historical source omitted the stated MPC/immune restriction. On the
  # recovered data its unmodified top-12 rule selects NO MPC interactions.
  # Complete that missing filter for the submitted panel's stated population.
  # Retain the original all-pair behavior as a separately labelled option.
  if (!figure5_e_scope %in% c("mpc_immune", "historical_all_pairs")) stop("Invalid E scope (formerly F).")
  if (figure5_e_scope == "mpc_immune") {
    comparison <- comparison %>% dplyr::filter(
      (source == figure5_mpc & target %in% figure5_immune_cells) |
        (target == figure5_mpc & source %in% figure5_immune_cells))
  }
  top <- comparison %>% dplyr::filter(pval_min < figure5_p_threshold) %>%
    dplyr::group_by(cellpair) %>%
    dplyr::summarise(meanFC = mean(abs(log2FC)), .groups = "drop") %>%
    dplyr::slice_max(meanFC, n = figure5_e_top_cellpairs) %>% dplyr::pull(cellpair)
  dat <- comparison %>% dplyr::filter(cellpair %in% top,
    pval_min < figure5_p_threshold, abs(log2FC) > figure5_e_min_abs_log2fc) %>%
    dplyr::group_by(LR) %>% dplyr::slice_max(abs(log2FC), n = 1) %>% dplyr::ungroup() %>%
    dplyr::mutate(sex_bias = ifelse(log2FC > 0, "Female-biased", "Male-biased"),
      pval_plot = pmax(pval_min, 1e-300), size_plot = -log10(pval_plot),
      selection_scope = figure5_e_scope)
  if (!nrow(dat)) stop("No interactions pass the historical E filters (formerly F).")
  dat
}

figure5_plot_E <- function(dat) {
  focused <- identical(unique(dat$selection_scope), "mpc_immune")
  title <- if (focused) {
    "Sex-biased ligand\u2013receptor interactions\nMalignant plasma cells vs immune subsets"
  } else "Sex-biased ligand\u2013receptor interactions\nAll cell-type pairs (historical selection)"
  figure5_tag(ggplot2::ggplot(dat, ggplot2::aes(x = cellpair, y = LR)) +
    ggplot2::geom_point(ggplot2::aes(color = log2FC, size = size_plot, shape = sex_bias), alpha = 0.85) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "log2FC (Female / Male)") +
    ggplot2::scale_shape_manual(values = c("Female-biased" = 16, "Male-biased" = 17)) +
    ggplot2::scale_size(range = c(1, 3), name = "-log10(p)") +
    ggplot2::theme_classic(base_size = 9) + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = ggplot2::element_text(size = 6),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 9),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 6),
      legend.title = ggplot2::element_text(size = 7), legend.text = ggplot2::element_text(size = 6),
      legend.key.height = grid::unit(3, "mm"), legend.key.width = grid::unit(3, "mm")) +
    ggplot2::labs(title = title,
      subtitle = if (focused) "Reconstructed MPC/immune selection" else NULL,
      x = "Cell-type pair", y = "Ligand\u2013receptor pair"), "E")
}

# ======================== 7. F (formerly G): within-OR/SOR comparisons ========================
figure5_F_data <- function(comparison) {
  comparison %>% dplyr::filter(
    (source == figure5_mpc & target %in% figure5_immune_cells) |
      (target == figure5_mpc & source %in% figure5_immune_cells)) %>%
    dplyr::filter(pval_min < figure5_p_threshold, log2FC != 0) %>%
    dplyr::mutate(sex_bias = ifelse(log2FC > 0, "Female-biased", "Male-biased"),
      pval_plot = ifelse(pval_min == 0, 1e-10, pval_min), size_plot = -log10(pval_plot)) %>%
    dplyr::group_by(group.a, sex_bias, cellpair) %>%
    dplyr::slice_max(abs(log2FC), n = figure5_f_top_per_cellpair) %>% dplyr::ungroup() %>%
    dplyr::mutate(LR = factor(LR, levels = unique(LR[order(log2FC)])),
      cellpair = factor(cellpair, levels = unique(cellpair)),
      group.a = factor(group.a, levels = c("OR", "SOR")))
}

figure5_plot_F <- function(dat) {
  if (!nrow(dat)) stop("No interactions pass the historical F filters (formerly G).")
  if (!all(c("OR", "SOR") %in% as.character(dat$group.a))) stop("F requires both OR and SOR.")
  figure5_tag(ggplot2::ggplot(dat, ggplot2::aes(x = cellpair, y = LR)) +
    ggplot2::geom_point(ggplot2::aes(color = log2FC, size = size_plot, shape = sex_bias), alpha = 0.85) +
    ggplot2::facet_wrap(~group.a, ncol = 1, scales = "free_y") +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "log2FC\n(Female / Male)") +
    ggplot2::scale_shape_manual(values = c("Female-biased" = 16, "Male-biased" = 17), name = "Sex bias") +
    ggplot2::scale_size(range = c(1, 3), name = "-log10(p)") + ggplot2::coord_flip() +
    ggplot2::theme_classic(base_size = 9) + ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 5),
      axis.text.y = ggplot2::element_text(size = 7,
        margin = ggplot2::margin(r = 3)),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = figure5_f_y_title_margin_pt)),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
      legend.title = ggplot2::element_text(size = 7), legend.text = ggplot2::element_text(size = 6),
      legend.key.height = grid::unit(3, "mm"), legend.key.width = grid::unit(3, "mm")) +
    ggplot2::labs(title = paste0("Sex-biased ligand\u2013receptor interactions between malignant plasma cells\n",
      "and immune subsets (OR vs SOR)"), x = "Cell-type pair", y = "Ligand\u2013receptor pair"), "F") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 12))
}

# ======================== 8. Export devices ========================
figure5_export <- function(plot, stem, width_mm, height_mm, formats = figure5_formats) {
  # Assemble lazy plot arguments before opening the export device. Otherwise
  # grob construction can initialize a page before grid.newpage(), yielding a
  # leading blank PDF page when export is called with figure5_assemble(...).
  force(plot)
  for (format in formats) {
    filename <- paste0(stem, ".", format)
    export_one <- function() {
      switch(format,
        pdf = grDevices::cairo_pdf(filename, width = width_mm / 25.4, height = height_mm / 25.4,
          family = "sans", bg = "white"),
        png = grDevices::png(filename, width = width_mm / 25.4, height = height_mm / 25.4,
          units = "in", res = figure5_dpi, type = "cairo", bg = "white"),
        svg = grDevices::svg(filename, width = width_mm / 25.4, height = height_mm / 25.4,
          family = "sans", bg = "white"),
        tiff = grDevices::tiff(filename, width = width_mm / 25.4, height = height_mm / 25.4,
          units = "in", res = figure5_dpi, type = "cairo", compression = "lzw", bg = "white"),
        stop("Unsupported export format: ", format))
      device <- grDevices::dev.cur()
      on.exit({if (device %in% grDevices::dev.list()) grDevices::dev.off(device)}, add = TRUE)
      if (inherits(plot, c("grob", "gTree", "gtable"))) {
        grid::grid.newpage()
        grid::grid.draw(plot)
      } else print(plot)
    }
    tryCatch(export_one(), error = function(e) stop("Export failed [", filename, "]: ", conditionMessage(e)))
  }
}


# ======================== 9. Load verified networks ========================
figure5_network_keys <- function() {
  c("pooled_Female", "pooled_Male",
    "ESR2_Male_High", "ESR2_Male_Low", "ESR2_Female_High", "ESR2_Female_Low",
    "response_Female_OR", "response_Male_OR", "response_Female_SOR", "response_Male_SOR")
}

figure5_load_networks <- function() {
  keys <- figure5_network_keys()
  if (figure5_use_cached_networks && !figure5_rebuild_networks) {
    # Hash the input once, without loading the 1.7-GB Seurat object into RAM.
    if (figure5_verify_cached_input_md5) {
      if (!file.exists(figure5_input_file)) stop("Input atlas missing: ", figure5_input_file)
      md5 <- unname(tools::md5sum(figure5_input_file))
      if (!identical(md5, figure5_expected_input_md5)) {
        stop("Input atlas differs from the verified caches. Set figure5_use_cached_networks = FALSE to recompute.")
      }
    } else md5 <- figure5_expected_input_md5
    common <- list(schema = 1L, input_md5 = md5, group_by = "anno",
      database = "CellChatDB.human", type = "triMean", trim = 0.1,
      raw_use = TRUE, population_size = FALSE, nboot = figure5_nboot,
      seed = figure5_seed, p_threshold = figure5_p_threshold)
    nets <- stats::setNames(lapply(keys, function(key) {
      path <- file.path(figure5_cache_dir, paste0(key, ".rds"))
      if (!file.exists(path)) stop("Network cache missing: ", path,
        ". Restore it or set figure5_use_cached_networks = FALSE.")
      net <- readRDS(path)
      required <- c("settings", "key", "edges", "net", "group_sizes")
      if (!all(required %in% names(net)) || !identical(net$key, key)) {
        stop("Invalid network cache: ", key)
      }
      expected <- common
      expected["min_cells"] <- list(if (startsWith(key, "ESR2_")) NULL else figure5_min_cells)
      compatible <- all(vapply(names(expected), function(nm) {
        identical(net$settings[[nm]], expected[[nm]])
      }, logical(1)))
      if (!compatible) stop("Cache settings differ: ", key,
        ". Set figure5_use_cached_networks = FALSE to recompute.")
      for (metric in c("count", "weight")) {
        mat <- net$net[[metric]]
        if (!is.matrix(mat) || nrow(mat) != ncol(mat) || any(!is.finite(mat)) ||
            !identical(rownames(mat), colnames(mat)) ||
            !identical(rownames(mat), names(net$group_sizes))) {
          stop("Invalid network matrix: ", key, "/", metric)
        }
      }
      if (sum(net$group_sizes) != net$settings$n_cells) stop("Cell counts differ: ", key)
      figure5_write_csv(net$edges, paste0("network_", key, "_LR_edges"))
      message("Verified cache: ", key)
      net
    }), keys)
    versions <- lapply(nets, function(x) x$settings$package_versions)
    if (!all(vapply(versions, identical, logical(1), versions[[1L]]))) {
      stop("The caches were inferred with inconsistent package versions.")
    }
    return(list(networks = nets, input_md5 = md5,
      mode = "verified_cache", esr2_cutoff = NA_real_))
  }

  input <- figure5_read_input(figure5_input_file)
  meta <- input$meta
  expr <- input$expr
  valid <- !is.na(meta$gender) & meta$gender %in% c("Female", "Male")
  figure5_write_csv(meta[valid, , drop = FALSE] %>%
    dplyr::count(gender, group.a, anno, name = "n_cells"), "atlas_cell_counts")
  figure5_write_csv(unique(meta[valid, c("orig.ident", "gender", "group.a")]),
    "atlas_donor_metadata")
  nets <- list()
  for (sex in c("Female", "Male")) {
    use <- valid & meta$gender == sex
    key <- paste0("pooled_", sex)
    nets[[key]] <- figure5_run_network(expr, meta, rownames(meta)[use], key, input$input_md5)
    for (level in c("High", "Low")) {
      key <- paste0("ESR2_", sex, "_", level)
      cells <- rownames(meta)[use & meta$ESR2_group == level]
      nets[[key]] <- figure5_run_network(expr, meta, cells, key, input$input_md5,
        filter_small_groups = FALSE)
    }
    for (response in c("OR", "SOR")) {
      key <- paste0("response_", sex, "_", response)
      cells <- rownames(meta)[use & !is.na(meta$group.a) & meta$group.a == response]
      nets[[key]] <- figure5_run_network(expr, meta, cells, key, input$input_md5)
    }
  }
  list(networks = nets[keys], input_md5 = input$input_md5,
    mode = "atlas_pipeline", esr2_cutoff = input$esr2_cutoff)
}

# ======================== 10. Six-panel assembly ========================
figure5_assemble <- function(panels) {
  rows <- c(figure5_layout_top_mm, figure5_layout_middle_mm, figure5_layout_bottom_mm)
  if (any(!is.finite(rows)) || any(rows <= 0) ||
      abs(sum(rows) - figure5_height_mm) > 1e-8) {
    stop("Positive row heights must sum to figure5_height_mm.")
  }
  if (!all(LETTERS[1:6] %in% names(panels))) stop("Panels A-F are required.")
  top <- rows[1] / figure5_height_mm
  middle <- rows[2] / figure5_height_mm
  bottom <- rows[3] / figure5_height_mm
  left <- figure5_left_column_mm / figure5_width_mm
  ab <- figure5_ab_height_mm / figure5_height_mm
  if (!is.finite(left) || left <= 0 || left >= 1 ||
      !is.finite(ab) || ab <= 0 || ab > top) stop("Invalid panel width/height.")
  positions <- list(
    A = c(left / 4, 1 - ab / 2, left / 2, ab),
    B = c(3 * left / 4, 1 - ab / 2, left / 2, ab),
    C = c((1 + left) / 2, 1 - top / 2, 1 - left, top),
    D = c(left / 2, bottom + middle / 2, left, middle),
    E = c((1 + left) / 2, bottom + middle / 2, 1 - left, middle),
    F = c(0.5, bottom / 2, 1, bottom))
  children <- lapply(names(positions), function(letter) {
    g <- patchwork::patchworkGrob(patchwork::wrap_plots(list(panels[[letter]])))
    pos <- positions[[letter]]
    grid::grobTree(g, vp = grid::viewport(x = pos[1], y = pos[2],
      width = pos[3], height = pos[4],
      clip = if (letter %in% c("A", "B", "D")) "off" else "on"),
      name = paste0("panel_", letter))
  })
  do.call(grid::grobTree, children)
}

# ======================== 11. Run and export ========================
figure5_main <- function() {
  figure5_check_packages()
  stopifnot(figure5_c_mode %in% c("legacy", "corrected"),
    figure5_e_scope %in% c("mpc_immune", "historical_all_pairs"),
    length(figure5_seed) == 1L, figure5_nboot > 0, figure5_min_cells > 0,
    figure5_p_threshold > 0, figure5_p_threshold < 1,
    figure5_ab_edge_width_max > 0, figure5_d_edge_width_max > 0,
    figure5_edge_alpha > 0, figure5_edge_alpha <= 1,
    figure5_arrow_width > 0, figure5_arrow_size > 0,
    figure5_width_mm > 0, figure5_height_mm > 0, figure5_dpi > 0,
    figure5_c_x_label_angle > 0, figure5_c_x_label_angle < 90,
    figure5_c_standalone_width_mm > 0, figure5_c_standalone_height_mm > 0,
    figure5_f_y_title_margin_pt > 0, figure5_f_standalone_width_mm > 0,
    figure5_e_top_cellpairs > 0, figure5_f_top_per_cellpair > 0,
    figure5_e_min_abs_log2fc >= 0,
    figure5_network_label_cex > 0, figure5_network_label_wrap > 0,
    figure5_network_label_y_extent > 0, figure5_network_label_y_extent < 1.5,
    figure5_network_title_cex > 0, figure5_network_leader_width > 0,
    all(figure5_formats %in% c("pdf", "png", "svg", "tiff")))
  # Validate layout before any potentially expensive inference.
  if (abs(figure5_layout_top_mm + figure5_layout_middle_mm +
      figure5_layout_bottom_mm - figure5_height_mm) > 1e-8) {
    stop("Row heights do not sum to the page height.")
  }
  for (path in c(figure5_output_dir, figure5_cache_dir,
      file.path(figure5_output_dir, "source_data"), dirname(figure5_figure_stem))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan("sequential")
  set.seed(figure5_seed)
  loaded <- figure5_load_networks()
  nets <- loaded$networks
  pooled <- list(Female = nets$pooled_Female, Male = nets$pooled_Male)
  panels <- list(
    A = figure5_network_plot(pooled$Female, "Female Network", "A"),
    B = figure5_network_plot(pooled$Male, "Male Network", "B"))
  C_data <- figure5_strength_data(pooled)
  panels$C <- figure5_plot_C(C_data)
  figure5_write_csv(C_data, paste0("Figure5C_", figure5_c_mode))

  d_panels <- list()
  for (sex in c("Male", "Female")) for (level in c("High", "Low")) {
    key <- paste0("ESR2_", sex, "_", level)
    d_panels[[key]] <- figure5_network_plot(nets[[key]], paste(sex, "ESR2", level),
      esr2_panel = TRUE)
  }
  panels$D <- figure5_tag(patchwork::wrap_elements(full = patchwork::patchworkGrob(
    patchwork::wrap_plots(d_panels, ncol = 2)), clip = FALSE), "D")

  # New E is the former pooled-sex ligand-receptor panel F.
  E_comparison <- figure5_compare_sexes(pooled$Male$edges, pooled$Female$edges, absent = "zero")
  E_data <- figure5_E_data(E_comparison)
  panels$E <- figure5_plot_E(E_data)
  figure5_write_csv(E_comparison, "Figure5E_all_LR_comparisons")
  figure5_write_csv(E_data, "Figure5E_plotted_interactions")

  # New F is the former OR/SOR ligand-receptor panel G.
  F_comparison <- dplyr::bind_rows(lapply(c("OR", "SOR"), function(response) {
    figure5_compare_sexes(nets[[paste0("response_Male_", response)]]$edges,
      nets[[paste0("response_Female_", response)]]$edges, absent = "drop") %>%
      dplyr::mutate(group.a = response)
  }))
  F_data <- figure5_F_data(F_comparison)
  panels$F <- figure5_plot_F(F_data)
  figure5_write_csv(F_comparison, "Figure5F_all_LR_comparisons")
  figure5_write_csv(F_data, "Figure5F_plotted_interactions")
  panels <- panels[LETTERS[1:6]]
  figure <- figure5_assemble(panels)

  if (figure5_export_panels) {
    dir.create(file.path(figure5_output_dir, "panels"), showWarnings = FALSE)
    sizes <- list(A = c(100, 85), B = c(100, 85),
      C = c(figure5_c_standalone_width_mm, figure5_c_standalone_height_mm),
      D = c(200, 170), E = c(160, 175), F = c(figure5_f_standalone_width_mm, 150))
    for (letter in names(panels)) {
      figure5_export(panels[[letter]],
        file.path(figure5_output_dir, "panels", paste0("Figure5", letter)),
        sizes[[letter]][1], sizes[[letter]][2])
    }
  }
  figure5_export(figure, figure5_figure_stem, figure5_width_mm, figure5_height_mm)
  mapping <- data.frame(new_panel = LETTERS[1:6], old_panel = c("A", "B", "C", "D", "F", "G"),
    content = c("Female count network", "Male count network",
      "Plasma-cell communication strength", "Sex-by-ESR2 weight networks",
      "Pooled-sex ligand-receptor interactions", "OR/SOR ligand-receptor interactions"))
  figure5_write_csv(mapping, "Figure5_panel_mapping")
  manifest <- dplyr::bind_rows(lapply(nets, function(net) {
    data.frame(key = net$key, input_md5 = net$settings$input_md5,
      cell_ids_xxhash64 = net$settings$cell_ids, n_cells = net$settings$n_cells,
      n_LR_edges = nrow(net$edges), cache_md5 = unname(tools::md5sum(
        file.path(figure5_cache_dir, paste0(net$key, ".rds")))))
  }))
  figure5_write_csv(manifest, "Figure5_network_manifest")
  provenance <- c(
    "Figure5.new.R: six panels after removing the former E pathway panel.",
    "Mapping: new A/B/C/D/E/F = old A/B/C/D/F/G; old E is omitted entirely.",
    "No PDF import, raster reference panel or original E pathway calculation is used.",
    paste("Input:", figure5_input_file), paste("Input MD5:", loaded$input_md5),
    paste("Network loading mode:", loaded$mode),
    paste("Cached-input MD5 verification requested:", figure5_verify_cached_input_md5),
    paste("Global ESR2 median:", if (is.na(loaded$esr2_cutoff))
      "defined in cached analysis; not recalculated during replot" else loaded$esr2_cutoff),
    paste("C mode:", figure5_c_mode),
    paste("E (formerly F) selection scope:", figure5_e_scope),
    "E retains the reconstructed MPC/immune focus filter; exact historical filter source remains unknown.",
    "F (formerly G): absent-sex edges dropped; historical top-5-per-response/sex-bias/cellpair rule, ties retained.",
    "P values are within-network CellChat interaction support, not between-sex tests.",
    "A-C/E pool HD/MM cells by sex; D uses the global ESR2 split; F uses OR/SOR cells only.",
    "Network matrices, edge sets and C/E/F comparison/selection methods are retained.",
    paste("Page size (mm):", figure5_width_mm, "x", figure5_height_mm),
    paste("Top/middle/bottom row heights (mm):", figure5_layout_top_mm,
      figure5_layout_middle_mm, figure5_layout_bottom_mm),
    paste("Left-column width (mm):", figure5_left_column_mm),
    paste("C x-label angle (degrees):", figure5_c_x_label_angle),
    paste("F y-axis title right margin (pt):", figure5_f_y_title_margin_pt),
    paste("C/E/F plotted records:", nrow(C_data), nrow(E_data), nrow(F_data)),
    paste("Seed:", figure5_seed, "nboot:", figure5_nboot, "min.cells:", figure5_min_cells),
    "", utils::capture.output(utils::sessionInfo()))
  writeLines(provenance, file.path(figure5_output_dir, "Figure5.new_provenance.txt"))
  message("Completed: ", figure5_figure_stem, ".pdf")
  invisible(list(figure = figure, panels = panels,
    source_data = list(C = C_data, E = E_data, F = F_data)))
}

if (!identical(Sys.getenv("FIGURE5_RUN", unset = "1"), "0")) {
  figure5_result <- figure5_main()
}
