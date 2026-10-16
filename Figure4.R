# Figure 4
# ======================== 1. Parameters ========================
# Dataset audit 2026-10-14: identical cells and ESR2 counts, but the 2026-05-16 UMAP
# in MM18.MM.anno matches the historical H-J; the earlier 2026-04-21 UMAP does not.
figure4_analysis_file <- "/Users/tf/e/ndmm/HM38/scRNA_harmony_MM18.MM.anno.rds"
figure4_cyto_anno_file <- "/Users/tf/e/ndmm/HM38/CytoTRACE/CytoTRACE_plot_table_anno3.txt"
figure4_cyto_response_file <- "/Users/tf/e/ndmm/HM38/CytoTRACE/CytoTRACE_plot_table_group.a.txt"
figure4_cyto_gender_file <- "/Users/tf/e/ndmm/HM38/CytoTRACE_plot_table_gender.txt"
figure4_output_dir <- "/Users/tf/e/ndmm/ESR2/submit/Communications Biology/Figure4_code/outputs_manual"
figure4_formats <- c("pdf", "png")  # Optional: "svg", "tiff".
figure4_width_mm <- 183
figure4_height_mm <- 274
figure4_dpi <- 300
figure4_font_increase_pt <- 2  # Add 2 pt to every existing A-M text size.
figure4_seed <- 20261014L
figure4_cores <- 4L

# These three positions were read only from original H, before inspecting K-M.
# They are newly selected approximate manual roots, not recovered original IDs.
# Set root_positions to NULL to select interactively; supply only one root input.
figure4_root_mode <- "interactive"
figure4_root_cells <- NULL       # Example: c("MM012_AAACCCAGTAGAAACT-1").
figure4_root_pr_nodes <- NULL    # Example: c("Y_1", "Y_50"); graph-dependent IDs.
figure4_root_positions <- data.frame(
  reference_label = c("1", "2", "3"),
  UMAP_1 = c(5.75, 0.70, 0.95), UMAP_2 = c(-6.95, -5.45, 9.75))
# Reuse an already learned graph only after validating its inputs and settings.
figure4_graph_cache_file <- "/Users/tf/e/ndmm/ESR2/submit/Communications Biology/Figure4_code/outputs_auto_corrected/cache/Figure4_trajectory.rds"
figure4_cache_file <- file.path(figure4_output_dir, "cache", "Figure4_trajectory.rds")
figure4_rebuild_trajectory <- FALSE

# ======================== 2. Input and export helpers ========================
figure4_check_packages <- function() {
  required <- c("SeuratObject", "Seurat", "Matrix", "ggplot2", "patchwork",
                "monocle3", "SingleCellExperiment", "SummarizedExperiment", "igraph", "RColorBrewer")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". On this machine use /Library/Frameworks/R.framework/Resources/bin/Rscript.")
  if (!capabilities("cairo")) stop("Cairo graphics support is required for export.")
}

figure4_read_cyto <- function(path) {
  if (!file.exists(path)) stop("CytoTRACE source table not found: ", path)
  # Historical tables have an unnamed leading cell-ID column.
  x <- utils::read.delim(path, row.names = 1, check.names = FALSE,
                         stringsAsFactors = FALSE)
  required <- c("CytoTRACE", "Component1", "Component2", "Phenotype")
  if (!all(required %in% names(x))) stop("Invalid CytoTRACE table: ", path)
  if (anyDuplicated(rownames(x)) || anyNA(x[, required]) ||
      !all(vapply(x[, required[1:3]], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(x[, required[1:3]]))) ||
      any(x$CytoTRACE < 0 | x$CytoTRACE > 1)) stop("Invalid CytoTRACE values: ", path)
  x$cell <- rownames(x)
  x
}

figure4_get_counts <- function(obj) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
  } else {
    SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts")
  }
}

figure4_export <- function(plot, stem, width_mm, height_mm, formats, dpi) {
  for (format in formats) {
    filename <- paste0(stem, ".", format)
    export_one <- function() {
      # Built-in devices avoid the ragg/svglite/systemfonts lookup issue observed
      # in this project's interactive R sessions. PDF text remains selectable.
      switch(format,
             pdf = grDevices::cairo_pdf(filename, width = width_mm / 25.4,
                                        height = height_mm / 25.4, family = "sans", bg = "white"),
             png = grDevices::png(filename, width = width_mm / 25.4,
                                  height = height_mm / 25.4, units = "in", res = dpi,
                                  type = "cairo", family = "sans", bg = "white"),
             svg = grDevices::svg(filename, width = width_mm / 25.4,
                                  height = height_mm / 25.4, family = "sans", bg = "white"),
             tiff = grDevices::tiff(filename, width = width_mm / 25.4,
                                    height = height_mm / 25.4, units = "in", res = dpi,
                                    type = "cairo", compression = "lzw", family = "sans", bg = "white")
      )
      device <- grDevices::dev.cur()
      on.exit({
        if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
      }, add = TRUE)
      print(plot)
    }
    tryCatch(export_one(), error = function(e) {
      stop("Export failed [", filename, "]: ", conditionMessage(e), call. = FALSE)
    })
  }
}

# ======================== 3. Monocle graph and explicit rooting ========================
figure4_prepare_trajectory <- function(obj, meta, cyto_scores, root_mode,
                                       root_cells, root_pr_nodes, cache_file, rebuild, provenance, seed, cores,
                                       root_positions = NULL, graph_cache_file = NULL) {
  if (!root_mode %in% c("interactive", "cytotrace")) {
    stop("root_mode must be 'interactive' or 'cytotrace'.")
  }
  if (sum(c(!is.null(root_cells), !is.null(root_pr_nodes), !is.null(root_positions))) > 1L) {
    stop("Supply only one of root_cells, root_pr_nodes, or root_positions.")
  }
  if (!is.null(root_positions) && (!is.data.frame(root_positions) ||
                                   !all(c("reference_label", "UMAP_1", "UMAP_2") %in% names(root_positions)) ||
                                   nrow(root_positions) < 1L || anyDuplicated(root_positions$reference_label) ||
                                   any(!is.finite(as.matrix(root_positions[, c("UMAP_1", "UMAP_2")]))))) {
    stop("root_positions must contain unique labels and finite UMAP_1/UMAP_2 coordinates.")
  }
  settings <- list(input = provenance, seed = seed, cores = cores,
                   num_dim = 100L, root_mode = root_mode, root_cells = root_cells,
                   root_pr_nodes = root_pr_nodes, root_positions = root_positions, method_version = 1L,
                   package_versions = vapply(c("monocle3", "SeuratObject", "igraph", "batchelor"),
                                             function(p) as.character(utils::packageVersion(p)), character(1)))
  settings$graph_cache_signature <- if (is.null(graph_cache_file)) NULL else list(
    path = normalizePath(graph_cache_file, mustWork = FALSE),
    md5 = if (file.exists(graph_cache_file)) unname(tools::md5sum(graph_cache_file)) else NA_character_)
  if (!rebuild && !is.null(cache_file) && file.exists(cache_file)) {
    saved <- readRDS(cache_file)
    if (is.list(saved) && identical(saved$settings, settings) &&
        inherits(saved$cds, "cell_data_set")) {
      message("Reusing the saved trajectory and recorded roots: ", cache_file)
      return(saved)
    }
    message("Trajectory cache settings differ; rebuilding the graph.")
  }
  if (root_mode == "interactive" && is.null(root_cells) &&
      is.null(root_pr_nodes) && is.null(root_positions) && !interactive()) {
    stop("Original roots were not saved. Run source() in interactive R to select roots, ",
         "supply figure4_root_cells / figure4_root_pr_nodes, or explicitly select ",
         "Figure4_auto.R for a new reproducible automatic ordering.")
  }
  set.seed(seed)
  umap <- SeuratObject::Embeddings(obj, reduction = "umap")[rownames(meta), 1:2, drop = FALSE]
  cds <- NULL
  graph_source <- data.frame(path = NA_character_, md5 = NA_character_, reuse = "newly learned graph")
  if (!rebuild && !is.null(graph_cache_file) && file.exists(graph_cache_file)) {
    candidate <- readRDS(graph_cache_file)
    graph_fields <- c("input", "seed", "cores", "num_dim", "method_version", "package_versions")
    compatible <- is.list(candidate) && inherits(candidate$cds, "cell_data_set") &&
      identical(candidate$settings[graph_fields], settings[graph_fields]) &&
      identical(colnames(candidate$cds), rownames(meta))
    if (compatible) {
      cached_meta <- as.data.frame(SummarizedExperiment::colData(candidate$cds))
      compatible <- all(vapply(names(meta), function(nm) {
        nm %in% names(cached_meta) && identical(as.character(cached_meta[[nm]]), meta[[nm]])
      }, logical(1))) && isTRUE(all.equal(
        unname(SingleCellExperiment::reducedDim(candidate$cds, "UMAP")),
        unname(umap), tolerance = 1e-12))
    }
    if (compatible) {
      cds <- candidate$cds
      graph_source <- data.frame(path = normalizePath(graph_cache_file),
                                 md5 = unname(tools::md5sum(graph_cache_file)), reuse = "learned graph only; roots replaced")
      message("Reusing the validated learned graph; replacing its roots and pseudotime.")
    } else message("Graph cache does not match the current input/settings; learning a new graph.")
    rm(candidate)
  }
  if (is.null(cds)) {
    counts <- figure4_get_counts(obj)
    counts <- counts[, rownames(meta), drop = FALSE]
    genes <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
    cds <- monocle3::new_cell_data_set(counts, cell_metadata = meta, gene_metadata = genes)
    message("Preprocessing the Monocle 3 object...")
    cds <- monocle3::preprocess_cds(cds, num_dim = 100L)
    cds <- monocle3::align_cds(cds, num_dim = 100L, alignment_group = "orig.ident")
    SingleCellExperiment::reducedDim(cds, "UMAP") <- umap
    cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP",
                                   random_seed = seed, nn_control = list(cores = cores))
    message("Learning the trajectory graph on the saved UMAP...")
    cds <- monocle3::learn_graph(cds, use_partition = TRUE,
                                 learn_graph_control = list(nn.cores = cores))
  }
  graph <- monocle3::principal_graph(cds)[["UMAP"]]
  closest <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)
  node_names <- igraph::V(graph)$name
  # Monocle stores nearest principal-vertex IDs as numbers and names them Y_<ID>.
  nearest_node <- paste0("Y_", as.integer(closest[colnames(cds), 1]))
  if (any(!nearest_node %in% node_names)) stop("Unknown nearest principal-vertex IDs.")
  names(nearest_node) <- colnames(cds)
  membership <- igraph::components(graph)$membership
  cell_component <- unname(membership[nearest_node])
  names(cell_component) <- colnames(cds)
  if (anyNA(cell_component)) stop("Could not map cells to trajectory graph components.")
  if (!is.null(root_positions)) {
    vertices <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
    root_pr_nodes <- vapply(seq_len(nrow(root_positions)), function(i) {
      distances <- (vertices[, 1] - root_positions$UMAP_1[i])^2 +
        (vertices[, 2] - root_positions$UMAP_2[i])^2
      rownames(vertices)[which.min(distances)]
    }, character(1))
    if (anyDuplicated(root_pr_nodes)) stop("Reference positions map to the same principal node.")
    root_distance <- sqrt((vertices[root_pr_nodes, 1] - root_positions$UMAP_1)^2 +
                            (vertices[root_pr_nodes, 2] - root_positions$UMAP_2)^2)
    if (any(root_distance > 0.3)) stop("A reference root is >0.3 UMAP units from the learned graph; select roots interactively.")
    cds <- monocle3::order_cells(cds, root_pr_nodes = root_pr_nodes)
    root_method <- "manual reference-guided approximate positions from original H"
  } else if (!is.null(root_cells)) {
    if (!length(root_cells) || any(!root_cells %in% colnames(cds))) stop("Unknown/empty root_cells.")
    cds <- monocle3::order_cells(cds, root_cells = root_cells)
    root_method <- "user-specified cells"
  } else if (!is.null(root_pr_nodes)) {
    if (!length(root_pr_nodes) || any(!root_pr_nodes %in% node_names)) stop("Unknown/empty root_pr_nodes.")
    cds <- monocle3::order_cells(cds, root_pr_nodes = root_pr_nodes)
    root_method <- "user-specified principal nodes"
  } else if (root_mode == "cytotrace") {
    # One highest-score cell per connected component ensures each disconnected
    # trajectory has a root; ties are resolved lexicographically by cell ID.
    root_cells <- vapply(sort(unique(cell_component)), function(component) {
      cells <- names(cell_component)[cell_component == component]
      cells[order(-cyto_scores[cells], cells)][1]
    }, character(1))
    cds <- monocle3::order_cells(cds, root_cells = root_cells)
    root_method <- "highest saved CytoTRACE cell per connected graph component"
  } else {
    message("Select the earliest root node(s) in Monocle's interactive window.")
    cds <- monocle3::order_cells(cds)
    root_method <- "interactive manual principal-node selection"
  }
  roots <- cds@principal_graph_aux[["UMAP"]]$root_pr_nodes
  if (!length(roots)) stop("No trajectory roots were selected.")
  root_table <- data.frame(root_pr_node = roots,
                           component = unname(membership[roots]), method = root_method,
                           stringsAsFactors = FALSE)
  if (!is.null(root_cells)) {
    root_table$root_cell <- vapply(roots, function(node) {
      paste(root_cells[nearest_node[root_cells] == node], collapse = ";")
    }, character(1))
  } else root_table$root_cell <- NA_character_
  vertices <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  root_table$node_UMAP_1 <- vertices[roots, 1]
  root_table$node_UMAP_2 <- vertices[roots, 2]
  if (!is.null(root_positions)) {
    idx <- match(roots, root_pr_nodes)
    root_table$reference_label <- root_positions$reference_label[idx]
    root_table$requested_UMAP_1 <- root_positions$UMAP_1[idx]
    root_table$requested_UMAP_2 <- root_positions$UMAP_2[idx]
    root_table$distance_to_requested <- root_distance[idx]
  }
  time <- monocle3::pseudotime(cds)
  if (!any(is.finite(time))) stop("All pseudotime values are non-finite; check trajectory roots.")
  if (any(!is.finite(time))) warning(sum(!is.finite(time)),
                                     " cells are unreachable from the chosen roots and are excluded from K-M.")
  saved <- list(cds = cds, settings = settings, roots = root_table, graph_source = graph_source)
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(saved, cache_file)
  }
  saved
}

# ======================== 4. Complete A-M figure ========================
build_figure4 <- function(analysis_file, cyto_anno_file, cyto_response_file,
                          cyto_gender_file, output_dir, formats = c("pdf", "png"),
                          width_mm = 183, height_mm = 274, dpi = 300, seed = 20261014L, cores = 4L,
                          root_mode = "interactive", root_cells = NULL, root_pr_nodes = NULL,
                          cache_file = NULL, rebuild_trajectory = FALSE,
                          root_positions = NULL, graph_cache_file = NULL,
                          font_increase_pt = 2) {
  figure4_check_packages()
  if (!length(formats) || !all(formats %in% c("pdf", "png", "svg", "tiff"))) {
    stop("formats must contain pdf, png, svg, and/or tiff.")
  }
  if (any(!is.finite(c(width_mm, height_mm, dpi, cores))) ||
      any(c(width_mm, height_mm, dpi, cores) <= 0)) stop("Invalid export dimensions or cores.")
  if (length(font_increase_pt) != 1L || !is.finite(font_increase_pt) || font_increase_pt < 0) {
    stop("font_increase_pt must be one finite, nonnegative point size.")
  }
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  input_files <- c(Seurat = analysis_file, CytoTRACE_anno3 = cyto_anno_file,
                   CytoTRACE_response = cyto_response_file, CytoTRACE_gender = cyto_gender_file)
  if (any(!file.exists(input_files))) stop("Input file not found: ",
                                           paste(input_files[!file.exists(input_files)], collapse = ", "))
  message("Reading the annotated 18-patient malignant plasma-cell object...")
  obj <- readRDS(analysis_file)
  if (!inherits(obj, "Seurat")) stop("analysis_file must contain a Seurat object.")
  required_meta <- c("orig.ident", "anno3", "gender", "group.a")
  meta <- obj[[]]
  if (!all(required_meta %in% names(meta))) stop("Required metadata columns are missing.")
  if (!"umap" %in% SeuratObject::Reductions(obj)) stop("Input has no saved UMAP.")
  meta <- meta[, required_meta, drop = FALSE]
  meta[] <- lapply(meta, as.character)
  if (anyNA(meta) || any(as.matrix(meta) == "")) stop("Required metadata contain missing values.")
  anno_order <- c("PBM_CCND1", "MM_IL5RA", "PBM", "MM_ANXA1", "MM",
                  "MM_IGLL5", "MM_CRIP1", "MM_AZGP1", "MM_IGHGP")
  if (!setequal(unique(meta$anno3), anno_order) ||
      !setequal(unique(meta$gender), c("Female", "Male")) ||
      !setequal(unique(meta$group.a), c("OR", "SOR"))) stop("Unexpected subcluster, sex, or response labels.")
  patients <- unique(meta[, c("orig.ident", "gender", "group.a")])
  if (anyDuplicated(patients$orig.ident)) stop("A patient has inconsistent sex/response labels.")
  anno <- figure4_read_cyto(cyto_anno_file)
  response <- figure4_read_cyto(cyto_response_file)
  gender <- figure4_read_cyto(cyto_gender_file)
  if (!setequal(rownames(anno), rownames(meta))) stop("CytoTRACE and Seurat cell sets differ.")
  anno <- anno[rownames(meta), , drop = FALSE]
  for (nm in c("response", "gender")) {
    tab <- get(nm)
    if (!setequal(rownames(tab), rownames(meta))) stop("CytoTRACE cell sets differ: ", nm)
    tab <- tab[rownames(meta), , drop = FALSE]
    if (!isTRUE(all.equal(as.matrix(tab[, 1:3]), as.matrix(anno[, 1:3]), tolerance = 1e-12))) {
      stop("Historical CytoTRACE scores/coordinates differ across grouping tables: ", nm)
    }
    assign(nm, tab)
  }
  if (!identical(as.character(anno$Phenotype), meta$anno3) ||
      !identical(as.character(response$Phenotype), meta$group.a) ||
      !identical(as.character(gender$Phenotype), meta$gender)) stop("CytoTRACE phenotype labels differ from metadata.")
  plot_data <- data.frame(cell = rownames(meta), meta,
                          CytoTRACE = anno$CytoTRACE, Component1 = anno$Component1, Component2 = anno$Component2,
                          row.names = NULL, check.names = FALSE)
  cyto_scores <- stats::setNames(anno$CytoTRACE, rownames(meta))
  manifest <- data.frame(role = names(input_files), path = normalizePath(input_files),
                         md5 = unname(tools::md5sum(input_files)), row.names = NULL)
  provenance <- manifest[, c("role", "md5")]
  source_dir <- file.path(output_dir, "source_data")
  panel_dir <- file.path(output_dir, "panels")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  write_source <- function(x, name) utils::write.csv(x, file.path(source_dir, name), row.names = FALSE)
  
  # A-G use the original CytoTRACE components, not the Seurat UMAP used in H-J.
  font_pt <- function(original_pt) original_pt + font_increase_pt
  theme4 <- function() theme_classic(base_size = font_pt(6.8), base_family = "sans") + theme(
    axis.text = element_text(size = font_pt(6), colour = "black"),
    axis.title = element_text(size = font_pt(7)),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.25),
    legend.text = element_text(size = font_pt(5.5)), legend.title = element_text(size = font_pt(6)),
    legend.key.size = grid::unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(1.5, "mm"),
    legend.spacing.x = grid::unit(1, "mm"),
    strip.background = element_blank(), strip.text = element_text(size = font_pt(7)),
    panel.spacing.x = grid::unit(3.5, "mm"),
    plot.title = element_text(size = font_pt(8), hjust = 0.5),
    plot.tag = element_text(size = font_pt(12)), plot.tag.position = c(0, 1),
    plot.margin = margin(6, 6, 6, 6))
  spectral <- RColorBrewer::brewer.pal(11, "Spectral")
  spectral[6] <- "gold"  # Same palette as the historical plotCytoTRACE().
  cyto_palette <- grDevices::colorRampPalette(spectral)
  anno_colors <- stats::setNames(cyto_palette(length(anno_order)), anno_order)
  response_colors <- c(SOR = "#9E0142", OR = "#5E4FA2")
  gender_colors <- c(Male = "#9E0142", Female = "#5E4FA2")
  curve_colors <- c(OR = "#F8766D", SOR = "#00BFC4")
  cyto_colors <- rev(cyto_palette(50))
  phenotype_plot <- function(column, levels, colors, tag) {
    dat <- plot_data
    dat$Phenotype <- factor(dat[[column]], levels = levels)
    ggplot(dat, aes(Component1, Component2, colour = Phenotype)) +
      geom_point(size = 0.12) + scale_colour_manual(values = colors, drop = FALSE) +
      coord_fixed() + labs(title = "Phenotype", tag = tag, colour = NULL) +
      theme4() + guides(colour = guide_legend(
        nrow = if (column == "anno3") NULL else 1,
        override.aes = list(size = 1.6))) +
      theme(legend.position = if (column == "anno3") "right" else "bottom")
  }
  cyto_boxplot <- function(column, levels, colors, tag) {
    dat <- plot_data
    # plotCytoTRACE orders boxplots by the median and phenotype maps by the mean.
    # In the nine-subcluster plot AZGP1/CRIP1 therefore exchange positions.
    levels <- names(sort(tapply(dat$CytoTRACE, dat[[column]], stats::median), decreasing = TRUE))
    if (column == "anno3") colors <- stats::setNames(cyto_palette(length(levels)), levels)
    dat$Phenotype <- factor(dat[[column]], levels = levels)
    if (column != "anno3") {
      # E/G: all cells contribute to the bounded density and the box statistics.
      # A light violin replaces the opaque stack of jittered observations;
      # the white box keeps the median and interquartile range visible.
      return(ggplot(dat, aes(Phenotype, CytoTRACE, colour = Phenotype, fill = Phenotype)) +
               geom_violin(width = 0.82, scale = "width", trim = TRUE, bounds = c(0, 1),
                           alpha = 0.20, linewidth = 0.35) +
               geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA,
                            linewidth = 0.45, median.linewidth = 0.75,
                            whisker.linewidth = 0.35, staplewidth = 0.55) +
               scale_colour_manual(values = colors, drop = FALSE) +
               scale_fill_manual(values = colors, drop = FALSE) +
               scale_x_discrete(expand = expansion(add = 0.55)) +
               scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                                  expand = expansion(mult = c(0.025, 0.045))) +
               labs(x = NULL, y = "CytoTRACE score", tag = tag) + theme4() +
               theme(legend.position = "none",
                     panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2),
                     axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.text.x = element_text(angle = 0, hjust = 0.5, margin = margin(t = 5)),
                     plot.margin = margin(6, 7, 7, 6)))
    }
    # C uses the same all-cell distribution view as E/G. Neutral box outlines
    # keep the median and quartiles readable in the pale yellow/green groups.
    ggplot(dat, aes(Phenotype, CytoTRACE, fill = Phenotype)) +
      geom_violin(width = 0.82, scale = "width", trim = TRUE, bounds = c(0, 1),
                  alpha = 0.35, colour = "grey55", linewidth = 0.25) +
      geom_boxplot(width = 0.20, fill = "white", colour = "grey25", outlier.shape = NA,
                   linewidth = 0.35, median.linewidth = 0.60,
                   whisker.linewidth = 0.25, staplewidth = 0.50) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_x_discrete(expand = expansion(add = 0.55)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                         expand = expansion(mult = c(0.025, 0.045))) +
      labs(x = "Cell phenotypes", y = "CytoTRACE score", tag = tag) +
      theme4() + theme(legend.position = "none",
                       panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2),
                       axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                       axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 4)),
                       plot.margin = margin(6, 7, 7, 6))
  }
  panels <- list(
    A = phenotype_plot("anno3", anno_order, anno_colors, "A"),
    B = ggplot(plot_data[order(plot_data$CytoTRACE), ], aes(Component1, Component2, colour = CytoTRACE)) +
      geom_point(size = 0.12) + scale_colour_gradientn(colours = cyto_colors,
                                                       limits = c(0, 1), breaks = c(0, 0.5, 1),
                                                       labels = c("0\n(More diff.)", "0.5", "1\n(Less diff.)"),
                                                       name = "Predicted\norder") + coord_fixed() +
      labs(title = "CytoTRACE", tag = "B") + theme4() +
      guides(colour = guide_colourbar(theme = theme(
        legend.key.height = grid::unit(18, "mm"),
        legend.title = element_text(margin = margin(b = 8))))),
    C = cyto_boxplot("anno3", anno_order, anno_colors, "C"),
    D = phenotype_plot("group.a", c("SOR", "OR"), response_colors, "D"),
    E = cyto_boxplot("group.a", c("SOR", "OR"), response_colors, "E"),
    F = phenotype_plot("gender", c("Male", "Female"), gender_colors, "F"),
    G = cyto_boxplot("gender", c("Male", "Female"), gender_colors, "G")
  )
  
  # H-J share exactly the same CDS, graph, roots, UMAP and pseudotime scale.
  trajectory <- figure4_prepare_trajectory(obj, meta, cyto_scores,
                                           root_mode, root_cells, root_pr_nodes, cache_file, rebuild_trajectory,
                                           provenance, seed, cores, root_positions, graph_cache_file)
  cds <- trajectory$cds
  rm(obj)
  trajectory_plot <- function(tag, facet = NULL) {
    p <- monocle3::plot_cells(cds, color_cells_by = "pseudotime",
                              label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
                              # Monocle's text-size unit is mm; convert the added pt.
                              label_roots = TRUE, show_trajectory_graph = TRUE,
                              graph_label_size = 2 + font_increase_pt / ggplot2::.pt,
                              cell_size = 0.12, rasterize = FALSE) +
      labs(tag = tag, x = "UMAP 1", y = "UMAP 2") + theme4() +
      scale_x_continuous(breaks = c(-10, 0, 10))
    if (!is.null(facet)) {
      p <- p + facet_wrap(stats::reformulate(facet), nrow = 1) +
        guides(colour = guide_colourbar(direction = "horizontal", theme = theme(
          legend.key.width = grid::unit(36, "mm"),
          legend.key.height = grid::unit(2.5, "mm"), legend.title.position = "top"))) +
        theme(legend.position = "bottom")
    }
    p
  }
  # Factor levels control facet order; they do not alter the inferred trajectory.
  SummarizedExperiment::colData(cds)$gender <- factor(meta$gender, levels = c("Female", "Male"))
  SummarizedExperiment::colData(cds)$group.a <- factor(meta$group.a, levels = c("OR", "SOR"))
  panels$H <- trajectory_plot("H")
  panels$I <- trajectory_plot("I", "gender")
  panels$J <- trajectory_plot("J", "group.a")
  
  # K-M: raw counts, same CDS ordering, OR/SOR, LOESS with original default span.
  rna_counts <- SummarizedExperiment::assay(cds, "counts")
  if (!"ESR2" %in% rownames(rna_counts)) stop("ESR2 is missing from RNA counts.")
  umap <- SingleCellExperiment::reducedDim(cds, "UMAP")
  time <- monocle3::pseudotime(cds)[colnames(cds)]
  trajectory_data <- data.frame(cell = colnames(cds), meta[colnames(cds), , drop = FALSE],
                                UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], pseudotime = as.numeric(time),
                                ESR2 = as.numeric(rna_counts["ESR2", ]), CytoTRACE = cyto_scores[colnames(cds)],
                                finite_pseudotime = is.finite(time), row.names = NULL)
  curve_data <- trajectory_data[trajectory_data$finite_pseudotime &
                                  is.finite(trajectory_data$ESR2), , drop = FALSE]
  if (any(curve_data$ESR2 < 0)) stop("Raw ESR2 counts cannot be negative.")
  curve_plot <- function(dat, tag, title = NULL, legend = FALSE) {
    n_unique <- tapply(dat$pseudotime, dat$group.a, function(v) length(unique(v)))
    if (!all(c("OR", "SOR") %in% names(n_unique)) || any(n_unique < 4)) {
      stop("Insufficient finite pseudotime values for both groups in panel ", tag)
    }
    # Literal original LOESS call; verified defaults are y~x, span=.75, degree=2, n=80.
    ggplot(dat, aes(pseudotime, ESR2, colour = group.a)) +
      geom_smooth(method = "loess", se = FALSE, linewidth = 1.3) +
      scale_colour_manual(values = curve_colors, breaks = c("OR", "SOR")) +
      labs(x = if (tag == "K") "Pseudotime" else "pseudotime",
           y = if (tag == "K") "ESR2 expression" else "ESR2",
           colour = NULL, title = title, tag = tag) +
      theme4() + guides(colour = guide_legend(nrow = 1)) +
      theme(legend.position = if (legend) "bottom" else "none")
  }
  panels$K <- curve_plot(curve_data, "K")
  panels$L <- curve_plot(curve_data[curve_data$gender == "Female", ], "L", "Female")
  panels$M <- curve_plot(curve_data[curve_data$gender == "Male", ], "M", "Male", legend = TRUE)
  panels <- panels[LETTERS[1:13]]
  
  # Keep the original reading order: A/B left; C/D/E right; F/G/H; I/J; K/L/M.
  # Separate row layouts keep fixed-aspect CytoTRACE maps from changing the
  # widths of the unrelated boxplots and ESR2 curves in other rows.
  rows <- list(
    patchwork::wrap_plots(panels[c("A", "C")], nrow = 1, widths = c(5.25, 6.75)),
    patchwork::wrap_plots(panels[c("B", "D", "E")], nrow = 1, widths = c(4.8, 3.8, 3.4)),
    patchwork::wrap_plots(panels[c("F", "G", "H")], nrow = 1, widths = c(3.9, 3.3, 4.8)),
    patchwork::wrap_plots(panels[c("I", "J")], nrow = 1, widths = c(1, 1),
                          guides = "collect") & theme(legend.position = "bottom"),
    patchwork::wrap_plots(panels[c("K", "L", "M")], nrow = 1, widths = c(1, 1, 1),
                          guides = "collect") & theme(legend.position = "bottom"))
  # Keep the first row independent of the right-side pseudotime guide below.
  rows[[1]] <- patchwork::free(rows[[1]], type = "panel", side = "r")
  figure <- patchwork::wrap_plots(rows, ncol = 1, heights = c(1.12, 1.08, 1.10, 1.04, 0.90))
  write_source(plot_data, "Figure4A_G_CytoTRACE_source.csv")
  write_source(trajectory_data, "Figure4H_M_trajectory_ESR2_source.csv")
  write_source(trajectory$roots, "Figure4_trajectory_roots.csv")
  if (!is.null(trajectory$graph_source)) {
    write_source(trajectory$graph_source, "Figure4_graph_provenance.csv")
  }
  write_source(manifest, "Figure4_input_manifest.csv")
  write_source(patients, "Figure4_patient_metadata.csv")
  coverage <- aggregate(finite_pseudotime ~ gender + group.a, trajectory_data,
                        function(x) c(total = length(x), finite = sum(x)))
  write_source(coverage, "Figure4_pseudotime_coverage.csv")
  curve_source <- do.call(rbind, lapply(c("K", "L", "M"), function(tag) {
    z <- ggplot2::ggplot_build(panels[[tag]])$data[[1]]
    data.frame(panel = tag, response = names(curve_colors)[match(z$colour, unname(curve_colors))],
               pseudotime = z$x, ESR2_LOESS = z$y)
  }))
  write_source(curve_source, "Figure4K_M_LOESS_curves.csv")
  # Document the cell counts and unmodified box statistics for C/E/G.
  distribution_summary <- do.call(rbind, lapply(c(C = "anno3", E = "group.a", G = "gender"), function(column) {
    tag <- switch(column, anno3 = "C", group.a = "E", gender = "G")
    group_order <- names(sort(tapply(plot_data$CytoTRACE, plot_data[[column]],
                                     stats::median), decreasing = TRUE))
    do.call(rbind, lapply(group_order, function(group) {
      values <- plot_data$CytoTRACE[plot_data[[column]] == group]
      q <- stats::quantile(values, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
      data.frame(panel = tag, group = group, n_cells = length(values),
                 minimum = q[1], q25 = q[2], median = q[3], q75 = q[4], maximum = q[5])
    }))
  }))
  write_source(distribution_summary[distribution_summary$panel %in% c("E", "G"), , drop = FALSE],
               "Figure4E_G_distribution_summary.csv")
  write_source(distribution_summary[distribution_summary$panel == "C", , drop = FALSE],
               "Figure4C_distribution_summary.csv")
  message("Exporting the integrated Figure 4 and all 13 panels...")
  figure4_export(figure, file.path(output_dir, "Figure4_integrated"),
                 width_mm, height_mm, formats, dpi)
  panel_sizes <- list(A = c(90, 70), B = c(90, 70), C = c(112, 68),
                      D = c(73, 70), E = c(60, 62), F = c(73, 70), G = c(60, 62), H = c(83, 70),
                      I = c(105, 75), J = c(105, 75), K = c(65, 57), L = c(65, 57), M = c(65, 57))
  for (tag in names(panels)) {
    dims <- panel_sizes[[tag]]
    figure4_export(panels[[tag]], file.path(panel_dir, paste0("Figure4", tag)),
                   dims[1], dims[2], formats, dpi)
  }
  methods <- c(
    "Figure 4 | CytoTRACE, inferred trajectories and ESR2 expression in malignant plasma cells.",
    "A/B/C: subcluster phenotypes, saved CytoTRACE scores, and score distributions across nine subclusters.",
    "D/E: OR/SOR phenotypes and CytoTRACE distributions. F/G: sex phenotypes and distributions.",
    "C/E/G: all-cell violin distributions (bounded to 0-1, each scaled to equal maximum width) with white boxplots showing the median and interquartile range; whiskers follow the default 1.5-IQR rule. No cell subsampling is used.",
    "CytoTRACE scores and Component1/2 coordinates are reused from the historical exported tables.",
    "H: Monocle 3 graph on the saved Seurat UMAP. I: same graph faceted by sex. J: same graph faceted by response.",
    paste("Root selection:", paste(unique(trajectory$roots$method), collapse = "; ")),
    "K/L/M: ESR2 raw RNA counts versus finite pseudotime in all/female/male cells, coloured by OR/SOR; LOESS span 0.75, no confidence ribbon.",
    "Negative LOESS predictions are not clipped, matching the historical smoothing method; raw counts remain nonnegative.",
    "The original interactive roots/CDS were not saved. H-M are reconstructed and cannot be claimed as an exact numerical reproduction of the submitted PDF.",
    if (!is.null(root_positions)) "Manual reference run: roots are approximate positions read from H before comparing K-M; no curve-based root optimization was performed." else "Roots were selected using the configured mode.",
    "No cell-level P values or patient-level inference are added. Pseudotime is an ordering of cross-sectional states, not elapsed time.",
    sprintf("Input: %d cells from %d patients; %d cells have finite pseudotime.",
            nrow(meta), nrow(patients), nrow(curve_data)),
    "The composite layout is newly assembled; submitted Figure 4.pdf and all original data are preserved."
  )
  writeLines(methods, file.path(output_dir, "Figure4_legend_and_methods.txt"))
  original_font_pt <- c(base = 6.8, axis_text = 6, axis_title = 7,
                        legend_text = 5.5, legend_title = 6, strip_text = 7, plot_title = 8,
                        panel_tag = 12, root_number = 2 * ggplot2::.pt)
  utils::write.csv(data.frame(element = names(original_font_pt),
                              original_pt = unname(original_font_pt), exported_pt = unname(font_pt(original_font_pt))),
                   file.path(output_dir, "Figure4_typography.csv"), row.names = FALSE)
  capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  message("Completed Figure 4 A-M: ", output_dir)
  invisible(list(figure = figure, panels = panels, cds = cds,
                 cytotrace_data = plot_data, trajectory_data = trajectory_data,
                 roots = trajectory$roots, input_manifest = manifest))
}

# ======================== 5. Run and expose editable plots ========================
if (!identical(Sys.getenv("FIGURE4_RUN"), "0")) {
  figure4_result <- build_figure4(
    analysis_file = figure4_analysis_file,
    cyto_anno_file = figure4_cyto_anno_file,
    cyto_response_file = figure4_cyto_response_file,
    cyto_gender_file = figure4_cyto_gender_file,
    output_dir = figure4_output_dir, formats = figure4_formats,
    width_mm = figure4_width_mm, height_mm = figure4_height_mm, dpi = figure4_dpi,
    seed = figure4_seed, cores = figure4_cores,
    root_mode = figure4_root_mode, root_cells = figure4_root_cells,
    root_pr_nodes = figure4_root_pr_nodes,
    root_positions = figure4_root_positions, graph_cache_file = figure4_graph_cache_file,
    cache_file = figure4_cache_file, rebuild_trajectory = figure4_rebuild_trajectory,
    font_increase_pt = figure4_font_increase_pt
  )
  Figure4 <- figure4_result$figure
  p4A <- figure4_result$panels$A
  p4B <- figure4_result$panels$B
  p4C <- figure4_result$panels$C
  p4D <- figure4_result$panels$D
  p4E <- figure4_result$panels$E
  p4F <- figure4_result$panels$F
  p4G <- figure4_result$panels$G
  p4H <- figure4_result$panels$H
  p4I <- figure4_result$panels$I
  p4J <- figure4_result$panels$J
  p4K <- figure4_result$panels$K
  p4L <- figure4_result$panels$L
  p4M <- figure4_result$panels$M
}
