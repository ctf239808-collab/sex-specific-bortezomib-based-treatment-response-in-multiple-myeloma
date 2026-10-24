 #Figure 6
# ======================== 1. Parameters / parameter ========================
figure6_config <- list(
  female_file = "/Users/tf/e/ndmm/HM38/FMM.hdWGCNA/female.hdWGCNA_object.rds",
  male_file = "/Users/tf/e/ndmm/HM38/Male.MM.hdWGCNA/Male.MM.hdWGCNA_object.rds",
  female_umap_file = "/Users/tf/e/ndmm/HM38/FMM.hdWGCNA/FMM.umap_df.csv",
  male_umap_file = "/Users/tf/e/ndmm/HM38/Male.MM.hdWGCNA/Male.MM.umap_df.csv",
  female_tom_file = "/Users/tf/e/ndmm/HM38/FMM.hdWGCNA/female.MM_TOM.rda",
  male_tom_file = "/Users/tf/e/ndmm/HM38/TOM/male.MM_TOM.rda",
  output_dir = "/Users/tf/e/ndmm/ESR2/submit/NC/Figure6_code/outputs",
  pathway_cache_dir = "/Users/tf/e/ndmm/ESR2/submit/NC/Figure6_code/cache",
  # Optional explicit local database tables; NULL uses a cache or the database API.
  # MSigDB: data.frame with gs_name, gene_symbol and ncbi_gene/entrez_gene.
  # KEGG: see f6_read_live_kegg() for the serialized download_KEGG cache schema.
  msigdb_legacy_file = NULL,
  kegg_cache_file = NULL,
  allow_pathway_download = TRUE,
  patient_col = "orig.ident",
  male_module_j = "Male.MM1",
  formats = c("pdf", "png"),
  width_mm = 230,
  height_mm = 270,
  dpi = 300,
  seed = 12345L,
  network_hubs_inner = 10L,
  network_hubs_outer = 15L,
  network_connections = 500L,
  topology_edge_fraction = 0.1,
  topology_label_hubs = 5L,
  # Module order copied from the submitted C/E and D/F panels.
  female_dot_order = c("FMM8", "FMM6", "FMM10", "FMM5", "FMM9", "FMM12",
                       "FMM11", "FMM4", "FMM7", "FMM3", "FMM2", "FMM1", "FMM13"),
  male_dot_order = c("Male.MM1", "Male.MM3", "Male.MM2")
)

# ======================== 2. Common helpers and validation ========================
figure6_check_packages <- function() {
  required <- c("Seurat", "SeuratObject", "Matrix", "hdWGCNA", "ggplot2",
                "patchwork", "ggrepel", "scales", "clusterProfiler", "enrichplot",
                "GSVA", "org.Hs.eg.db", "AnnotationDbi", "BiocParallel", "gson", "aplot")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". On this machine use /Library/Frameworks/R.framework/Resources/bin/Rscript.")
  if (!capabilities("cairo")) stop("Cairo graphics support is required for export.")
}

figure6_theme <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black", size = base_size - 1),
      axis.line = ggplot2::element_line(linewidth = 0.25),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1.5),
      legend.key.size = grid::unit(2.5, "mm"),
      plot.title = ggplot2::element_text(size = base_size + 1, hjust = 0.5),
      plot.margin = ggplot2::margin(3, 3, 3, 3))
}

figure6_write_table <- function(x, filename, config) {
  utils::write.csv(x, file.path(config$output_dir, "source_data", filename),
                   row.names = FALSE, na = "")
}

figure6_export <- function(plot, stem, width_mm, height_mm, config) {
  for (format in config$formats) {
    filename <- paste0(stem, ".", format)
    write_one <- function() {
      switch(format,
             pdf = grDevices::cairo_pdf(filename, width = width_mm / 25.4,
                                        height = height_mm / 25.4, family = "sans", bg = "white"),
             png = grDevices::png(filename, width = width_mm / 25.4,
                                  height = height_mm / 25.4, units = "in", res = config$dpi,
                                  type = "cairo", bg = "white"),
             svg = grDevices::svg(filename, width = width_mm / 25.4,
                                  height = height_mm / 25.4, family = "sans", bg = "white"),
             tiff = grDevices::tiff(filename, width = width_mm / 25.4,
                                    height = height_mm / 25.4, units = "in", res = config$dpi,
                                    type = "cairo", compression = "lzw", bg = "white"))
      device <- grDevices::dev.cur()
      on.exit({
        if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
      }, add = TRUE)
      print(plot)
    }
    tryCatch(write_one(), error = function(e) stop("Export failed [", filename,
                                                   "]: ", conditionMessage(e), call. = FALSE))
  }
}

figure6_rna_data <- function(obj) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    # Legacy Assay treats `search` as exact names, while Assay5 accepts regex.
    layers <- SeuratObject::Layers(obj[["RNA"]], search = NA)
    if (!"data" %in% layers || any(grepl("^data\\.", layers))) {
      stop("RNA requires a single joined normalized data layer.")
    }
    SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
  } else {
    SeuratObject::GetAssayData(obj, assay = "RNA", slot = "data")
  }
}

figure6_prepare_object <- function(path, sex, dot_order, config) {
  message("Reading ", sex, " historical hdWGCNA object...")
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("Not a Seurat object: ", path)
  meta <- obj[[]]
  required <- c("gender", "group.a", "orig.ident", "ESR2_group")
  if (!all(required %in% names(meta))) stop(sex, " metadata require: ",
                                            paste(setdiff(required, names(meta)), collapse = ", "))
  if (anyNA(meta[, required]) || any(as.matrix(meta[, required]) == "")) {
    stop(sex, " metadata contain missing labels.")
  }
  if (!all(as.character(meta$gender) == sex)) stop("Sex mismatch in ", path)
  if (!all(as.character(meta$group.a) %in% c("OR", "SOR"))) {
    stop("Unexpected response labels in ", path)
  }
  active <- obj@misc$active_wgcna
  if (length(active) != 1L || is.null(obj@misc[[active]]$hMEs)) {
    stop(sex, " object has no saved harmonized module eigengenes.")
  }
  modules <- hdWGCNA::GetModules(obj, wgcna_name = active)
  if (!all(c("gene_name", "module", "color") %in% names(modules)) ||
      anyDuplicated(modules$gene_name)) stop("Invalid module table for ", sex)
  hmes <- as.data.frame(hdWGCNA::GetMEs(obj, harmonized = TRUE, wgcna_name = active))
  cells <- colnames(obj)
  if (is.null(rownames(hmes)) || !all(cells %in% rownames(hmes))) {
    stop("Saved eigengene cell names do not match ", sex, " metadata.")
  }
  hmes <- hmes[cells, , drop = FALSE]
  if (!all(dot_order %in% names(hmes))) stop("Expected modules missing for ", sex)
  if (!setequal(dot_order, setdiff(as.character(unique(modules$module)), "grey"))) {
    stop("Module identities differ from the submitted Figure 6 for ", sex)
  }
  if (any(!is.finite(as.matrix(hmes[, dot_order, drop = FALSE])))) {
    stop("Nonfinite harmonized module eigengenes in ", sex)
  }
  # Preserve the saved historical labels exactly; do not compute a new median split.
  labels <- as.character(meta$ESR2_group)
  labels[labels %in% c("Low", "low", "ESR2_low")] <- "ESR2_Low"
  labels[labels %in% c("High", "high", "ESR2_high")] <- "ESR2_High"
  if (!setequal(unique(labels), c("ESR2_Low", "ESR2_High"))) {
    stop("Both saved ESR2_Low/ESR2_High labels are required for ", sex)
  }
  obj$ESR2_group <- factor(labels, levels = c("ESR2_High", "ESR2_Low"))
  obj$group.a <- factor(as.character(meta$group.a), levels = c("OR", "SOR"))
  # Replace any stale metadata columns with aligned saved hMEs, without cbind duplicates.
  obj <- SeuratObject::AddMetaData(obj, metadata = hmes[, dot_order, drop = FALSE])
  SeuratObject::DefaultAssay(obj) <- "RNA"
  rna <- figure6_rna_data(obj)
  if (!"ESR2" %in% rownames(rna)) stop("ESR2 is missing from RNA/data in ", sex)
  esr2 <- as.numeric(rna["ESR2", cells])
  # Record how the historical group relates to expression, rather than silently relabelling.
  group_audit <- data.frame(cell = cells, patient = as.character(meta$orig.ident),
                            sex = sex, response = as.character(meta$group.a), ESR2_log_normalized = esr2,
                            ESR2_group = labels, agrees_with_RNA_positive =
                              labels == ifelse(esr2 > 0, "ESR2_High", "ESR2_Low"))
  figure6_write_table(group_audit, paste0(sex, "_cell_metadata.csv"), config)
  figure6_write_table(data.frame(gene = modules$gene_name, modules[,
                                                                   setdiff(names(modules), "gene_name"), drop = FALSE]),
                      paste0(sex, "_module_assignments_kME.csv"), config)
  figure6_write_table(data.frame(cell = cells, hmes[, dot_order, drop = FALSE]),
                      paste0(sex, "_hMEs.csv"), config)
  obj
}

# ======================== 3. A-F: Cell UMAP and module DotPlots ========================
figure6_feature_grid <- function(obj, module_names, ncol) {
  args <- list(seurat_obj = obj, module_names = module_names, features = "hMEs",
               reduction = "umap", point_size = 0.12, raster = FALSE)
  # hdWGCNA renamed the historical `order` argument to `order_points`.
  api <- names(formals(hdWGCNA::ModuleFeaturePlot))
  if ("order_points" %in% api) args$order_points <- TRUE else args$order <- TRUE
  plots <- do.call(hdWGCNA::ModuleFeaturePlot, args)
  if (inherits(plots, "ggplot")) plots <- list(plots)
  plots <- lapply(plots, function(p) p +
                    ggplot2::theme(plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
                                   plot.margin = ggplot2::margin(1, 1, 1, 1), legend.text = ggplot2::element_text(size = 5),
                                   legend.key.height = grid::unit(1.7, "mm")))
  patchwork::wrap_plots(plots, ncol = ncol)
}

figure6_dotplot <- function(obj, modules, grouping, panel, config) {
  p <- Seurat::DotPlot(obj, features = modules, assay = "RNA", group.by = grouping,
                       dot.scale = 3.5, scale = TRUE) + Seurat::RotatedAxis() +
    ggplot2::scale_colour_gradient2(low = "blue", mid = "grey95", high = "red") +
    figure6_theme(7) + ggplot2::labs(x = "Features", y = "Identity") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5.5),
                   legend.position = "right")
  # Seurat DotPlot retains expm1/scale handling for metadata features, as in the
  # original scripts. This is an exploratory cell-level display, not a patient test.
  figure6_write_table(p$data, paste0("panel_", panel, "_DotPlot.csv"), config)
  p
}

# ======================== 4. G-I: TOM networks and saved gene UMAP ========================
figure6_read_tom <- function(path, modules) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  candidates <- intersect(c("consTomDS", "TOM", "tom"), loaded)
  if (!length(candidates)) stop("No consTomDS/TOM matrix found in ", path)
  tom <- as.matrix(env[[candidates[1]]])
  genes <- as.character(modules$gene_name)
  if (!identical(dim(tom), c(length(genes), length(genes)))) {
    stop("TOM dimensions do not match saved module gene order: ", path)
  }
  generic_names <- function(x) is.null(x) || identical(x, as.character(seq_along(genes)))
  if (generic_names(rownames(tom)) && generic_names(colnames(tom))) {
    # Historical consTomDS is an unnamed dist; as.matrix generates "1", "2", ... .
    # hdWGCNA::GetTOM likewise assigns the saved module-table gene order here.
    attr(tom, "figure6_gene_order") <- "positional saved module-table order"
  } else if (setequal(rownames(tom), genes) && setequal(colnames(tom), genes) &&
             !anyDuplicated(rownames(tom)) && !anyDuplicated(colnames(tom))) {
    tom <- tom[genes, genes, drop = FALSE]
  } else stop("TOM gene labels do not match the saved module table: ", path)
  dimnames(tom) <- list(genes, genes)
  tom
}

figure6_module_genes <- function(modules, module, n = Inf) {
  score <- paste0("kME_", module)
  if (!score %in% names(modules)) stop("Missing kME column: ", score)
  x <- modules[as.character(modules$module) == module, , drop = FALSE]
  x <- x[is.finite(x[[score]]), , drop = FALSE]
  x <- x[order(-x[[score]], x$gene_name), , drop = FALSE]
  if (!nrow(x)) stop("No finite hub scores for ", module)
  utils::head(x, min(n, nrow(x)))
}

figure6_circle_network <- function(modules, tom, module, config) {
  nodes <- figure6_module_genes(modules, module,
                                config$network_hubs_inner + config$network_hubs_outer)
  n_inner <- min(config$network_hubs_inner, nrow(nodes) - 1L)
  circle <- function(n, radius) {
    angle <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
    cbind(x = radius * cos(angle), y = radius * sin(angle))
  }
  xy <- rbind(circle(n_inner, 0.48), circle(nrow(nodes) - n_inner, 1))
  nodes$x <- xy[, 1]; nodes$y <- xy[, 2]
  a <- tom[nodes$gene_name, nodes$gene_name, drop = FALSE]
  idx <- which(upper.tri(a) & is.finite(a) & a > 0, arr.ind = TRUE)
  edges <- data.frame(from = nodes$gene_name[idx[, 1]], to = nodes$gene_name[idx[, 2]],
                      weight = a[idx], x = nodes$x[idx[, 1]], y = nodes$y[idx[, 1]],
                      xend = nodes$x[idx[, 2]], yend = nodes$y[idx[, 2]])
  # The historical n_conns=500 counted directed matrix entries. The explicit
  # reconstruction keeps up to 500 unique, non-self undirected connections.
  edges <- edges[order(-edges$weight, edges$from, edges$to), , drop = FALSE]
  edges <- utils::head(edges, config$network_connections)
  color <- as.character(nodes$color[1])
  figure6_write_table(nodes, paste0("panel_G_", module, "_nodes.csv"), config)
  figure6_write_table(edges, paste0("panel_G_", module, "_edges.csv"), config)
  # aes() captures this function's environment; unbind the full TOM first so a
  # returned ggplot cannot retain a multi-GB matrix after the network helper ends.
  rm(tom, a, idx)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x, y, xend = xend, yend = yend,
                                                     alpha = weight), colour = color, linewidth = 0.16) +
    ggplot2::scale_alpha_continuous(range = c(0.03, 0.35), guide = "none") +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x, y), shape = 21,
                        fill = color, colour = "grey35", stroke = 0.2, size = 1.2) +
    ggrepel::geom_text_repel(data = nodes, ggplot2::aes(x, y, label = gene_name),
                             size = 1.7, fontface = "italic", seed = config$seed, max.overlaps = Inf,
                             min.segment.length = 0, segment.size = 0.15, box.padding = 0.1) +
    ggplot2::coord_equal(clip = "off", xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3)) +
    ggplot2::theme_void() + ggplot2::labs(title = module) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
                   plot.margin = ggplot2::margin(3, 3, 3, 3))
}

figure6_topology <- function(modules, tom, csv, panel, config) {
  nodes <- utils::read.csv(csv, row.names = 1, check.names = FALSE,
                           stringsAsFactors = FALSE)
  required <- c("gene", "module", "color", "hub", "kME", "UMAP1", "UMAP2")
  if (!all(required %in% names(nodes)) || anyDuplicated(nodes$gene) ||
      !all(nodes$gene %in% rownames(tom)) ||
      any(!is.finite(as.matrix(nodes[, c("UMAP1", "UMAP2", "kME")])))) {
    stop("Invalid saved module UMAP table: ", csv)
  }
  match_module <- match(nodes$gene, modules$gene_name)
  if (!all(as.character(nodes$module) == as.character(modules$module[match_module])) ||
      !all(as.character(nodes$color) == as.character(modules$color[match_module]))) {
    stop("Saved gene UMAP assignments do not match this network: ", csv)
  }
  mods <- setdiff(unique(as.character(nodes$module)), "grey")
  edges <- lapply(mods, function(mod) {
    all_genes <- nodes$gene[nodes$module == mod]
    hubs <- nodes$gene[nodes$module == mod & nodes$hub == "hub"]
    if (!length(hubs)) stop("No saved UMAP hubs for module ", mod)
    a <- tom[all_genes, hubs, drop = FALSE]
    idx <- which(is.finite(a) & a > 0, arr.ind = TRUE)
    x <- data.frame(from = all_genes[idx[, 1]], to = hubs[idx[, 2]], weight = a[idx])
    x <- x[x$from != x$to, , drop = FALSE]
    # Deduplicate undirected edges, then sample 10% within each module, with seed.
    pair <- paste(pmin(x$from, x$to), pmax(x$from, x$to), sep = "\r")
    x <- x[!duplicated(pair), , drop = FALSE]
    take <- round(nrow(x) * config$topology_edge_fraction)
    if (!take) return(NULL)
    x <- x[sample.int(nrow(x), take), , drop = FALSE]
    a_idx <- match(x$from, nodes$gene); b_idx <- match(x$to, nodes$gene)
    x$x <- nodes$UMAP1[a_idx]; x$y <- nodes$UMAP2[a_idx]
    x$xend <- nodes$UMAP1[b_idx]; x$yend <- nodes$UMAP2[b_idx]
    x$color <- as.character(nodes$color[a_idx])
    x
  })
  edges <- do.call(rbind, edges)
  if (is.null(edges) || !nrow(edges)) stop("No drawable topology edges in ", panel)
  labels <- unlist(lapply(mods, function(mod) figure6_module_genes(modules, mod,
                                                                   config$topology_label_hubs)$gene_name), use.names = FALSE)
  nodes$label <- ifelse(nodes$gene %in% labels, nodes$gene, "")
  centroids <- stats::aggregate(nodes[, c("UMAP1", "UMAP2")],
                                list(module = as.character(nodes$module)), stats::median)
  centroids <- centroids[centroids$module != "grey", , drop = FALSE]
  centroids$color <- nodes$color[match(centroids$module, nodes$module)]
  figure6_write_table(nodes, paste0("panel_", panel, "_topology_nodes.csv"), config)
  figure6_write_table(edges, paste0("panel_", panel, "_topology_edges.csv"), config)
  rm(tom)
  # Draw the saved coordinates; no stochastic UMAP recomputation occurs.
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x, y, xend = xend, yend = yend,
                                                     colour = color, alpha = weight), linewidth = 0.12) +
    ggplot2::scale_alpha_continuous(range = c(0.01, 0.16), guide = "none") +
    ggplot2::geom_point(data = nodes, ggplot2::aes(UMAP1, UMAP2, colour = color,
                                                   size = pmax(kME, 0)), alpha = 0.9) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_size_continuous(range = c(0.1, 0.65), guide = "none") +
    ggrepel::geom_text_repel(data = nodes[nodes$label != "", , drop = FALSE],
                             ggplot2::aes(UMAP1, UMAP2, label = label), size = 1.8, seed = config$seed,
                             max.overlaps = Inf, segment.size = 0.12, box.padding = 0.15) +
    ggrepel::geom_text_repel(data = centroids,
                             ggplot2::aes(UMAP1, UMAP2, label = module, colour = color),
                             fontface = "bold", size = 2.2, seed = config$seed,
                             max.overlaps = Inf, min.segment.length = Inf) +
    ggplot2::coord_equal() + figure6_theme() +
    ggplot2::labs(x = "TOM UMAP-1", y = "TOM UMAP-2") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
}

figure6_network_panels <- function(obj, tom_file, umap_file, sex, config) {
  message("Building ", sex, " TOM network panels...")
  modules <- hdWGCNA::GetModules(obj)
  tom <- figure6_read_tom(tom_file, modules)
  selected <- if (sex == "Female") c("FMM7", "FMM8") else "Male.MM1"
  circles <- stats::setNames(lapply(selected, function(mod)
    figure6_circle_network(modules, tom, mod, config)), selected)
  topology <- figure6_topology(modules, tom, umap_file,
                               if (sex == "Female") "H" else "I", config)
  list(circles = circles, topology = topology)
}

# ======================== 5. J-L: Enrichment, GSVA and GSEA ========================
# All pathway helpers are included below; this file requires no companion script.

# -----------------------------------------------------------------------------
# Figure 6 J-L: historical KEGG enrichment, GSVA correlation and Proteasome GSEA.
# Source: Fig.5-new.R, lines 191-275, 859-1065, and 1067-1482.
# Pathway databases are frozen to local caches when available. New downloads must
# be enabled explicitly with config$allow_pathway_download = TRUE.
# -----------------------------------------------------------------------------

f6_pathway_option <- function(config, name, default = NULL) {
  value <- config[[name]]
  if (is.null(value)) default else value
}

f6_pathway_write_csv <- function(x, filename, config, row.names = FALSE) {
  utils::write.csv(x, file.path(config$output_dir, "source_data", filename),
                   row.names = row.names, na = "NA")
}

f6_read_legacy_kegg <- function(config) {
  cached <- file.path(config$pathway_cache_dir, "msigdb_legacy_human.rds")
  explicit <- f6_pathway_option(config, "msigdb_legacy_file")
  if (!is.null(explicit) && !file.exists(explicit))
    stop("Configured MSigDB file is missing: ", explicit, call. = FALSE)
  candidates <- unique(c(explicit, cached,
                         path.expand("~/Library/Caches/org.R-project.R/R/msigdbr/msigdb.2025.1.Hs.rds")))
  candidates <- candidates[!is.na(candidates) & file.exists(candidates)]
  if (length(candidates)) {
    input <- candidates[[1L]]
    raw <- readRDS(input)
    provenance <- list(source = normalizePath(input), md5 = unname(tools::md5sum(input)))
    if (is.list(raw) && !is.data.frame(raw) && !is.null(raw$table)) {
      provenance <- raw$provenance
      raw <- raw$table
    }
  } else {
    if (!isTRUE(f6_pathway_option(config, "allow_pathway_download", FALSE)))
      stop("J/K require human MSigDB C2:CP:KEGG_LEGACY. Set config$msigdb_legacy_file ",
           "to a local RDS, or enable config$allow_pathway_download.", call. = FALSE)
    if (!requireNamespace("msigdbr", quietly = TRUE))
      stop("Downloading MSigDB requires the msigdbr package.", call. = FALSE)
    args <- list(species = "Homo sapiens")
    if ("collection" %in% names(formals(msigdbr::msigdbr))) {
      args$collection <- "C2"
      args$subcollection <- "CP:KEGG_LEGACY"
    } else {
      args$category <- "C2"
      args$subcategory <- "CP:KEGG_LEGACY"
    }
    raw <- do.call(msigdbr::msigdbr, args)
    provenance <- list(source = "msigdbr human C2:CP:KEGG_LEGACY",
                       accessed_date = as.character(Sys.Date()))
  }
  raw <- as.data.frame(raw)
  collection_col <- intersect(c("gs_collection", "gs_cat"), names(raw))
  subcollection_col <- intersect(c("gs_subcollection", "gs_subcat"), names(raw))
  if (length(collection_col)) raw <- raw[raw[[collection_col[[1L]]]] == "C2", , drop = FALSE]
  if (length(subcollection_col))
    raw <- raw[raw[[subcollection_col[[1L]]]] == "CP:KEGG_LEGACY", , drop = FALSE]
  symbol_col <- intersect(c("gene_symbol", "db_gene_symbol"), names(raw))
  entrez_col <- intersect(c("ncbi_gene", "entrez_gene", "db_ncbi_gene"), names(raw))
  if (!nrow(raw) || !"gs_name" %in% names(raw) ||
      !length(symbol_col) || !length(entrez_col))
    stop("MSigDB input does not contain human KEGG_LEGACY gene sets with SYMBOL and ENTREZ IDs.",
         call. = FALSE)
  tbl <- unique(data.frame(gs_name = as.character(raw$gs_name),
                           gene_symbol = as.character(raw[[symbol_col[[1L]]]]),
                           entrez_gene = as.character(raw[[entrez_col[[1L]]]]),
                           stringsAsFactors = FALSE))
  tbl <- tbl[complete.cases(tbl) & nzchar(tbl$gene_symbol) & nzchar(tbl$entrez_gene), ]
  provenance$db_version <- if ("db_version" %in% names(raw))
    paste(unique(raw$db_version), collapse = ";") else if (!is.null(provenance$db_version))
      provenance$db_version else "not recorded in supplied cache"
  provenance$collection <- "C2:CP:KEGG_LEGACY"
  if (!file.exists(cached)) saveRDS(list(table = tbl, provenance = provenance), cached)
  list(table = tbl, provenance = provenance)
}

f6_read_live_kegg <- function(config) {
  default_file <- file.path(config$pathway_cache_dir, "kegg_hsa_download_KEGG.rds")
  input <- f6_pathway_option(config, "kegg_cache_file", default_file)
  if (!is.null(config$kegg_cache_file) && !file.exists(input)) {
    stop("Configured KEGG cache file is missing: ", input, call. = FALSE)
  }
  if (file.exists(input)) {
    raw <- readRDS(input)
    if (is.list(raw) && !is.null(raw$data)) {
      provenance <- raw$provenance
      raw <- raw$data
    } else provenance <- list(source = normalizePath(input),
                              md5 = unname(tools::md5sum(input)),
                              accessed_date = "not recorded in supplied cache")
  } else {
    if (!isTRUE(f6_pathway_option(config, "allow_pathway_download", FALSE)))
      stop("K (FoxO) / L require a download_KEGG('hsa') cache: ", input,
           ". Supply config$kegg_cache_file or enable config$allow_pathway_download. ",
           "KEGG_LEGACY is not substituted for these panels.", call. = FALSE)
    raw <- clusterProfiler::download_KEGG("hsa")
    provenance <- list(source = "clusterProfiler::download_KEGG('hsa')",
                       accessed_date = as.character(Sys.Date()),
                       clusterProfiler_version = as.character(utils::packageVersion("clusterProfiler")))
    saveRDS(list(data = raw, provenance = provenance), default_file)
  }
  needed <- c("KEGGPATHID2EXTID", "KEGGPATHID2NAME")
  if (!all(needed %in% names(raw)))
    stop("KEGG cache must contain KEGGPATHID2EXTID and KEGGPATHID2NAME from download_KEGG.",
         call. = FALSE)
  for (field in needed) {
    raw[[field]] <- as.data.frame(raw[[field]], stringsAsFactors = FALSE)
    if (ncol(raw[[field]]) != 2L)
      stop("Invalid KEGG cache component: ", field, call. = FALSE)
    names(raw[[field]]) <- c("from", "to")
    raw[[field]]$from <- sub("^path:", "", as.character(raw[[field]]$from))
    raw[[field]]$to <- as.character(raw[[field]]$to)
  }
  raw$KEGGPATHID2EXTID$to <- sub("^hsa:", "", raw$KEGGPATHID2EXTID$to)
  if (!all(c("hsa04068", "hsa03050") %in% raw$KEGGPATHID2EXTID$from))
    stop("KEGG cache lacks FoxO (hsa04068) or Proteasome (hsa03050).", call. = FALSE)
  list(data = raw, provenance = provenance)
}

f6_patient_pathway_input <- function(obj, sex, config) {
  patient_col <- f6_pathway_option(config, "patient_col", "orig.ident")
  meta <- obj[[]]
  if (!patient_col %in% names(meta))
    stop(sex, " object lacks patient metadata column: ", patient_col, call. = FALSE)
  patients <- as.character(meta[[patient_col]])
  if (anyNA(patients) || any(!nzchar(patients)))
    stop(sex, " object contains missing patient IDs.", call. = FALSE)
  patient_ids <- unique(patients)
  if (length(patient_ids) < 3L)
    stop(sex, " pathway correlation requires at least 3 patients.", call. = FALSE)
  log_data <- tryCatch(
    SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) stop(sex, " RNA data layer could not be read. ",
                             "For split Seurat v5 layers, first join RNA data layers. ", conditionMessage(e), call. = FALSE))
  if (!"ESR2" %in% rownames(log_data))
    stop(sex, " RNA assay has no ESR2 feature.", call. = FALSE)
  # Safe temporary patient keys prevent Seurat name sanitization (e.g. g-prefix).
  keys <- sprintf("F6P%05d", seq_along(patient_ids))
  key_for_cell <- keys[match(patients, patient_ids)]
  obj[["Figure6_patient_group"]] <- key_for_cell
  avg_args <- list(object = obj, group.by = "Figure6_patient_group", assays = "RNA",
                   verbose = FALSE)
  if ("layer" %in% names(formals(Seurat::AverageExpression))) avg_args$layer <- "data"
  else avg_args$slot <- "data"
  # Historical AverageExpression(RNA, data) averages expm1(log-normalized data).
  average_expression <- as.matrix(do.call(Seurat::AverageExpression, avg_args)$RNA)
  if (!all(keys %in% colnames(average_expression)))
    stop("Could not align ", sex, " patient means returned by AverageExpression.", call. = FALSE)
  average_expression <- average_expression[, keys, drop = FALSE]
  colnames(average_expression) <- patient_ids
  log_esr2 <- as.numeric(log_data["ESR2", rownames(meta)])
  mean_log_esr2 <- vapply(patient_ids, function(p)
    mean(log_esr2[patients == p], na.rm = TRUE), numeric(1))
  names(mean_log_esr2) <- patient_ids
  if (stats::sd(mean_log_esr2) == 0 || stats::sd(average_expression["ESR2", ]) == 0)
    stop(sex, " patient-level ESR2 has zero variance; correlation is undefined.", call. = FALSE)
  f6_pathway_write_csv(data.frame(patient = patient_ids,
                                  n_cells = as.integer(table(factor(patients, levels = patient_ids))),
                                  ESR2_mean_log_normalized = unname(mean_log_esr2),
                                  ESR2_AverageExpression = as.numeric(average_expression["ESR2", ])),
                       paste0("Figure6_KL_", sex, "_patient_ESR2.csv"), config)
  saveRDS(average_expression, file.path(config$output_dir, "source_data",
                                        paste0("Figure6_KL_", sex, "_AverageExpression_RNA_data.rds")))
  list(expr = average_expression, ESR2_mean_log = mean_log_esr2)
}

f6_symbol_mapping <- function(symbols) {
  map <- clusterProfiler::bitr(unique(symbols), fromType = "SYMBOL", toType = "ENTREZID",
                               OrgDb = get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db")), drop = TRUE)
  map <- map[!duplicated(map$SYMBOL), , drop = FALSE]
  if (!nrow(map)) stop("No gene symbols could be mapped to ENTREZ IDs.", call. = FALSE)
  map
}

f6_run_gsva <- function(expr, sets) {
  expr <- as.matrix(expr)
  if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
    param <- GSVA::gsvaParam(exprData = expr, geneSets = sets, kcdf = "Gaussian")
    GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
  } else {
    # Same historical GSVA method, using the pre-1.50 package API.
    GSVA::gsva(expr, sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
  }
}

f6_safe_pearson <- function(x, y, label) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0)
    stop("Undefined Pearson correlation for ", label, ".", call. = FALSE)
  unname(stats::cor(x[ok], y[ok], method = "pearson"))
}

f6_pathway_correlations <- function(input, sex, legacy, kegg, config) {
  required <- c("KEGG_PROTEASOME", "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_CELL_CYCLE",
                "KEGG_APOPTOSIS", "KEGG_DNA_REPLICATION", "KEGG_BASE_EXCISION_REPAIR",
                "KEGG_TGF_BETA_SIGNALING_PATHWAY", "KEGG_VEGF_SIGNALING_PATHWAY",
                "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS", "KEGG_CELL_ADHESION_MOLECULES_CAMS",
                "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION")
  # Keep the full historical gene-set collection during GSVA, then select rows.
  legacy_sets <- split(legacy$table$gene_symbol, legacy$table$gs_name)
  scores <- f6_run_gsva(input$expr, legacy_sets)
  if (!all(required %in% rownames(scores)))
    stop(sex, " GSVA lacks expected KEGG_LEGACY pathways: ",
         paste(setdiff(required, rownames(scores)), collapse = ", "), call. = FALSE)
  scores <- scores[, colnames(input$expr), drop = FALSE]
  correlation <- vapply(required, function(path)
    f6_safe_pearson(scores[path, ], input$ESR2_mean_log[colnames(scores)],
                    paste(sex, path)), numeric(1))
  # Historical exception: FoxO came from current KEGG/ENTREZ GSVA and used
  # AverageExpression ESR2, rather than the mean log-normalized ESR2 above.
  map <- f6_symbol_mapping(rownames(input$expr))
  expr_entrez <- input$expr[map$SYMBOL, , drop = FALSE]
  rownames(expr_entrez) <- map$ENTREZID
  # Remove duplicate ENTREZ rows deterministically (first SYMBOL mapping), needed
  # by GSVA's unique-feature validation; the historical code omitted this guard.
  expr_entrez <- expr_entrez[!duplicated(rownames(expr_entrez)), , drop = FALSE]
  kegg_sets <- split(kegg$data$KEGGPATHID2EXTID$to, kegg$data$KEGGPATHID2EXTID$from)
  foxo_scores <- f6_run_gsva(expr_entrez, kegg_sets)
  if (!"hsa04068" %in% rownames(foxo_scores))
    stop(sex, " KEGG GSVA did not return FoxO hsa04068.", call. = FALSE)
  foxo_scores <- foxo_scores[, colnames(input$expr), drop = FALSE]
  correlation <- c(correlation, KEGG_FOXO_SIGNALING_PATHWAY =
                     f6_safe_pearson(foxo_scores["hsa04068", ], input$expr["ESR2", ], paste(sex, "FoxO")))
  f6_pathway_write_csv(as.data.frame(scores), paste0("Figure6_K_", sex,
                                                     "_GSVA_KEGG_LEGACY.csv"), config, row.names = TRUE)
  f6_pathway_write_csv(as.data.frame(foxo_scores), paste0("Figure6_K_", sex,
                                                          "_GSVA_KEGG_current_Entrez.csv"), config, row.names = TRUE)
  f6_pathway_write_csv(map, paste0("Figure6_K_", sex, "_SYMBOL_to_ENTREZ.csv"), config)
  correlation
}

f6_proteasome_gsea <- function(input, sex, kegg, config) {
  # Pearson correlation of every gene with ESR2, using patient AverageExpression.
  correlation <- apply(input$expr, 1L, function(x)
    suppressWarnings(stats::cor(x, input$expr["ESR2", ], use = "pairwise.complete.obs",
                                method = "pearson")))
  correlation <- sort(correlation[is.finite(correlation)], decreasing = TRUE)
  map <- f6_symbol_mapping(names(correlation))
  ranks <- correlation[map$SYMBOL]
  names(ranks) <- map$ENTREZID
  keep <- !duplicated(names(ranks))
  rank_source <- data.frame(SYMBOL = map$SYMBOL[keep], ENTREZID = names(ranks)[keep],
                            ESR2_Pearson_r = unname(ranks[keep]))
  ranks <- sort(ranks[keep], decreasing = TRUE)
  rank_source <- rank_source[match(names(ranks), rank_source$ENTREZID), ]
  f6_pathway_write_csv(rank_source, paste0("Figure6_L_", sex, "_ranked_genes.csv"), config)
  # gseKEGG accepts a GSON object, allowing the exact cached download_KEGG gene
  # sets to be used without an unrecorded database refresh or method substitution.
  pathways <- kegg$data$KEGGPATHID2EXTID
  names(pathways) <- c("gsid", "gene")
  path_names <- kegg$data$KEGGPATHID2NAME
  names(path_names) <- c("gsid", "name")
  frozen_kegg <- gson::gson(gsid2gene = pathways, gsid2name = path_names,
                            species = "hsa", gsname = "KEGG", keytype = "kegg",
                            accessed_date = if (is.null(kegg$provenance$accessed_date)) "not recorded" else
                              as.character(kegg$provenance$accessed_date))
  set.seed(config$seed)
  result <- clusterProfiler::gseKEGG(geneList = ranks, organism = frozen_kegg,
                                     minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)
  if (is.null(result) || !nrow(as.data.frame(result)))
    stop(sex, " gseKEGG returned no pathways at the historical pvalueCutoff=0.05.", call. = FALSE)
  tbl <- as.data.frame(result)
  f6_pathway_write_csv(tbl, paste0("Figure6_L_", sex, "_gseKEGG.csv"), config)
  saveRDS(result, file.path(config$output_dir, "source_data",
                            paste0("Figure6_L_", sex, "_gseKEGG.rds")))
  if (!"hsa03050" %in% tbl$ID)
    stop(sex, " Proteasome hsa03050 was not returned at pvalueCutoff=0.05. ",
         "The script does not relax this cutoff silently.", call. = FALSE)
  plot <- enrichplot::gseaplot2(result, geneSetID = "hsa03050",
                                title = paste0("GSEA - Proteasome (", sex, ")"))
  if (inherits(plot, "gglist")) plot <- aplot::as.patchwork(plot)
  else if (inherits(plot, c("ggplot", "patchwork"))) plot <- plot
  else if (is.list(plot)) plot <- patchwork::wrap_plots(plot, ncol = 1)
  else stop("Unsupported gseaplot2 plot class: ", paste(class(plot), collapse = ", "),
            call. = FALSE)
  plot & ggplot2::theme(text = ggplot2::element_text(size = 6.5),
                        axis.text = ggplot2::element_text(size = 5.5),
                        axis.title = ggplot2::element_text(size = 6),
                        plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
                        plot.margin = ggplot2::margin(2, 2, 2, 2))
}

build_figure6_pathways <- function(female_obj, male_obj, config) {
  required_packages <- c("Seurat", "SeuratObject", "GSVA", "clusterProfiler", "org.Hs.eg.db",
                         "enrichplot", "ggplot2", "patchwork", "aplot", "gson", "hdWGCNA")
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Figure 6 J-L missing R packages: ", paste(missing, collapse = ", "),
                            call. = FALSE)
  if (is.null(config$seed)) config$seed <- 12345L
  if (is.null(config$pathway_cache_dir))
    config$pathway_cache_dir <- file.path(config$output_dir, "pathway_cache")
  dir.create(file.path(config$output_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(config$pathway_cache_dir, recursive = TRUE, showWarnings = FALSE)
  legacy <- f6_read_legacy_kegg(config)
  kegg <- f6_read_live_kegg(config)
  f6_pathway_write_csv(legacy$table, "Figure6_JK_KEGG_LEGACY_TERM2GENE.csv", config)
  f6_pathway_write_csv(kegg$data$KEGGPATHID2EXTID, "Figure6_KL_KEGG_hsa_TERM2GENE.csv", config)
  f6_pathway_write_csv(kegg$data$KEGGPATHID2NAME, "Figure6_KL_KEGG_hsa_TERM2NAME.csv", config)
  saveRDS(list(MSigDB_KEGG_LEGACY = legacy$provenance, KEGG_current = kegg$provenance),
          file.path(config$output_dir, "source_data", "Figure6_JKL_database_provenance.rds"))
  female <- f6_patient_pathway_input(female_obj, "Female", config)
  male <- f6_patient_pathway_input(male_obj, "Male", config)
  
  # J: Male.MM1 (historical turquoise module), default enricher universe.
  modules <- hdWGCNA::GetModules(male_obj)
  module_col <- intersect(c("module", "module_name"), names(modules))
  module_name <- f6_pathway_option(config, "male_module_j", "Male.MM1")
  if (!length(module_col) || !"gene_name" %in% names(modules))
    stop("Male hdWGCNA modules require module and gene_name columns.", call. = FALSE)
  module_genes <- unique(modules$gene_name[as.character(modules[[module_col[[1L]]]]) == module_name])
  if (!length(module_genes))
    stop("J module not found: ", module_name, ". Available: ",
         paste(unique(modules[[module_col[[1L]]]]), collapse = ", "), call. = FALSE)
  module_map <- f6_symbol_mapping(module_genes)
  f6_pathway_write_csv(module_map, "Figure6_J_Male_MM1_genes_ENTREZ.csv", config)
  enrichment <- clusterProfiler::enricher(gene = module_map$ENTREZID,
                                          TERM2GENE = legacy$table[, c("gs_name", "entrez_gene")], pvalueCutoff = 0.05)
  if (is.null(enrichment) || !nrow(as.data.frame(enrichment)))
    stop("J Male.MM1 has no KEGG_LEGACY enrichment at pvalueCutoff=0.05.", call. = FALSE)
  f6_pathway_write_csv(as.data.frame(enrichment), "Figure6_J_Male_MM1_KEGG_enrichment.csv", config)
  J <- graphics::barplot(enrichment, showCategory = 20, label_format = 60) +
    ggplot2::ggtitle("Male.MM1 KEGG pathways") +
    ggplot2::scale_fill_gradient(low = "#D1FAE5", high = "#0D9488") +
    figure6_theme(6) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4.7))
  f6_pathway_write_csv(J$data, "Figure6_J_plotted_pathways.csv", config)
  
  # K: original row order; no clustering.
  ordered_paths <- c("KEGG_PROTEASOME", "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_CELL_CYCLE",
                     "KEGG_APOPTOSIS", "KEGG_DNA_REPLICATION", "KEGG_BASE_EXCISION_REPAIR",
                     "KEGG_TGF_BETA_SIGNALING_PATHWAY", "KEGG_VEGF_SIGNALING_PATHWAY",
                     "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS", "KEGG_CELL_ADHESION_MOLECULES_CAMS",
                     "KEGG_FOXO_SIGNALING_PATHWAY", "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION")
  correlations <- cbind(Female = f6_pathway_correlations(female, "Female", legacy, kegg, config)[ordered_paths],
                        Male = f6_pathway_correlations(male, "Male", legacy, kegg, config)[ordered_paths])
  f6_pathway_write_csv(as.data.frame(correlations), "Figure6_K_ESR2_pathway_correlations.csv",
                       config, row.names = TRUE)
  heat_data <- data.frame(Pathway = rep(sub("^KEGG_", "", ordered_paths), 2L),
                          Sex = rep(colnames(correlations), each = length(ordered_paths)),
                          Pearson_r = as.vector(correlations))
  heat_data$Pathway <- factor(heat_data$Pathway, levels = rev(sub("^KEGG_", "", ordered_paths)))
  heat_data$Sex <- factor(heat_data$Sex, levels = c("Female", "Male"))
  # Replaces only the pheatmap drawing API so K participates in a patchwork layout.
  K <- ggplot2::ggplot(heat_data, ggplot2::aes(x = Sex, y = Pathway, fill = Pearson_r)) +
    ggplot2::geom_tile(colour = "grey80", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                                  midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = 6) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 4.5),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5),
                   legend.text = ggplot2::element_text(size = 4.5),
                   legend.title = ggplot2::element_text(size = 5),
                   legend.key.height = grid::unit(2.5, "mm"))
  L <- patchwork::wrap_plots(f6_proteasome_gsea(female, "Female", kegg, config),
                             f6_proteasome_gsea(male, "Male", kegg, config), ncol = 2)
  methods <- c(
    "J: Male.MM1 genes -> SYMBOL/ENTREZ (org.Hs.eg.db) -> enricher, human MSigDB C2:CP:KEGG_LEGACY, pvalueCutoff=0.05; historical default background; top 20 barplot.",
    "K: GSVA(kcdf=Gaussian) on Seurat AverageExpression(RNA,data), which averages inverse-log normalized RNA expression. 11 LEGACY pathways correlate with patient mean log-normalized ESR2.",
    "K FoxO exception retained from history: download_KEGG hsa04068 / ENTREZ GSVA correlated with ESR2 from AverageExpression, not mean log-normalized ESR2.",
    "K: Pearson correlation; same 12 historical rows, no row/column clustering. ggplot tile drawing replaces pheatmap drawing for assembly. Duplicate ENTREZ expression rows are retained once (first SYMBOL mapping) to satisfy current GSVA feature validation.",
    "L: Gene-wise Pearson correlation with ESR2 using patient AverageExpression(RNA,data), nonfinite values removed, SYMBOL and ENTREZ duplicates handled in historical order, decreasing ranks.",
    "L: gseKEGG(minGSSize=10,pvalueCutoff=0.05), cached human KEGG supplied as GSON; gseaplot2(hsa03050) for both sexes; Gaussian GSVA is not substituted for GSEA.",
    "Upstream historical female/male expr_patient initialization was incomplete; the explicit AverageExpression definition at Fig.5-new.R:1072 and 1242 is used consistently.")
  writeLines(methods, file.path(config$output_dir, "source_data", "Figure6_JKL_methods.txt"))
  list(J = J, K = K, L = L)
}

# ======================== 6. A-L assembly and execution ========================
figure6_label <- function(plot, letter) {
  # wrap_elements makes nested module grids a single labelled panel.
  patchwork::wrap_elements(full = plot) + patchwork::plot_annotation(tag_levels = list(letter)) &
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 13, face = "bold"),
                   plot.tag.position = c(0, 1), plot.margin = ggplot2::margin(2, 2, 2, 2))
}

build_figure6 <- function(config = figure6_config) {
  figure6_check_packages()
  suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
  if (!length(config$formats) || any(!config$formats %in% c("pdf", "png", "svg", "tiff"))) {
    stop("formats must contain pdf, png, svg, or tiff.")
  }
  if (!is.numeric(config$topology_edge_fraction) ||
      config$topology_edge_fraction <= 0 || config$topology_edge_fraction > 1) {
    stop("topology_edge_fraction must be in (0, 1].")
  }
  paths <- unlist(config[c("female_file", "male_file", "female_umap_file", "male_umap_file",
                           "female_tom_file", "male_tom_file")], use.names = TRUE)
  if (any(!file.exists(paths))) stop("Input files not found: ",
                                     paste(paths[!file.exists(paths)], collapse = ", "))
  for (d in c(config$output_dir, file.path(config$output_dir, "panels"),
              file.path(config$output_dir, "source_data"), config$pathway_cache_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  # Fail early on a missing database/cache, before reading GB-sized network inputs.
  message("Checking MSigDB and KEGG pathway database caches...")
  invisible(f6_read_legacy_kegg(config))
  invisible(f6_read_live_kegg(config))
  set.seed(config$seed)
  saveRDS(config, file.path(config$output_dir, "Figure6_parameters.rds"))
  provenance <- data.frame(input = names(paths), path = normalizePath(paths),
                           bytes = file.info(paths)$size, modified = as.character(file.info(paths)$mtime))
  figure6_write_table(provenance, "input_manifest.csv", config)
  female <- figure6_prepare_object(config$female_file, "Female", config$female_dot_order, config)
  male <- figure6_prepare_object(config$male_file, "Male", config$male_dot_order, config)
  panels <- list(
    A = figure6_feature_grid(female, paste0("FMM", 1:13), 5),
    B = figure6_feature_grid(male, paste0("Male.MM", 1:3), 3),
    C = figure6_dotplot(female, config$female_dot_order, "ESR2_group", "C", config),
    D = figure6_dotplot(male, config$male_dot_order, "ESR2_group", "D", config),
    E = figure6_dotplot(female, config$female_dot_order, "group.a", "E", config),
    F = figure6_dotplot(male, config$male_dot_order, "group.a", "F", config))
  female_net <- figure6_network_panels(female, config$female_tom_file,
                                       config$female_umap_file, "Female", config)
  # Each TOM is local to its helper and released before loading the next one.
  invisible(gc())
  male_net <- figure6_network_panels(male, config$male_tom_file,
                                     config$male_umap_file, "Male", config)
  panels$G <- patchwork::wrap_plots(c(female_net$circles, male_net$circles), ncol = 3)
  panels$H <- female_net$topology
  panels$I <- male_net$topology
  pathway_panels <- build_figure6_pathways(female, male, config)
  panels[c("J", "K", "L")] <- pathway_panels[c("J", "K", "L")]
  stopifnot(identical(sort(names(panels)), LETTERS[1:12]))
  # Assemble the top 9 panels using the submitted figure's left/right arrangement.
  # G spans the left side below C-F; H/I form the right-hand network column.
  labelled <- stats::setNames(lapply(names(panels), function(letter)
    figure6_label(panels[[letter]], letter)), names(panels))
  layout <- c(
    patchwork::area(t = 1, l = 1, b = 2, r = 3),  # A
    patchwork::area(t = 1, l = 4, b = 1, r = 6),  # B
    patchwork::area(t = 3, l = 1, b = 3, r = 2),  # C
    patchwork::area(t = 3, l = 3, b = 3, r = 3),  # D
    patchwork::area(t = 4, l = 1, b = 4, r = 2),  # E
    patchwork::area(t = 4, l = 3, b = 4, r = 3),  # F
    patchwork::area(t = 5, l = 1, b = 6, r = 3),  # G
    patchwork::area(t = 2, l = 4, b = 4, r = 6),  # H
    patchwork::area(t = 5, l = 4, b = 6, r = 6))  # I
  top <- patchwork::wrap_plots(labelled[LETTERS[1:9]], design = layout,
                               heights = c(1, 0.7, 1, 1, 0.7, 0.7), widths = c(1, 1, 1, 1, 1, 1))
  bottom <- patchwork::wrap_plots(labelled[c("J", "K", "L")],
                                  nrow = 1, widths = c(1.1, 1.25, 2.6))
  figure <- patchwork::wrap_plots(list(top, bottom), ncol = 1, heights = c(4.5, 1.25))
  message("Exporting A-L panels and combined Figure 6...")
  # Generous standalone sizes keep the individual panel text readable.
  sizes <- list(A = c(120, 90), B = c(120, 45), C = c(100, 70), D = c(65, 70),
                E = c(100, 70), F = c(65, 70), G = c(180, 65), H = c(125, 135),
                I = c(125, 120), J = c(100, 100), K = c(110, 100), L = c(180, 100))
  for (letter in names(panels)) figure6_export(panels[[letter]],
                                               file.path(config$output_dir, "panels", paste0("Figure6", letter)),
                                               sizes[[letter]][1], sizes[[letter]][2], config)
  figure6_export(figure, file.path(config$output_dir, "Figure6_A-L"),
                 config$width_mm, config$height_mm, config)
  writeLines(capture.output(utils::sessionInfo()), file.path(config$output_dir, "sessionInfo.txt"))
  writeLines(c(
    "Figure6.R consolidates historical A-L drawing and analysis code.",
    "C/D: reconstructed drawing calls using saved cell-level ESR2_group labels.",
    "High/Low labels in these inputs correspond to positive/zero RNA/data ESR2.",
    "A/B use saved hMEs and cell UMAP; G-I use saved modules, TOM and gene UMAP.",
    "G-I use explicit seeded edge selection and new plotting code; layout is reconstructed.",
    "K/L retain patient-level summaries; C-F are pooled-cell descriptive displays.",
    "J/K use KEGG LEGACY sets; L uses KEGG pathway definitions cached at run time.",
    "Database version changes can change J-L; cache and source tables record the run.",
    "No exact reproduction of the original PDF or validation of upstream WGCNA is claimed."),
    file.path(config$output_dir, "reconstruction_notes.txt"))
  message("Completed: ", file.path(config$output_dir, "Figure6_A-L.pdf"))
  invisible(list(figure = figure, panels = panels, config = config))
}

if (!identical(Sys.getenv("FIGURE6_RUN", unset = "1"), "0")) {
  figure6_result <- build_figure6(figure6_config)
}
