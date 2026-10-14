#Figure 7
# ======================== 1. Editable parameters ========================
figure7_project_root <- "/Users/tf/e/ndmm"
figure7_b_metadata_file <- "/Users/tf/e/MMspactial/ESR2/analysis/scrna/cell_metadata_esr2.csv.gz"
figure7_ctrp_dir <- file.path(figure7_project_root,
                              "ESR2/gender/Manuscript_ESR2_gender/Sex-specific ESR2 expression patterns")
figure7_b_reference_file <- file.path(figure7_ctrp_dir,
                                      "Figures1_8/Figure7B_pairwise_pvalues.csv")
figure7_output_dir <- file.path(figure7_project_root,
                                "ESR2/submit/Communications Biology/Figure7_code/outputs")
figure7_formats <- c("pdf", "png")  # Optional: "svg", "tiff".
figure7_width_mm <- 183
figure7_height_mm <- 174
figure7_dpi <- 300
figure7_strict_reference <- TRUE  # Stop if the original B-D results do not match.
figure7_b_axis_label <- "Bortezomib Response Rate (%)"  # Historical cell-level label.

# ======================== 2. Shared helpers ========================
figure7_check_packages <- function() {
  needed <- c("ggplot2", "patchwork", "ggpubr", "binom")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". On this machine use /Library/Frameworks/R.framework/Resources/bin/Rscript.")
  if (!capabilities("cairo")) stop("An R installation with Cairo graphics is required.")
}

figure7_theme <- function() {
  ggplot2::theme_classic(base_size = 7.5, base_family = "sans") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 8, colour = "black"),
      axis.text = ggplot2::element_text(size = 7, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.3),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
      axis.ticks.length = grid::unit(1, "mm"),
      legend.title = ggplot2::element_text(size = 6.5),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = grid::unit(2.5, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.title = ggplot2::element_text(size = 9, hjust = 0.5),
      plot.tag = ggplot2::element_text(size = 13),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(5, 4, 4, 4)
    )
}

figure7_export <- function(plot, stem, width_mm, height_mm, formats, dpi) {
  for (format in formats) {
    filename <- paste0(stem, ".", format)
    export_one <- function() {
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
    export_one()
  }
}

figure7_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

# ======================== 3. Panel A: English vector flowchart ========================
figure7_plot_a <- function(text_family = "serif", base_size = 6.6) {
  blue <- "#4279D5"
  border_blue <- "#91BDE9"
  light_blue <- "#B9D4EB"
  title_blue <- "#96BEE4"
  connectors <- data.frame(
    x = c(0.50, 0.27, 0.27, 0.76, 0.27, 0.76, 0.27,
          0.50, 0.50, 0.17, 0.17, 0.50, 0.84),
    y = c(0.865, 0.835, 0.610, 0.610, 0.570, 0.835, 0.835,
          0.570, 0.385, 0.335, 0.335, 0.335, 0.335),
    xend = c(0.50, 0.76, 0.27, 0.76, 0.76, 0.76, 0.27,
             0.50, 0.50, 0.84, 0.17, 0.50, 0.84),
    yend = c(0.835, 0.835, 0.570, 0.570, 0.570, 0.790, 0.790,
             0.535, 0.335, 0.335, 0.280, 0.280, 0.280),
    arrow = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
              TRUE, FALSE, FALSE, TRUE, TRUE, TRUE)
  )
  nodes <- data.frame(
    x = c(0.50, 0.27, 0.76, 0.50, 0.17, 0.50, 0.84),
    y = c(0.91, 0.70, 0.70, 0.46, 0.17, 0.17, 0.17),
    width = c(0.48, 0.46, 0.43, 0.57, 0.30, 0.30, 0.33),
    height = c(0.09, 0.18, 0.18, 0.15, 0.22, 0.22, 0.22),
    label = c("scRNA-seq data", "ESR2 expression\nstratification",
              "Gender\nstratification", "Treatment response\nanalysis",
              "Pathway\nenrichment\nanalysis", "Correlation\nwith drug\nsensitivity",
              "Experimental\nValidation"),
    emphasis = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = connectors[!connectors$arrow, ],
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          colour = blue, linewidth = 0.16, lineend = "round") +
    ggplot2::geom_segment(data = connectors[connectors$arrow, ],
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          colour = blue, linewidth = 0.16, lineend = "round",
                          arrow = grid::arrow(length = grid::unit(0.85, "mm"), angle = 24, type = "closed")) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_size = base_size, base_family = text_family) +
    ggplot2::theme(plot.margin = grid::unit(rep(0.4, 4), "mm"))
  for (i in seq_len(nrow(nodes))) {
    node <- nodes[i, ]
    node_grob <- grid::grobTree(
      grid::roundrectGrob(r = grid::unit(if (node$emphasis) 1.2 else 1.6, "mm"),
                          gp = grid::gpar(fill = if (node$emphasis) title_blue else light_blue,
                                          col = border_blue, lwd = 0.5)),
      grid::textGrob(node$label, gp = grid::gpar(fontfamily = text_family,
                                                 fontsize = if (node$emphasis) base_size + 0.25 else base_size,
                                                 fontface = if (node$emphasis) "bold" else "plain",
                                                 col = if (node$emphasis) blue else "black", lineheight = 1.06)))
    p <- p + ggplot2::annotation_custom(node_grob,
                                        xmin = node$x - node$width / 2, xmax = node$x + node$width / 2,
                                        ymin = node$y - node$height / 2, ymax = node$y + node$height / 2)
  }
  p
}

# ======================== 4. Panel B: reproduce the original cell counts ========================
figure7_read_b <- function(path) {
  if (!file.exists(path)) stop("B metadata file not found: ", path)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  x <- utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("cell_id", "sample", "sex", "disease", "response_group", "celltype", "ESR2_count")
  if (!all(needed %in% names(x))) stop("B metadata are missing required columns: ",
                                       paste(setdiff(needed, names(x)), collapse = ", "))
  x <- x[!is.na(x$disease) & !is.na(x$celltype) & x$disease == "MM" &
           x$celltype == "Malignant plasma cells", needed, drop = FALSE]
  if (!nrow(x) || anyNA(x) || anyDuplicated(x$cell_id) ||
      !all(x$sex %in% c("Female", "Male")) ||
      !all(x$response_group %in% c("OR", "SOR")) ||
      !is.numeric(x$ESR2_count) || any(!is.finite(x$ESR2_count)) || any(x$ESR2_count < 0)) {
    stop("Invalid malignant-cell metadata for B.")
  }
  donor <- unique(x[, c("sample", "sex", "response_group")])
  if (anyDuplicated(donor$sample)) stop("A patient has inconsistent sex or response labels.")
  threshold <- stats::median(x$ESR2_count)
  # Original RNA/data expression has median zero. Raw counts and log-normalized
  # expression have identical zero/positive membership, so this is equivalent.
  if (threshold != 0) stop("B count median is no longer zero; normalized-expression grouping must be reviewed.")
  x$ESR2_group <- ifelse(x$ESR2_count > threshold, "ESR2_High", "ESR2_Low")
  x$gender <- x$sex
  x$group.a <- x$response_group
  x$patient_group <- paste(x$gender, x$ESR2_group, sep = "/")
  attr(x, "ESR2_threshold") <- threshold
  x
}

figure7_summarise_b <- function(meta) {
  groups <- expand.grid(gender = c("Female", "Male"),
                        ESR2_group = c("ESR2_High", "ESR2_Low"), stringsAsFactors = FALSE)
  summary <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    d <- meta[meta$gender == groups$gender[i] & meta$ESR2_group == groups$ESR2_group[i], ]
    if (!nrow(d)) stop("A Figure 7B group is empty.")
    total <- nrow(d)
    n_or <- sum(d$group.a == "OR")
    ci <- binom::binom.confint(n_or, total, methods = "wilson")
    data.frame(gender = groups$gender[i], ESR2_group = groups$ESR2_group[i],
               patient_group = paste(groups$gender[i], groups$ESR2_group[i], sep = "/"),
               n_total = total, n_OR = n_or, n_SOR = total - n_or,
               n_patients_represented = length(unique(d$sample)), rate = n_or / total,
               ci_lower = ci$lower, ci_upper = ci$upper,
               label = sprintf("%.1f%%", 100 * n_or / total))
  }))
  comparisons <- utils::combn(seq_len(nrow(summary)), 2, simplify = FALSE)
  pstar <- function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  pairs <- do.call(rbind, lapply(seq_along(comparisons), function(i) {
    g <- summary[comparisons[[i]], ]
    p <- stats::fisher.test(cbind(g$n_OR, g$n_SOR))$p.value
    data.frame(group1 = g$patient_group[1], group2 = g$patient_group[2],
               pstar = pstar(p), p_value_numeric = p, y.position = 0.5 + i * 0.035)
  }))
  summary$patient_group <- factor(summary$patient_group, levels = summary$patient_group)
  list(summary = summary, pairwise = pairs)
}

figure7_plot_b <- function(summary, pairwise, axis_label) {
  ggplot2::ggplot(summary, ggplot2::aes(x = patient_group, y = rate, fill = gender)) +
    ggplot2::geom_col(width = 0.5, colour = "black", linewidth = 0.25) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                           width = 0.15, linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(y = ci_upper + 0.025, label = label), size = 2.2) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                                breaks = seq(0, 0.6, by = 0.2), limits = c(0, 0.75), expand = ggplot2::expansion(mult = c(0, 0.025))) +
    ggplot2::scale_fill_manual(values = c(Female = "#FF9999", Male = "#9999FF")) +
    ggplot2::labs(x = NULL, y = axis_label, fill = "Gender") +
    figure7_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 5.2, face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 6.7), legend.position = "top") +
    ggpubr::stat_pvalue_manual(pairwise, label = "pstar", tip.length = 0.01,
                               bracket.size = 0.25, label.size = 2.4, inherit.aes = FALSE)
}

# ======================== 5. Panels C/D: stable archived sample membership ========================
figure7_group_ids <- list(
  Overall = c("AMO1", "JJN3", "KARPAS620", "EJM", "KHM1B", "KMS11", "KMS12BM",
              "KMS18", "KMS20", "KMS21BM", "KMS26", "KMS27", "KMS28BM", "KMS34", "U266",
              "SKMM2", "RPMI8226", "OPM2", "NCIH929", "MOLP8", "MOLP2", "MM1S", "LP1", "L363"),
  Female_All = c("AMO1", "JJN3", "KARPAS620", "EJM", "KMS11", "KMS12BM", "KMS20",
                 "KMS28BM", "KMS34", "NCIH929", "MOLP8", "LP1", "L363"),
  # Legacy female_ESR2_low[-6, ] excludes NCIH929 only from Female Low.
  Female_Low = c("EJM", "KMS11", "KMS20", "KMS28BM", "KMS34", "MOLP8", "LP1", "L363"),
  Female_High = c("AMO1", "JJN3", "KARPAS620", "KMS12BM"),
  Male_All = c("KHM1B", "KMS18", "KMS26", "KMS27", "U266", "RPMI8226", "OPM2", "MOLP2", "MM1S"),
  Male_Low = c("KHM1B", "KMS18", "KMS26", "OPM2", "MOLP2", "MM1S"),
  Male_High = c("KMS27", "U266", "RPMI8226")
)

figure7_read_ctrp <- function(directory) {
  filenames <- c(Overall = "df_ctrp_unique_MM.csv", Female = "df_ctrp_unique_MM_female.csv",
                 Male = "df_ctrp_unique_MM_male.csv")
  tables <- lapply(filenames, function(filename) {
    path <- file.path(directory, filename)
    if (!file.exists(path)) stop("CTRP input not found: ", path)
    d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    needed <- c("cell_line", "ESR2", "bortezomib")
    if (!all(needed %in% names(d))) stop("Missing CTRP columns in ", path)
    d <- d[, needed, drop = FALSE]
    if (anyNA(d$cell_line) || anyDuplicated(d$cell_line) ||
        !is.numeric(d$ESR2) || !is.numeric(d$bortezomib)) stop("Invalid CTRP table: ", path)
    d
  })
  all <- tables$Overall
  if (!setequal(all$cell_line, figure7_group_ids$Overall) ||
      !setequal(tables$Female$cell_line, figure7_group_ids$Female_All) ||
      !setequal(tables$Male$cell_line, figure7_group_ids$Male_All)) {
    stop("CTRP cell-line IDs do not match the archived Figure 7 groups.")
  }
  for (sex in c("Female", "Male")) {
    d <- tables[[sex]]
    matched <- all[match(d$cell_line, all$cell_line), c("ESR2", "bortezomib"), drop = FALSE]
    if (!isTRUE(all.equal(unname(as.matrix(matched)),
                          unname(as.matrix(d[, c("ESR2", "bortezomib")]))))) {
      stop("Sex-specific and overall CTRP values differ for ", sex)
    }
  }
  membership <- do.call(rbind, lapply(names(figure7_group_ids), function(group) {
    d <- all[match(figure7_group_ids[[group]], all$cell_line), , drop = FALSE]
    d$group <- group
    d$included_in_fit <- is.finite(d$ESR2) & is.finite(d$bortezomib)
    d$exclusion_reason <- ifelse(d$included_in_fit, "", "Missing/nonfinite ESR2 or bortezomib AUC")
    d
  }))
  membership <- membership[, c("group", "cell_line", "ESR2", "bortezomib", "included_in_fit", "exclusion_reason")]
  groups <- lapply(names(figure7_group_ids), function(group) {
    d <- membership[membership$group == group & membership$included_in_fit, , drop = FALSE]
    if (nrow(d) < 3 || stats::sd(d$ESR2) == 0 || stats::sd(d$bortezomib) == 0) {
      stop("Insufficient variation or observations for ", group)
    }
    d
  })
  names(groups) <- names(figure7_group_ids)
  statistics <- do.call(rbind, lapply(names(groups), function(group) {
    d <- groups[[group]]
    ct <- stats::cor.test(d$ESR2, d$bortezomib, method = "pearson", alternative = "two.sided")
    # cor.test does not provide a correlation CI when n=3 (Male High).
    ci <- if (length(ct$conf.int) == 2L) ct$conf.int else c(NA_real_, NA_real_)
    data.frame(group = group, n = nrow(d), r = unname(ct$estimate), p = ct$p.value,
               r_ci_lower = ci[1], r_ci_upper = ci[2])
  }))
  list(groups = groups, membership = membership, statistics = statistics,
       files = stats::setNames(file.path(directory, filenames), names(filenames)))
}

figure7_plot_correlation <- function(data, title = NULL) {
  # Place r/p above the regression line using the full confidence-band range.
  # stat_cor's default NPC position only uses the observed y range and overlaps
  # the nearly flat Male High regression line in the historical artwork.
  fit <- stats::lm(bortezomib ~ ESR2, data = data)
  prediction <- stats::predict(fit, newdata = data.frame(
    ESR2 = seq(min(data$ESR2), max(data$ESR2), length.out = 80)),
    interval = "confidence", level = 0.95)
  y_range <- range(c(data$bortezomib, prediction), finite = TRUE)
  label_y <- max(y_range) - 0.035 * diff(y_range)
  ggplot2::ggplot(data, ggplot2::aes(x = ESR2, y = bortezomib)) +
    ggplot2::geom_point(alpha = 0.8, size = 1.15, colour = "#FE00FA") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         level = 0.95, colour = "blue", linewidth = 0.55) +
    ggpubr::stat_cor(method = "pearson", label.x.npc = 0.02, label.y = label_y,
                     r.digits = 3, p.digits = 3, size = 2.55) +
    ggplot2::labs(x = "ESR2 expression", y = "Bortezomib AUC", title = title) +
    figure7_theme()
}

figure7_row_label <- function(label) {
  ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5,
                                        label = label, size = 2.7) + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
}

# ======================== 6. Reference verification ========================
figure7_validate <- function(b, cd, reference_file, strict) {
  checks <- data.frame(check = character(), passed = logical(), detail = character())
  add_check <- function(name, passed, detail) {
    checks <<- rbind(checks, data.frame(check = name, passed = isTRUE(passed), detail = detail))
  }
  add_check("B four OR counts", identical(as.integer(b$summary$n_OR), c(137L, 171L, 6233L, 6517L)),
            paste(b$summary$n_OR, collapse = ", "))
  add_check("B four total cell counts", identical(as.integer(b$summary$n_total), c(486L, 869L, 13592L, 22244L)),
            paste(b$summary$n_total, collapse = ", "))
  add_check("B four percentages", identical(b$summary$label, c("28.2%", "19.7%", "45.9%", "29.3%")),
            paste(b$summary$label, collapse = ", "))
  if (file.exists(reference_file)) {
    ref <- utils::read.csv(reference_file, stringsAsFactors = FALSE)
    needed <- c("group1", "group2", "pstar", "p_value_numeric")
    if (!all(needed %in% names(ref))) stop("Invalid B pairwise reference file.")
    key <- function(d) paste(d$group1, d$group2, sep = "|")
    idx <- match(key(b$pairwise), key(ref))
    matched <- !anyNA(idx) && nrow(ref) == 6
    add_check("B six historical Fisher tests", matched &&
                all(abs(b$pairwise$p_value_numeric / ref$p_value_numeric[idx] - 1) < 1e-8) &&
                identical(b$pairwise$pstar, ref$pstar[idx]), "Compared with the archived pairwise CSV")
  } else {
    add_check("B six historical Fisher tests", FALSE, paste("Reference not found:", reference_file))
  }
  expected <- data.frame(group = names(figure7_group_ids),
                         n = c(23L, 13L, 8L, 4L, 8L, 5L, 3L),
                         r = c(0.737, 0.806, -0.708, 0.979, 0.0378, 0.135, 0.395),
                         p = c(5.97e-5, 0.000870, 0.0496, 0.0208, 0.929, 0.829, 0.742))
  for (i in seq_len(nrow(expected))) {
    actual <- cd$statistics[match(expected$group[i], cd$statistics$group), ]
    add_check(paste("C/D", expected$group[i]), actual$n == expected$n[i] &&
                abs(signif(actual$r, 3) - expected$r[i]) < 1e-12 &&
                abs(signif(actual$p, 3) - expected$p[i]) < 1e-12,
              sprintf("n=%d, r=%.10g, p=%.10g", actual$n, actual$r, actual$p))
  }
  if (strict && !all(checks$passed)) {
    stop("Reference checks failed: ", paste(checks$check[!checks$passed], collapse = "; "))
  }
  checks
}

# ======================== 7. Build, save panels, and export source data ========================
build_figure7 <- function(b_metadata_file, ctrp_dir, b_reference_file, output_dir,
                          formats = c("pdf", "png"), width_mm = 183, height_mm = 174, dpi = 300,
                          strict_reference = TRUE, b_axis_label = "Bortezomib Response Rate (%)") {
  figure7_check_packages()
  if (!length(formats) || any(!formats %in% c("pdf", "png", "svg", "tiff"))) {
    stop("Supported export formats: pdf, png, svg, tiff.")
  }
  if (any(!is.finite(c(width_mm, height_mm, dpi))) || any(c(width_mm, height_mm, dpi) <= 0)) {
    stop("Figure dimensions and DPI must be positive finite numbers.")
  }
  message("Reading lightweight cell metadata and archived CTRP tables...")
  b_meta <- figure7_read_b(b_metadata_file)
  b <- figure7_summarise_b(b_meta)
  cd <- figure7_read_ctrp(ctrp_dir)
  checks <- figure7_validate(b, cd, b_reference_file, strict_reference)
  print(b$summary[, c("patient_group", "n_OR", "n_total", "label")], row.names = FALSE)
  print(cd$statistics[, c("group", "n", "r", "p")], row.names = FALSE)
  
  a_plot <- figure7_plot_a()
  b_plot <- figure7_plot_b(b$summary, b$pairwise, b_axis_label)
  c_plot <- figure7_plot_correlation(cd$groups$Overall)
  d_names <- setdiff(names(cd$groups), "Overall")
  d_plots <- lapply(d_names, function(group) figure7_plot_correlation(cd$groups[[group]]))
  names(d_plots) <- d_names
  headings <- c(Female_All = "ESR2-All", Female_Low = "ESR2-Low", Female_High = "ESR2-High")
  for (nm in names(headings)) d_plots[[nm]] <- d_plots[[nm]] + ggplot2::labs(title = headings[[nm]])
  d_plot <- patchwork::wrap_plots(
    figure7_row_label("female"), d_plots$Female_All, d_plots$Female_Low, d_plots$Female_High,
    figure7_row_label("male"), d_plots$Male_All, d_plots$Male_Low, d_plots$Male_High,
    nrow = 2, byrow = TRUE, widths = c(0.19, 1, 1, 1))
  tag_style <- ggplot2::theme(plot.tag = ggplot2::element_text(size = 13, family = "sans"),
                              plot.tag.position = c(0, 1))
  top <- patchwork::wrap_plots(a_plot + ggplot2::labs(tag = "A") + tag_style,
                               b_plot + ggplot2::labs(tag = "B"),
                               patchwork::wrap_elements(full = c_plot + ggplot2::labs(tag = "C")),
                               nrow = 1, widths = c(1.12, 1.08, 1))
  bottom <- patchwork::wrap_elements(full = d_plot) + ggplot2::labs(tag = "D") + tag_style
  integrated <- patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(1, 1.72))
  
  panel_dir <- file.path(output_dir, "panels")
  source_dir <- file.path(output_dir, "source_data")
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  figure7_write_csv(b$summary, file.path(source_dir, "Figure7B_cell_counts_Wilson_CI.csv"))
  figure7_write_csv(b$pairwise, file.path(source_dir, "Figure7B_pairwise_Fisher_tests.csv"))
  figure7_write_csv(cd$membership, file.path(source_dir, "Figure7CD_group_membership.csv"))
  figure7_write_csv(cd$statistics, file.path(source_dir, "Figure7CD_Pearson_statistics.csv"))
  figure7_write_csv(checks, file.path(source_dir, "Figure7_reference_checks.csv"))
  donor_counts <- unique(b_meta[, c("sample", "gender", "group.a")])
  figure7_write_csv(donor_counts, file.path(source_dir, "Figure7B_patient_metadata.csv"))
  # Preserve the explicit legacy exclusion that does not appear in Female Low.
  excluded <- data.frame(group = "Female_Low", cell_line = "NCIH929",
                         reason = "Legacy female_ESR2_low[-6, ]; retained in Female_All, excluded only from Female_Low")
  figure7_write_csv(excluded, file.path(source_dir, "Figure7CD_legacy_exclusions.csv"))
  inputs <- c(B_cell_metadata = b_metadata_file, cd$files, B_pairwise_reference = b_reference_file)
  existing <- inputs[file.exists(inputs)]
  info <- file.info(existing)
  manifest <- data.frame(input = names(existing), path = normalizePath(existing),
                         bytes = info$size, md5 = unname(tools::md5sum(existing)), stringsAsFactors = FALSE)
  figure7_write_csv(manifest, file.path(source_dir, "Figure7_input_manifest.csv"))
  
  message("Exporting Figure7_AD and individual A-D panels...")
  figure7_export(integrated, file.path(output_dir, "Figure7_AD"), width_mm, height_mm, formats, dpi)
  for (nm in c("A", "B", "C")) {
    plot <- switch(nm, A = a_plot, B = b_plot, C = c_plot)
    figure7_export(plot, file.path(panel_dir, paste0("Figure7", nm)), 63, 63, formats, dpi)
  }
  figure7_export(d_plot, file.path(panel_dir, "Figure7D"), width_mm, 110, formats, dpi)
  for (nm in names(d_plots)) {
    figure7_export(d_plots[[nm]], file.path(panel_dir, paste0("Figure7D_", nm)), 58, 53, formats, dpi)
  }
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  writeLines(c(
    "# Figure 7 A-D reproduction", "",
    "Run: /Library/Frameworks/R.framework/Resources/bin/Rscript --vanilla Figure7.R", "",
    "Figure7_AD.pdf/png contains A-D only; panels/ contains the individual editable panels.",
    "A is a vector reconstruction of the submitted English flowchart, not recovered artwork.",
    "B reproduces the historical cell-level counts: 37191 malignant cells from 18 patients.",
    "ESR2 count median is zero; >0 membership equals >median of normalized expression.",
    "B uses Wilson 95% confidence intervals and six unadjusted Fisher tests.",
    "Its historical Bortezomib Response Rate label describes pooled OR-labelled cells;",
    "it must not be interpreted as patient response rates or patient-level inference.",
    "C/D use Pearson correlation and ordinary linear regression with 95% confidence bands.",
    "C has n=23. D uses Female All/Low/High n=13/8/4 and Male All/Low/High n=8/5/3.",
    "Archived group labels are preserved, including MOLP8 in Female and OPM2/MM1S in Male.",
    "MM1S has missing AUC and is omitted from fits. Female Low excludes NCIH929, as in the legacy code.",
    "KMS21BM and SKMM2 appear in C but not in D. No sex labels or exclusions are silently corrected.",
    "source_data/ records memberships, exclusions, statistics, input hashes, and reference checks.",
    "Original input files and the submitted Figure 7.pdf are unchanged. E/F are outside this script."
  ), file.path(output_dir, "README.txt"))
  message("Saved: ", file.path(output_dir, "Figure7_AD.pdf"))
  invisible(list(A = a_plot, B = b_plot, C = c_plot, D = d_plot, integrated = integrated,
                 B_summary = b$summary, B_pairwise = b$pairwise, CD_statistics = cd$statistics,
                 CD_membership = cd$membership, checks = checks, output_dir = output_dir))
}

# ======================== 8. Run all panels in a clean or interactive session ========================
if (!identical(Sys.getenv("FIGURE7_RUN", unset = "1"), "0")) {
  figure7_result <- build_figure7(
    b_metadata_file = figure7_b_metadata_file, ctrp_dir = figure7_ctrp_dir,
    b_reference_file = figure7_b_reference_file, output_dir = figure7_output_dir,
    formats = figure7_formats, width_mm = figure7_width_mm, height_mm = figure7_height_mm,
    dpi = figure7_dpi, strict_reference = figure7_strict_reference,
    b_axis_label = figure7_b_axis_label)
}
