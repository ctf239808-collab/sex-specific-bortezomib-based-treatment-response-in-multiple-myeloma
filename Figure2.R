# =============================================================================
# Figure21.new.R
# 根据 Figure.2.pdf 的 A-I 面板整理自 Figure2.new.R。
# 文件名中的 21 表示本次新版；图中内容仍对应 Figure 2。
#
# 运行（本机已安装依赖的 R）：
# /Library/Frameworks/R.framework/Resources/bin/Rscript '/Users/tf/e/ndmm/ESR2/submit/Communications Biology/Figure21.new.R'
# 或在 RStudio 中 source() 本文件；结束后可直接查看 p2A ... p2I、p2_all。
#
# 可选设置，在 source() 前执行：
# options(figure21.input = '/absolute/path/scRNA_harmony_HM38.anno.rds',
#         figure21.output = '/absolute/path/output', figure21.dpi = 300,
#         figure21.formats = c('pdf', 'png'))
# 排版参数：figure21.width=13.2、figure21.height=13.4、figure21.font_size=11.5。
# 组合图使用三行网格、A/B 共用色条、独立字母边距和对齐的底部绘图区。
#
# 定义保持与原 R 脚本一致：
# - 所有表达来自 RNA/data（log-normalized），保留 ESR2 == 0 的细胞。
# - D：各细胞类型内，按已有性别/疾病组平均表达计算 Z-score。
# - E：细胞级双侧 Wilcoxon；F/G：细胞级双侧 Welch t 检验。
# - H/I：患者级手工评分，OR=1、SOR=0，性别内标准化，ROC direction='auto'。
#   mono_n 是计数；CD8_cyto 是全对象四基因 AddModuleScore 后的患者 CD8 均值。
# - E/F/G 用于重现原图的探索性细胞级比较，并非独立患者重复检验。
# - H/I 为同队列探索性 ROC，并非 LASSO、交叉验证或独立验证结果。
# 不根据参考图修改数据、P 值或 AUC；参考标记仅用于导出的核对表。
# 本机实测：H/I 六个 AUC 均匹配参考图（两位小数）。
# E 中 MM 恶性浆细胞的男女 Wilcoxon P=0.0696273（ns），与 PDF 的 *** 不同。
# =============================================================================

# 01. 参数、依赖与公共函数 ------------------------------------------------------
figure21_cfg <- list(
  input = getOption('figure21.input',
                    '/Users/tf/e/ndmm/HM38/scRNA_harmony_HM38.anno.rds'),
  out_dir = getOption('figure21.output', paste0(
    '/Users/tf/e/ndmm/ESR2/submit/Communications Biology/', 'Figure21_output')),
  formats = getOption('figure21.formats', c('pdf', 'png')),
  dpi = getOption('figure21.dpi', 300),
  figure_width = getOption('figure21.width', 13.2),
  figure_height = getOption('figure21.height', 13.4),
  font_size = getOption('figure21.font_size', 11.5),
  font_family = 'sans',
  gene = 'ESR2', assay = 'RNA', seed = 20261013L,
  module_seed = 1L, module_nbin = 24L, module_ctrl = 100L,
  cyto_genes = c('GZMB', 'PRF1', 'NKG7', 'GNLY'),
  mono_types = c('Macrophages', 'CD14+ monocytes', 'CD16+ Monocytes',
                 'Intermediate monocytes'),
  roc_direction = 'auto'
)
figure21_packages <- c('Seurat', 'ggplot2', 'dplyr', 'tidyr', 'patchwork',
                       'scales', 'pROC')
figure21_missing <- figure21_packages[!vapply(
  figure21_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(figure21_missing)) {
  stop('Missing R packages: ', paste(figure21_missing, collapse = ', '),
       '. On this machine use /Library/Frameworks/R.framework/Resources/bin/Rscript.')
}
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})
if (!length(figure21_cfg$formats) ||
    any(!figure21_cfg$formats %in% c('pdf', 'png'))) {
  stop('figure21.formats must contain pdf and/or png.')
}
set.seed(figure21_cfg$seed)

figure21_theme <- function(base_size = figure21_cfg$font_size) {
  theme_classic(base_size = base_size, base_family = figure21_cfg$font_family) +
    theme(
      text = element_text(colour = '#252525'),
      axis.text = element_text(size = base_size - 1.5, colour = '#333333'),
      axis.title = element_text(size = base_size, face = 'plain'),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.line = element_line(colour = '#454545', linewidth = 0.35),
      axis.ticks = element_line(colour = '#454545', linewidth = 0.3),
      axis.ticks.length = grid::unit(2.4, 'pt'),
      plot.title = element_text(size = base_size + 0.5, hjust = 0.5,
                                face = 'bold', margin = margin(b = 7)),
      plot.margin = margin(6, 7, 6, 7),
      legend.position = 'bottom', legend.box = 'vertical',
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1.5),
      legend.key = element_rect(fill = 'white', colour = NA),
      legend.key.size = grid::unit(9, 'pt'),
      legend.spacing.y = grid::unit(3, 'pt'),
      legend.margin = margin(3, 0, 0, 0),
      strip.background = element_rect(fill = '#F1F3F5', colour = NA),
      strip.text = element_text(size = base_size - 1, colour = '#252525',
                               margin = margin(4, 3, 4, 3)),
      panel.spacing = grid::unit(8, 'pt')
    )
}
figure21_stars <- function(p, four_stars = FALSE) {
  if (!is.finite(p)) return('NA')
  if (four_stars && p <= 0.0001) return('****')
  if (p <= 0.001) return('***')
  if (p <= 0.01) return('**')
  if (p <= 0.05) return('*')
  'ns'
}
figure21_row_z <- function(x) {
  s <- stats::sd(x)
  if (length(x) < 2L || !is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
figure21_save <- function(plot, name, width, height) {
  for (ext in unique(figure21_cfg$formats)) {
    ggplot2::ggsave(filename = file.path(figure21_cfg$out_dir, paste0(name, '.', ext)),
                   plot = plot, width = width, height = height, units = 'in',
                   dpi = figure21_cfg$dpi, bg = 'white', limitsize = FALSE,
                   # Cairo PDF embeds fonts and preserves editable vector text.
                   device = if (ext == 'pdf' && capabilities('cairo'))
                     grDevices::cairo_pdf else ext)
  }
}

# 02. 只读取一次 Seurat 对象，统一表达来源与细胞顺序 ----------------------------
if (!file.exists(figure21_cfg$input)) stop('Input file not found: ', figure21_cfg$input)
message('Reading Seurat object: ', figure21_cfg$input)
scRNA_harmony <- readRDS(figure21_cfg$input)
if (!inherits(scRNA_harmony, 'Seurat')) stop('Input must be a Seurat object.')
if (!figure21_cfg$assay %in% Assays(scRNA_harmony)) stop('RNA assay is missing.')
DefaultAssay(scRNA_harmony) <- figure21_cfg$assay
if (!'umap' %in% Reductions(scRNA_harmony)) stop('UMAP reduction is missing.')

figure21_meta <- scRNA_harmony[[]]
figure21_required <- c('orig.ident', 'anno', 'gender', 'group', 'group.a')
if (!all(figure21_required %in% names(figure21_meta))) {
  stop('Missing metadata: ', paste(setdiff(figure21_required, names(figure21_meta)),
                                  collapse = ', '))
}
figure21_rna <- if (packageVersion('SeuratObject') >= '5.0.0') {
  GetAssayData(scRNA_harmony, assay = figure21_cfg$assay, layer = 'data')
} else {
  GetAssayData(scRNA_harmony, assay = figure21_cfg$assay, slot = 'data')
}
figure21_needed <- c(figure21_cfg$gene, figure21_cfg$cyto_genes)
if (!all(figure21_needed %in% rownames(figure21_rna))) {
  stop('Missing genes: ', paste(setdiff(figure21_needed, rownames(figure21_rna)),
                              collapse = ', '))
}
if (!all(rownames(figure21_meta) %in% colnames(figure21_rna))) {
  stop('RNA/data cell names do not match metadata.')
}
figure21_cells <- figure21_meta[, figure21_required, drop = FALSE] %>%
  mutate(across(everything(), as.character))
figure21_cells$cell <- rownames(figure21_meta)
figure21_cells$ESR2 <- as.numeric(figure21_rna[figure21_cfg$gene, figure21_cells$cell])
if (anyNA(figure21_cells[figure21_required]) ||
    any(!is.finite(figure21_cells$ESR2)) || any(figure21_cells$ESR2 < 0)) {
  stop('Missing metadata or invalid RNA/data ESR2 expression.')
}
if (any(!figure21_cells$gender %in% c('Female', 'Male')) ||
    any(!figure21_cells$group %in% c('HD', 'MM'))) {
  stop('Expected gender Female/Male and group HD/MM.')
}
figure21_sample_check <- figure21_cells %>%
  group_by(orig.ident) %>%
  summarise(across(c(gender, group, group.a), n_distinct), .groups = 'drop')
if (any(as.matrix(figure21_sample_check[-1]) != 1L)) {
  stop('Each orig.ident must have one gender, disease status and response label.')
}
if (any(!figure21_cells$group.a[figure21_cells$group == 'MM'] %in% c('OR', 'SOR'))) {
  stop('MM group.a must contain OR/SOR.')
}
figure21_celltypes <- sort(unique(figure21_cells$anno))
figure21_sexes <- c('Female', 'Male')
figure21_groups <- c('HD', 'MM')
figure21_umap <- Embeddings(scRNA_harmony, reduction = 'umap')
if (!all(figure21_cells$cell %in% rownames(figure21_umap)) ||
    any(!is.finite(figure21_umap))) stop('Invalid UMAP coordinates.')
figure21_umap <- figure21_umap[figure21_cells$cell, 1:2, drop = FALSE]
figure21_xlim <- range(figure21_umap[, 1]) + c(-1, 1)
figure21_ylim <- range(figure21_umap[, 2]) + c(-1, 1)
rm(figure21_rna, figure21_meta)

# 03. A：全部细胞的 ESR2 UMAP；B：HD/MM 分组，共用坐标范围和色标 -----------------
figure21_feature <- function(cells = NULL, title = 'ESR2') {
  # FeaturePlot 自带颜色标尺；下面用原图的五色标尺统一替换。
  p <- Seurat::FeaturePlot(
    scRNA_harmony, features = figure21_cfg$gene, cells = cells,
    reduction = 'umap', slot = 'data', pt.size = 0.35, order = TRUE,
    raster = FALSE, combine = FALSE
  )[[1]]
  p <- suppressMessages(p + scale_colour_gradientn(
    colours = c('#C9C9C9', '#FFFF56', '#EE9438', '#E73726', '#6E140C'),
    values = c(0, 0.25, 0.50, 0.73, 1), limits = c(0, 3), breaks = 0:3,
    oob = scales::squish, name = 'ESR2 expression',
    guide = guide_colourbar(direction = 'horizontal', title.position = 'top',
                            title.hjust = 0.5, barwidth = grid::unit(32, 'mm'),
                            barheight = grid::unit(2.5, 'mm'))))
  p + coord_cartesian(xlim = figure21_xlim, ylim = figure21_ylim) +
    labs(title = title, x = 'UMAP 1', y = 'UMAP 2') + figure21_theme() +
    theme(axis.text = element_text(size = figure21_cfg$font_size - 2))
}
p2A <- figure21_feature()
p2B_HD <- figure21_feature(figure21_cells$cell[figure21_cells$group == 'HD'], 'HD')
p2B_MM <- figure21_feature(figure21_cells$cell[figure21_cells$group == 'MM'], 'MM')
p2B <- (p2B_HD + p2B_MM + plot_layout(ncol = 2, guides = 'collect')) &
  theme(legend.position = 'bottom')

# 04. C：各细胞类型的 ESR2 分布，按 HD/MM 分组 ---------------------------------
scRNA_harmony$anno <- factor(scRNA_harmony$anno, levels = figure21_celltypes)
scRNA_harmony$group <- factor(scRNA_harmony$group, levels = figure21_groups)
figure21_vln_args <- list(
  object = scRNA_harmony, features = figure21_cfg$gene, assay = figure21_cfg$assay,
  group.by = 'anno', split.by = 'group', cols = c('#F8766D', '#00BFC4'),
  pt.size = 0.0001, split.plot = FALSE, combine = FALSE, raster = FALSE)
if ('layer' %in% names(formals(Seurat::VlnPlot))) {
  figure21_vln_args$layer <- 'data'
} else {
  figure21_vln_args$slot <- 'data'
}
p2C <- do.call(Seurat::VlnPlot, figure21_vln_args)[[1]] +
  labs(title = NULL, x = NULL, y = 'ESR2 expression', fill = NULL) +
  figure21_theme() +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1,
                                   size = figure21_cfg$font_size - 2.5),
        legend.box = 'horizontal', legend.margin = margin(5, 0, 0, 0))
rm(figure21_vln_args)

# 05. D：表达 Z-score 热图 + ESR2 检出率点图 -----------------------------------
figure21_conditions <- c('Female | HD', 'Male | HD', 'Female | MM', 'Male | MM')
figure21_D_observed <- figure21_cells %>%
  group_by(anno, gender, group) %>%
  summarise(n_cells = n(), mean_expr = mean(ESR2), detect_rate = mean(ESR2 > 0),
            .groups = 'drop') %>%
  mutate(condition = paste(gender, group, sep = ' | ')) %>%
  group_by(anno) %>% mutate(mean_expr_z = figure21_row_z(mean_expr)) %>% ungroup()
figure21_D_data <- tidyr::expand_grid(anno = figure21_celltypes,
                                     condition = figure21_conditions) %>%
  left_join(figure21_D_observed, by = c('anno', 'condition')) %>%
  mutate(anno = factor(anno, levels = figure21_celltypes),
         condition = factor(condition, levels = figure21_conditions))
p2D <- ggplot(figure21_D_data, aes(condition, anno)) +
  # 未观察到细胞的组合留白；已观察到但组间无变异的行以灰色显示。
  geom_tile(data = filter(figure21_D_data, !is.na(n_cells)), aes(fill = mean_expr_z),
            colour = 'grey85', linewidth = 0.3) +
  geom_point(data = filter(figure21_D_data, !is.na(n_cells)), aes(size = detect_rate),
             shape = 21, fill = 'black', colour = 'black', stroke = 0.2) +
  scale_fill_gradient2(low = '#3B6FB6', mid = 'white', high = '#C7362A',
                       midpoint = 0, limits = c(-1.5, 1.5),
                       breaks = c(-1.5, 0, 1.5), oob = scales::squish,
                       na.value = 'grey50', name = 'ESR2 expression (Z-score)') +
  scale_size_continuous(range = c(0.6, 4.2),
                        limits = c(0, max(0.05, figure21_D_observed$detect_rate)),
                        breaks = seq(0, 0.05, 0.01),
                        labels = scales::label_percent(accuracy = 1),
                        name = 'Detection rate') +
  scale_x_discrete(drop = FALSE, labels = function(x) gsub(' | ', '\n', x, fixed = TRUE),
                    expand = expansion(add = 0)) +
  scale_y_discrete(drop = FALSE, expand = expansion(add = 0)) +
  labs(x = NULL, y = NULL) + figure21_theme() +
  guides(fill = guide_colourbar(order = 1, direction = 'horizontal',
                                title.position = 'top', title.hjust = 0,
                                barwidth = grid::unit(42, 'mm'),
                                barheight = grid::unit(2.5, 'mm')),
         size = guide_legend(order = 2, nrow = 1, title.position = 'top',
                              override.aes = list(colour = 'black', fill = 'black'))) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = figure21_cfg$font_size - 1.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 1.1,
                                   margin = margin(t = 7)),
        legend.box.just = 'left', legend.justification = 'left',
        legend.key.width = grid::unit(10, 'pt'))

# 06. E：四类重点细胞，HD/MM 内的女性与男性比较 --------------------------------
figure21_key_cells <- c('CD4+ T cells', 'CD8+ T cells', 'Malignant plasma cells',
                        'NK cells')
figure21_E_data <- figure21_cells %>% filter(anno %in% figure21_key_cells) %>%
  mutate(anno = factor(anno, levels = figure21_key_cells),
         gender = factor(gender, levels = figure21_sexes),
         group = factor(group, levels = figure21_groups))
figure21_E_stats <- figure21_E_data %>% group_by(group, anno, .drop = TRUE) %>%
  group_modify(function(.x, .y) {
    x <- .x$ESR2[.x$gender == 'Female']
    y <- .x$ESR2[.x$gender == 'Male']
    # 原 stat_compare_means 的两组默认方法：非配对、双侧 Wilcoxon。
    # 单细胞表达有大量相同值，显式使用渐近检验和连续性校正。
    p <- if (length(x) && length(y) && length(unique(c(x, y))) > 1L) {
      stats::wilcox.test(x, y, alternative = 'two.sided', paired = FALSE,
                         exact = FALSE, correct = TRUE)$p.value
    } else NA_real_
    data.frame(n_Female_cells = length(x), n_Male_cells = length(y),
               p_value = p, label = figure21_stars(p, four_stars = TRUE))
  }) %>% ungroup()
figure21_E_limits <- figure21_E_data %>% group_by(group) %>%
  summarise(y_max = max(ESR2), .groups = 'drop') %>%
  mutate(y_bracket = y_max * 1.06, y_tip = y_max * 1.04, y_text = y_max * 1.10)
figure21_E_reference <- data.frame(
  group = c(rep('HD', 3), rep('MM', 4)),
  anno = c('CD4+ T cells', 'CD8+ T cells', 'NK cells', figure21_key_cells),
  reference_label = c('****', '****', '**', 'ns', 'ns', '***', 'ns'))
figure21_E_stats <- figure21_E_stats %>%
  left_join(figure21_E_limits, by = 'group') %>%
  left_join(figure21_E_reference, by = c('group', 'anno')) %>%
  mutate(matches_reference = label == reference_label)
if (any(!figure21_E_stats$matches_reference, na.rm = TRUE)) {
  warning('Panel E significance differs from the reference PDF; actual P values are retained. See Figure21_E_cell_level_Wilcoxon_tests.csv.',
           call. = FALSE)
}
p2E <- ggplot(figure21_E_data, aes(gender, ESR2, fill = gender)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.2, height = 0,
                                       seed = figure21_cfg$seed),
             size = 0.32, alpha = 0.4) +
  geom_segment(data = filter(figure21_E_stats, is.finite(p_value)),
                aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
                linewidth = 0.35, inherit.aes = FALSE) +
  geom_text(data = figure21_E_stats, aes(x = 1.5, y = y_text, label = label),
             inherit.aes = FALSE, size = 3.2) +
  facet_grid(group ~ anno, scales = 'free_y', drop = FALSE,
              labeller = labeller(anno = as_labeller(c(
                'CD4+ T cells' = 'CD4+ T cells', 'CD8+ T cells' = 'CD8+ T cells',
                'Malignant plasma cells' = 'Malignant plasma\ncells', 'NK cells' = 'NK cells')))) +
  scale_fill_manual(values = c(Female = '#F8766D', Male = '#00BFC4'), drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = NULL, y = 'ESR2 expression', fill = NULL) + figure21_theme() +
  guides(fill = guide_legend(nrow = 1)) +
  theme(strip.text.x = element_text(size = figure21_cfg$font_size - 2, lineheight = 1),
        strip.text.y = element_text(angle = 0, face = 'bold'),
        axis.text.x = element_text(size = figure21_cfg$font_size - 2),
        panel.spacing.x = grid::unit(5, 'pt'), panel.spacing.y = grid::unit(10, 'pt'))

# 07. F/G：MM 恶性浆细胞，OR/SOR 总体和按性别分层 ------------------------------
figure21_mpc <- figure21_cells %>%
  filter(group == 'MM', anno == 'Malignant plasma cells') %>%
  mutate(group.a = factor(group.a, levels = c('OR', 'SOR')),
         gender = factor(gender, levels = figure21_sexes))
if (!nrow(figure21_mpc)) stop('No MM malignant plasma cells.')
figure21_ttest <- function(d) {
  counts <- table(d$group.a)
  if (any(counts[c('OR', 'SOR')] < 2L)) stop('Welch test requires >=2 cells per group.')
  p <- stats::t.test(ESR2 ~ group.a, data = d, alternative = 'two.sided',
                     var.equal = FALSE)$p.value
  data.frame(n_OR_cells = sum(d$group.a == 'OR'),
             n_SOR_cells = sum(d$group.a == 'SOR'), p_value = p,
             label = figure21_stars(p))
}
figure21_F_stats <- figure21_ttest(figure21_mpc)
figure21_G_stats <- bind_rows(lapply(figure21_sexes, function(sex) {
  cbind(gender = sex, figure21_ttest(filter(figure21_mpc, gender == sex)))
}))
figure21_ymax <- max(3, ceiling(max(figure21_mpc$ESR2) * 1.1))
figure21_violin <- function(d, stats, facet = FALSE) {
  stats <- mutate(stats, y_bracket = figure21_ymax * 0.95,
                  y_tip = figure21_ymax * 0.93, y_text = figure21_ymax * 0.98)
  p <- ggplot(d, aes(group.a, ESR2)) +
    geom_violin(fill = '#F1F3F5', colour = '#404040', linewidth = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = 'white', linewidth = 0.35) +
    geom_point(position = position_jitter(width = 0.18, height = 0,
                                         seed = figure21_cfg$seed),
               size = 0.2, alpha = 0.4, colour = '#252525') +
    geom_segment(data = stats, aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
                  inherit.aes = FALSE) +
    geom_segment(data = stats, aes(x = 1, xend = 1, y = y_tip, yend = y_bracket),
                  inherit.aes = FALSE) +
    geom_segment(data = stats, aes(x = 2, xend = 2, y = y_tip, yend = y_bracket),
                  inherit.aes = FALSE) +
    geom_text(data = stats, aes(x = 1.5, y = y_text, label = label),
               inherit.aes = FALSE, size = 3.3) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = seq(0, figure21_ymax, 1),
                        expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(-0.1, figure21_ymax * 1.06)) +
    labs(x = NULL, y = 'ESR2 expression') + figure21_theme() +
    theme(panel.border = element_rect(colour = '#454545', fill = NA, linewidth = 0.35),
          axis.line = element_blank(), strip.background = element_blank())
  if (facet) p <- p + facet_wrap(~gender, nrow = 1, drop = FALSE)
  p
}
p2F <- figure21_violin(figure21_mpc, figure21_F_stats) + labs(title = 'All MM')
p2G <- figure21_violin(figure21_mpc, figure21_G_stats, facet = TRUE)
figure21_patients <- figure21_mpc %>%
  group_by(orig.ident, gender, group.a) %>%
  summarise(ESR2_mean = mean(ESR2), n_MPC = n(), .groups = 'drop')

# 08. H/I：患者级探索性三评分 ROC；OR = 1，SOR = 0，沿用原脚本自动方向。
# mono_n 是细胞计数；CD8_cyto 是完整对象 AddModuleScore 后的患者 CD8 均值。
# 此处不拟合 LASSO、不进行交叉验证；所有 AUC 均由实际数据计算。
# -----------------------------------------------------------------------------
if (anyDuplicated(figure21_patients$orig.ident)) stop("ROC 患者编号重复。")
if (any(!as.character(figure21_patients$gender) %in% c("Female", "Male")) ||
    any(!as.character(figure21_patients$group.a) %in% c("OR", "SOR"))) {
  stop("ROC 患者性别或 OR/SOR 标签无效。")
}
if (!figure21_cfg$roc_direction %in% c("auto", "<", ">")) {
  stop("roc_direction 必须是 auto、< 或 >。")
}
figure21_missing_genes <- setdiff(
  figure21_cfg$cyto_genes, rownames(scRNA_harmony[[figure21_cfg$assay]])
)
if (length(figure21_missing_genes)) {
  stop("CD8 模块缺少基因：", paste(figure21_missing_genes, collapse = ", "))
}
if (!any(figure21_cells$anno == "CD8+ T cells")) stop("未找到 CD8+ T cells。")
figure21_mono <- figure21_cells %>%
  filter(anno %in% figure21_cfg$mono_types) %>%
  count(orig.ident, name = "mono_n")

# 基因池、分箱和对照抽样作用于全部细胞，之后才提取 CD8 T 细胞。
figure21_scored <- Seurat::AddModuleScore(
  object = scRNA_harmony, features = list(figure21_cfg$cyto_genes),
  pool = rownames(scRNA_harmony[[figure21_cfg$assay]]),
  nbin = figure21_cfg$module_nbin, ctrl = figure21_cfg$module_ctrl,
  assay = figure21_cfg$assay, slot = "data", name = "Figure21_Cytotoxicity",
  seed = figure21_cfg$module_seed, search = FALSE
)
figure21_scored_meta <- figure21_scored[[]]
if (!all(figure21_cells$cell %in% rownames(figure21_scored_meta))) {
  stop("模块评分与细胞条形码不匹配。")
}
figure21_cd8_cells <- figure21_cells %>% filter(anno == "CD8+ T cells")
figure21_cd8_cells$module_score <- figure21_scored_meta[
  figure21_cd8_cells$cell, "Figure21_Cytotoxicity1"
]
if (any(!is.finite(figure21_cd8_cells$module_score))) {
  stop("CD8 模块评分包含非有限值。")
}
figure21_cd8 <- figure21_cd8_cells %>%
  group_by(orig.ident, gender, group.a) %>%
  summarise(CD8_cyto = mean(module_score), n_CD8 = n(), .groups = "drop")
rm(figure21_scored, figure21_scored_meta, figure21_cd8_cells)

# 保留原始缺失排除逻辑：未捕获相应细胞的患者不补成零。
figure21_roc_input <- figure21_patients %>%
  mutate(gender = as.character(gender), group.a = as.character(group.a)) %>%
  left_join(figure21_mono, by = "orig.ident") %>%
  left_join(figure21_cd8, by = c("orig.ident", "gender", "group.a")) %>%
  mutate(response = as.integer(group.a == "OR"))
if (anyDuplicated(figure21_roc_input$orig.ident)) stop("ROC 合并后出现重复患者。")
figure21_predictors <- c("ESR2_mean", "mono_n", "CD8_cyto")
if (!all(vapply(figure21_roc_input[figure21_predictors], is.numeric, logical(1)))) {
  stop("ROC 三个输入指标必须为数值。")
}
figure21_roc_input$included <- complete.cases(
  figure21_roc_input[c("orig.ident", "gender", "response", figure21_predictors)]
) & Reduce(`&`, lapply(figure21_roc_input[figure21_predictors], is.finite))
figure21_roc_input$exclusion_reason <- apply(
  figure21_roc_input[figure21_predictors], 1, function(x) {
    paste(figure21_predictors[!is.finite(x)], collapse = "; ")
  }
)
message("ROC 纳入患者：", sum(figure21_roc_input$included),
        "；缺失或非有限指标排除：", sum(!figure21_roc_input$included))

figure21_z <- function(x, label) {
  s <- stats::sd(x)
  if (length(x) < 2L || any(!is.finite(x)) || !is.finite(s) || s == 0) {
    stop("无法计算 ", label, " 的性别内 Z-score。")
  }
  as.numeric(scale(x))
}
figure21_make_roc <- function(sex) {
  d <- figure21_roc_input %>% filter(included, gender == sex)
  if (!all(c(0L, 1L) %in% d$response)) stop(sex, " 缺少 OR 或 SOR 患者。")
  d$ESR2 <- d$ESR2_mean
  d$Microenv <- figure21_z(d$mono_n, "mono_n") - figure21_z(d$CD8_cyto, "CD8_cyto")
  d$Combined <- figure21_z(d$ESR2_mean, "ESR2_mean") + d$Microenv
  models <- c("ESR2", "Microenv", "Combined")
  rocs <- setNames(lapply(models, function(model) {
    pROC::roc(d$response, d[[model]], levels = c(0, 1),
              direction = figure21_cfg$roc_direction, percent = FALSE, quiet = TRUE)
  }), models)
  auc <- data.frame(
    gender = sex, model = models, n = nrow(d),
    n_OR = sum(d$response == 1L), n_SOR = sum(d$response == 0L),
    AUC = vapply(rocs, function(r) as.numeric(pROC::auc(r)), numeric(1)),
    direction = vapply(rocs, function(r) r$direction, character(1)), row.names = NULL
  )
  coordinates <- bind_rows(lapply(models, function(model) {
    r <- rocs[[model]]
    data.frame(gender = sex, model = model, threshold = r$thresholds,
               FPR = 1 - r$specificities, sensitivity = r$sensitivities)
  }))
  # 按原图顺序画蓝、橙、红；重叠时 Combined 红线最后绘制。
  coordinates$model <- factor(coordinates$model, levels = models)
  colors <- c(ESR2 = "#4575B4", Microenv = "#FDAE61", Combined = "#D73027")
  legend_labels <- setNames(sprintf("%s (AUC = %.2f)", models, auc$AUC), models)
  p <- ggplot(coordinates, aes(FPR, sensitivity, colour = model, group = model)) +
    geom_abline(slope = 1, intercept = 0, colour = '#B8BEC6', linetype = 2,
                  linewidth = 0.35) +
    geom_path(linewidth = 0.65) +
    scale_colour_manual(values = colors, breaks = models, labels = legend_labels) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       labels = scales::label_number(accuracy = 0.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       labels = scales::label_number(accuracy = 0.1), expand = c(0, 0)) +
    coord_cartesian() + labs(title = sex, x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
    guides(colour = guide_legend(ncol = 1)) + figure21_theme() +
    theme(panel.border = element_rect(colour = '#454545', fill = NA, linewidth = 0.35),
          axis.line = element_blank(),
          legend.key.width = grid::unit(15, 'pt'),
          legend.key.height = grid::unit(10, 'pt'),
          legend.text = element_text(size = figure21_cfg$font_size - 2),
          legend.spacing.y = grid::unit(0, 'pt'))
  list(plot = p, rocs = rocs, auc = auc, scores = d, coordinates = coordinates)
}
figure21_roc <- setNames(lapply(c("Female", "Male"), figure21_make_roc), c("Female", "Male"))
p2H <- figure21_roc$Female$plot
p2I <- figure21_roc$Male$plot
figure21_auc <- bind_rows(lapply(figure21_roc, `[[`, "auc"))
figure21_scores <- bind_rows(lapply(figure21_roc, `[[`, "scores"))
figure21_coordinates <- bind_rows(lapply(figure21_roc, `[[`, "coordinates"))

# 参考图的 AUC 仅用于审计差异，不参与评分、ROC 或图例的生成。
figure21_reference_auc <- data.frame(
  gender = rep(c("Female", "Male"), each = 3),
  model = rep(c("ESR2", "Microenv", "Combined"), 2),
  reference_AUC = c(0.75, 0.42, 0.42, 0.60, 0.80, 0.70)
)
figure21_auc_check <- figure21_auc %>%
  left_join(figure21_reference_auc, by = c("gender", "model")) %>%
  mutate(matches_reference_2dp = sprintf("%.2f", AUC) == sprintf("%.2f", reference_AUC))
if (any(!figure21_auc_check$matches_reference_2dp)) {
  warning("部分 ROC AUC 与参考图不同；已保留真实重算结果，请查看 AUC 审计表。", call. = FALSE)
}


# 09. 三行排版：字母独立留白、A/B 共享色条、底部坐标区对齐 --------------------
figure21_tag <- function(p, letter) {
  p + labs(tag = letter) +
    theme(plot.tag = element_text(size = 16, face = 'bold', colour = '#151515'),
          plot.tag.position = 'topleft', plot.tag.location = 'margin')
}
figure21_AB <- (
  figure21_tag(p2A, 'A') + figure21_tag(p2B_HD, 'B') + p2B_MM +
    plot_layout(ncol = 3, widths = c(1, 1, 1), guides = 'collect')
) & theme(legend.position = 'bottom')
# 只包裹整行/区块，避免跨行热图长标签把其他面板的左边距挤大。
figure21_top <- wrap_elements(full = patchworkGrob(figure21_AB)) +
  wrap_elements(full = ggplotGrob(figure21_tag(p2C, 'C'))) +
  plot_layout(widths = c(1.1, 1))
figure21_middle <- figure21_tag(p2D, 'D') + figure21_tag(p2E, 'E') +
  plot_layout(widths = c(0.95, 1.35))
figure21_bottom <- figure21_tag(p2F, 'F') + figure21_tag(p2G, 'G') +
  figure21_tag(p2H, 'H') + figure21_tag(p2I, 'I') +
  plot_layout(widths = c(0.9, 1.7, 1.4, 1.4))
p2_all <- wrap_elements(full = patchworkGrob(figure21_top)) /
  wrap_elements(full = patchworkGrob(figure21_middle)) /
  wrap_elements(full = patchworkGrob(figure21_bottom)) +
  plot_layout(heights = c(3.3, 6.2, 3.9)) +
  plot_annotation(theme = theme(plot.margin = margin(8, 8, 8, 8)))

# 10. 统一导出：独立面板、组合图、统计结果、数值源数据和运行设置 ---------------
dir.create(figure21_cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(figure21_cfg$out_dir)) stop('Cannot create output directory.')
figure21_panels <- list(A = p2A, B = p2B, C = p2C, D = p2D, E = p2E,
                        F = p2F, G = p2G, H = p2H, I = p2I)
figure21_sizes <- list(A = c(3.6, 4), B = c(6.4, 4), C = c(8.5, 4),
                       D = c(5.6, 6.8), E = c(8.2, 6.8), F = c(2.8, 3.8),
                       G = c(5.2, 3.8), H = c(4, 4.2), I = c(4, 4.2))
for (panel in names(figure21_panels)) {
  message('Saving panel ', panel)
  figure21_save(figure21_panels[[panel]], paste0('Figure21_', panel),
                figure21_sizes[[panel]][1], figure21_sizes[[panel]][2])
}
figure21_save(p2_all, 'Figure21_combined', figure21_cfg$figure_width,
                figure21_cfg$figure_height)

figure21_write_csv <- function(x, name) {
  write.csv(x, file.path(figure21_cfg$out_dir, paste0(name, '.csv')), row.names = FALSE)
}
figure21_write_csv(figure21_D_data, 'Figure21_D_source_data')
figure21_write_csv(figure21_E_stats, 'Figure21_E_cell_level_Wilcoxon_tests')
figure21_write_csv(bind_rows(mutate(figure21_F_stats, comparison = 'All'),
                             mutate(figure21_G_stats, comparison = gender)),
                   'Figure21_FG_cell_level_Welch_tests')
figure21_write_csv(figure21_patients, 'Figure21_patient_ESR2')
figure21_write_csv(figure21_roc_input, 'Figure21_HI_patient_features_and_inclusion')
figure21_write_csv(figure21_auc, 'Figure21_HI_AUC')
figure21_write_csv(figure21_auc_check, 'Figure21_HI_AUC_reference_check')
figure21_write_csv(figure21_scores, 'Figure21_HI_patient_scores')
figure21_write_csv(figure21_coordinates, 'Figure21_HI_ROC_coordinates')
figure21_write_csv(
  figure21_cells %>% count(group, gender, group.a, anno, name = 'n_cells'),
  'Figure21_cell_counts')
saveRDS(lapply(figure21_roc, function(x) x[c('rocs', 'auc', 'scores', 'coordinates')]),
        file.path(figure21_cfg$out_dir, 'Figure21_HI_ROC_objects.rds'))
saveRDS(figure21_cfg, file.path(figure21_cfg$out_dir, 'Figure21_settings.rds'))
writeLines(c(
  paste('Input:', normalizePath(figure21_cfg$input)),
  'Expression: RNA/data; ESR2 == 0 cells are retained.',
  'D: Z-score of group mean expression within cell type; detection rate = mean(ESR2 > 0).',
  'E: unpaired two-sided cell-level Wilcoxon test, exact=FALSE, correct=TRUE.',
  'E reference significance symbols are audited in the test CSV; plotted labels use actual P values.',
  'F/G: two-sided cell-level Welch t-tests; exploratory, cells are not independent patients.',
  'orig.ident is treated as the patient identifier.',
  paste('mono_n: captured cell count of', paste(figure21_cfg$mono_types, collapse = '; ')),
  paste('CD8 module genes:', paste(figure21_cfg$cyto_genes, collapse = ', ')),
  paste('AddModuleScore on all cells: RNA/data, all RNA genes as pool, seed=',
        figure21_cfg$module_seed, ', nbin=', figure21_cfg$module_nbin,
        ', ctrl=', figure21_cfg$module_ctrl, '.', sep = ''),
  'CD8_cyto: patient mean of module scores across CD8+ T cells.',
  'All three ROC curves use the same complete cases within each sex.',
  'ESR2 = ESR2_mean; Microenv = z(mono_n)-z(CD8_cyto); Combined adds z(ESR2_mean).',
  paste('ROC: OR=1, SOR=0; within-sex Z-score; direction =', figure21_cfg$roc_direction),
  'Hand-built scores assessed on the same cohort; not LASSO or cross-validation.',
  'Reference AUCs are used only in the comparison CSV, never as plotted values.',
  'UMAP color limits and D Z-score limits affect display only; no cell filtering.',
  'Output layout follows the reference panel organization, not a pixel-identical copy.'
), file.path(figure21_cfg$out_dir, 'analysis_settings.txt'))
capture.output(sessionInfo(), file = file.path(figure21_cfg$out_dir, 'sessionInfo.txt'))
print(figure21_auc_check)
message('Done. Panels and source data saved in: ', figure21_cfg$out_dir)
if (interactive()) print(p2_all)

