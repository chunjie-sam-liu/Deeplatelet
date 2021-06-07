

# Select features
fn_se2data.frame <- function(.se) {
  .df <- as.data.frame(t(assay(.se)))
  .target <- colData(x = .se)[, 'platinum']
  .target <- ordered(x = .target, levels = c('sensitive', 'resistant'))
  .df.t <- cbind(platinum = .target, .df)
  .df.t
}

fn_se2task <- function(.se, .id = 'task') {
  .d <- cbind(
    t(assay(.se)),
    data.frame(platinum = factor(x = as.character(colData(.se)[, 'platinum']), levels = c('sensitive', 'resistant')))
  )
  mlr::makeClassifTask(
    id = .id, data = .d,
    target = 'platinum', positive = 'sensitive'
  )
}

fn_se2task_ca125 <- function(.se, .id) {
  .d <- as.data.frame(.se@colData[, c('CA125', 'oc', 'platinum')])
  .d$platinum <- factor(x = as.character(.d$platinum), levels = c('sensitive', 'resistant'))
  .rownames <- rownames(.d)
  .d %>%
    tibble::as_tibble() %>%
    dplyr::mutate(oc = ifelse(oc == 'OC521', 'OC521', 'others')) %>%
    dplyr::group_by(oc) %>%
    dplyr::mutate(CA125 = scale(log2(CA125))[, 1]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-oc) %>%
    as.data.frame() ->
    .d
  rownames(.d) <- .rownames
  identical(rownames(.d), rownames(t(assay(.se))))
  
  mlr::makeClassifTask(
    id = .id, data = .d,
    target = 'platinum', positive = 'sensitive'
  )
}
fn_se2task_panel_ca125 <- function(.se, .id) {
  .d <- as.data.frame(.se@colData[, c('CA125', 'oc', 'platinum')])
  # .d <- as.data.frame(.se@colData[, c('CA125', 'platinum')])
  .d$platinum <- factor(x = as.character(.d$platinum), levels = c('sensitive', 'resistant'))
  .rownames <- rownames(.d)
  .d %>%
    tibble::as_tibble() %>%
    dplyr::mutate(oc = ifelse(oc == 'OC521', 'OC521', 'others')) %>%
    dplyr::group_by(oc) %>%
    dplyr::mutate(CA125 = scale(log2(CA125))[, 1]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-oc) %>%
    as.data.frame() ->
    .d
  
  identical(.rownames, rownames(t(assay(.se))))
  .d <- cbind(t(assay(.se)), .d)
  mlr::makeClassifTask(
    id = .id, data = .d,
    target = 'platinum', positive = 'sensitive'
  )
}

fn_sd_median <- function(.df) {
  .sd <- apply(
    X = .df[, -1],
    MARGIN = 2,
    FUN = sd
  )
  .median <- apply(
    X = .df[, -1],
    MARGIN = 2,
    FUN = median
  )
  intersect(
    names(.sd[.sd > 0]),
    names(.median[.median > 0])
  ) -> .feats.sd.median
  .feats.sd.median
}



sprintf_transformer <- function(code, envir) {
  m <- regexpr("%.+$", code)
  if (m != -1) {
    format <- regmatches(code, m)
    regmatches(code, m) <- ""
    res <- evaluate::evaluate(code, envir)
    do.call(sprintf, list(format, res))
  } else {
    evaluate::evaluate(code, envir)
  }
}

glue_fmt <- function(..., .envir = parent.frame()) {
  glue::glue(..., .transformer = sprintf_transformer, .envir = .envir)
}


fn_plot_ca125 <- function(.d, .t) {
  .d %>%
    ggplot(aes(x = class, y = CA125, color = class)) +
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_boxplot(outlier.colour = NA, width = 0.5) +
    geom_point(position = position_jitter(width = 0.05), alpha = 0.6, size = 1.3) +
    scale_x_discrete(name = 'Class', limits = c('M', 'B'), labels = c('Manignant', 'Benign')) +
    scale_color_manual(values = c("#e31a1c", "#1f78b4")) +
    labs(x = 'Class', y = 'Normalized CA125 level', title = glue::glue('{.t}, CA125 differences')) +
    ggthemes::theme_few() +
    theme(legend.position = 'none')
}


fn_ifs <- function(.x, .feats, .tose) {
  .feats <- head(.feats, .x)
  # 1. split data
  oc521.feats.se <- .tose[.feats, .tose$oc == 'OC521']
  oc44.feats.se <- .tose[.feats, .tose$oc == 'OC44']
  oc79.feats.se <- .tose[.feats, .tose$oc == 'OC79']
  oc172.feats.se <- .tose[.feats, .tose$oc == 'OC172']

  # 2. build model
  oc521.feats.model <- fn_caret_train_model(.se = oc521.feats.se)

  # 3. ensembl predict
  oc521.auc <- fn_ensembl_predict(.model = oc521.feats.model, .se = oc521.feats.se, .title = 'OC521')
  oc44.auc <- fn_ensembl_predict(.model = oc521.feats.model, .se = oc44.feats.se, .title = 'OC44')
  oc79.auc <- fn_ensembl_predict(.model = oc521.feats.model, .se = oc79.feats.se, .title = 'OC79')
  oc172.auc <- fn_ensembl_predict(.model = oc521.feats.model, .se = oc172.feats.se, .title = 'OC172')

  # 4. return
  tibble::tibble(
    num = .x,
    oc521 = oc521.auc$alist$auc,
    oc172 = oc172.auc$alist$auc,
    oc79 = oc79.auc$alist$auc,
    oc44 = oc44.auc$alist$auc,
    oc79_error = list(rownames(oc79.auc$predict)[oc79.auc$predict$correct == FALSE]),
    oc172_error = list(rownames(oc172.auc$predict)[oc172.auc$predict$correct == FALSE]),
    oc44_error = list(rownames(oc44.auc$predict)[oc44.auc$predict$correct == FALSE])
  )
}


fn_dge <- function(.se) {
  # transform se to dds
  .ddsse <- DESeqDataSet(se = .se, design = ~ class)
  .dds_deseq <- DESeq(object = .ddsse)
  # .res <- results(object = .dds_deseq)
  .resLFC <- lfcShrink(dds = .dds_deseq, coef = resultsNames(object = .dds_deseq)[2], type = 'apeglm')
  .resLFC
}
fn_filter_genes <- function(.se) {

  message(glue::glue('Total gene number is {nrow(.se)}.'))
  # unorder se
  .se$class <- factor(x = as.character(x = .se$class), levels = c('M', 'B'))
  .keep <- apply(X = assay(.se), MARGIN = 1, FUN = function(x) {sum(x > 10) / length(x) > 0.9})
  .se_filter <- .se[.keep, ]
  message(glue::glue('After filter with frac 10 in 90 percent, the gene number is {nrow(.se_filter)}.'))

  # Using inequality to filter hypervariants
  .classes <- .se_filter@colData$class
  levels(.classes) %>%
    purrr::map(.f = function(.x) {
      assay(.se_filter)[, .classes == .x] %>%
        apply(MARGIN = 1, FUN = ineq::ineq, type = 'var') %>%
        tibble::enframe()
    }) %>%
    purrr::reduce(.f = dplyr::left_join, by = 'name') %>%
    dplyr::filter_if(.predicate = is.numeric, .vars_predicate = dplyr::all_vars(. < 3)) %>%
    dplyr::pull(name) ->
    .keep_genes
  .se_filter_ineq <- .se_filter[.keep_genes, ]
  message(glue::glue('With inequality filter, the gene number is {nrow(.se_filter_ineq)}.'))
  .se_filter_ineq
}
fn_vst_normal <- function(.se) {
  # DESeq vst transformation
  .vst <- varianceStabilizingTransformation(object = assay(.se), fitType = 'local')
  .vst.se <- SummarizedExperiment(assays = as.matrix(.vst), colData = as.data.frame(colData(.se)))
  .vst.se
}
fn_dge_vst <- function(.se) {

  # filter genes by counts
  .se_filter <- fn_filter_genes(.se = .se)
  .se_dge <- fn_dge(.se = .se_filter)
  .se_vst <- fn_vst_normal(.se = .se_filter)

  # return
  list(
    dge = .se_dge,
    vst = .se_vst
  )
}

fn_cluster <- function(.vst) {

  # assay(x = .vst) %>% apply(1, scale) %>% t() -> .vst_mat
  .vst_mat <- assay(x = .vst) %>% apply(1, scale) %>% t()
  colnames(.vst_mat) <- colnames(x = .vst)

  # heatmap top annotation
  cluster_col <- c("#458B74", "#68228B")
  names(cluster_col) <- c('B','M')
  hma_top = HeatmapAnnotation(
    df = data.frame(Type = .vst$class),
    gap = unit(c(2,2), "mm"),
    col = list(Type = cluster_col)
  )

  .hm <- Heatmap(
    # data and color
    matrix = .vst_mat,
    col = colorRamp2(c(-1.1, 0, 1.1), c("blue", "white", "red"), space = "RGB"),
    name = "Normalized counts",
    na_col = 'grey', color_space = 'LAB', rect_gp = gpar(col = NA),
    border = NA, cell_fun = NULL, layer_fun = NULL,

    # title
    # row_title = 'Selected genes', # OC44
    row_title = '', # OC521
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 14),
    row_title_rot = 90,
    # column_title = 'Samples', # OC44
    column_title = '',
    column_title_side = 'top',
    column_title_gp = gpar(fontsize = 14),
    column_title_rot = 0,

    # clustering of row
    cluster_rows = T,
    cluster_row_slices = T,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D",
    row_dend_side = 'left',
    row_dend_width = unit(10, 'mm'),
    # show_row_dend = T, # OC44
    show_row_dend = F, # OC521
    row_dend_reorder = T,
    row_dend_gp = gpar(),

    # clustering of column
    cluster_columns = T,
    cluster_column_slices = T,
    clustering_distance_columns = "pearson",
    clustering_method_columns = "ward.D",
    column_dend_side = 'top',
    column_dend_height = unit(10, 'mm'),
    # show_column_dend = T, # OC44
    show_column_dend = F, # OC521
    column_dend_gp = gpar(),
    column_dend_reorder = TRUE,

    row_order = NULL,
    column_order = NULL,

    # row labels
    row_labels = rownames(.vst_mat),
    row_names_side = 'left',
    show_row_names = F,
    row_names_max_width = unit(6, 'cm'),
    row_names_gp = gpar(fontsize = 12),
    row_names_rot = 0,
    row_names_centered = FALSE,

    # column labels
    column_labels = colnames(.vst_mat),
    column_names_side = 'bottom',
    # show_column_names = T, # OC44
    show_column_names = F, # OC521
    column_names_max_height = unit(6, 'cm'),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 90,
    column_names_centered = FALSE,

    # annotation
    top_annotation = hma_top,
    bottom_annotation = NULL,
    left_annotation = NULL,
    right_annotation = NULL,

    # kmeans cluster number
    # row cluster is 1
    # column cluster is 2 with 10 repeats
    km = 1,
    split = NULL,
    row_km = 1,
    row_km_repeats = 1,
    row_split = NULL,
    column_km = 2,
    column_km_repeats = 10,
    column_split = NULL,
    gap = unit(1, 'mm'),
    row_gap = unit(1, 'mm'),
    column_gap = unit(1, 'mm'),

    show_heatmap_legend = T,
    heatmap_legend_param = list(title = 'Normalized counts'),

    # raster_device = 'tiff',
    raster_quality = 2,
    raster_device_param = list(),
    raster_resize = F,

    post_fun = NULL
  )

  .cluster <- ComplexHeatmap::column_order(object = .hm)
  list(
    heatmap = .hm,
    cluster = .cluster
  )
}
fn_test_cluster <- function(.class, .cluster) {
  .class <- as.character(.class)
  .class_order <- .class[c(.cluster$`1`, .cluster$`2`)]
  .cluster_name <- rep(x = c('cluster1', 'cluster2'),  c(length(.cluster$`1`), length(.cluster$`2`)))
  .class_cluster <- data.frame(class = .class_order, clsuter = .cluster_name)
  .table <- table(.class_cluster)
  .fisher_test <- fisher.test(.table)

  list(
    table = .table,
    fisher_test = .fisher_test
  )
}
fn_save_hm_pdf <- function(.hm, .filename, .path) {
  .filename <- file.path(.path, .filename)
  pdf(file = .filename, width = 7, height = 8)
  ComplexHeatmap::draw(object = .hm)
  dev.off()
}

scientific_10 <- function(x) {
  parse(text=gsub("e[+|-]", " %*% 10^", scales::scientific_format()(x)))
}


human_read <- function(.x){
  .sign = ifelse(.x < 0 , TRUE, FALSE)
  .x <- abs(.x)

  if (.x >= 0.1) {
    .x %>% signif(digits = 2) %>% toString() -> .xx
  } else if (.x < 0.1 && .x >= 0.001 ) {
    .x %>% signif(digits = 2) %>% toString() -> .xx
  } else if (.x < 0.001 && .x > 0) {
    .x %>% format(digits = 3, scientific = TRUE) -> .xx
  } else {
    .xx <- '0'
  }

  ifelse(.sign, paste0('-',.xx), .xx)
}

human_read_latex_pval <- function(.x, .s = NA) {

  if (is.na(.s)) {
    if (grepl(pattern = "e", x = .x)) {
      sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
      latex2exp::TeX(glue::glue("$\\textit{P}=<<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>"))
    } else {
      latex2exp::TeX(glue::glue("$\\textit{P}=<<.x>>$", .open = "<<", .close = ">>"))
    }
  } else {
    if (grepl(pattern = "e", x = .x)) {
      sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
      latex2exp::TeX(glue::glue("<<.s>>, $\\textit{P}=<<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>"))
    } else {
      latex2exp::TeX(glue::glue("<<.s>>, $\\textit{P}=<<.x>>$", .open = "<<", .close = ">>"))
    }
  }

}
