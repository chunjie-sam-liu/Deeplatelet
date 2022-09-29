# Library -----------------------------------------------------------------

library(magrittr)
library(assertthat)
library(DESeq2)
library(doParallel)
library(propagate)
library(sva)
library(limma)
library(ggplot2)

# src ---------------------------------------------------------------------

source(file = "src/doparallel.R", local = TRUE)

# Load data ---------------------------------------------------------------

total351.platinum.se <- readr::read_rds(file = "data/rda/total351.platinum.se.rds.gz")





# Function ----------------------------------------------------------------

fn_filter_samples <- function(.se) {
  # .m <- apply(X = assay(.se), MARGIN = 2, FUN = function(x) {
  #   sum(x) >= 5e6
  # })
  .lgl <- .se$X__mapped_reads > 5e6 & .se$mapping_rate >= 0.3
  .se[, .lgl]
}
fn_filter_genes_deseq_normalize <- function(.se, minCounts = 10, fSample = 0.9, thres = 3) {

  .keep_genes <- apply(X = assay(.se), MARGIN = 1, FUN = function(x) {
    sum(x >= minCounts) / length(x) > fSample
  })
  .se_filter_genes <- .se[.keep_genes, ]

  classes <- .se_filter_genes@colData$platinum
  levels(classes) %>%
    purrr::map(.f = function(.x) {
      assay(.se_filter_genes)[, classes == .x] %>%
        apply(MARGIN = 1, FUN = ineq::ineq, type = "var") %>%
        tibble::enframe()
    }) %>%
    purrr::reduce(.f = dplyr::left_join, by = "name") %>%
    dplyr::filter_if(.predicate = is.numeric, .vars_predicate = dplyr::all_vars(. < thres)) %>%
    dplyr::pull(name) ->
  .keep_genes
  .se_filter_genes_ineq <- .se_filter_genes[.keep_genes, ]

  # DESeq vst transformation
  .vst <- varianceStabilizingTransformation(object = assay(.se_filter_genes_ineq), fitType = "local")

  # Return SE
  SummarizedExperiment(assays = as.matrix(.vst), colData = as.data.frame(colData(.se_filter_genes_ineq)))
}
fn_filter_samples_by_cor <- function(.se, coef = 0.4, type = "spearman") {
  # Filter the samples with lower correlation.
  assert_that(!missing(.se) && is(.se, "SummarizedExperiment"), msg = ".se must be provided with SummarizedExperiment.")

  .se_mat <- assay(.se)
  .total_number_sample <- ncol(.se_mat)

  classes <- .se@colData$platinum
  levels(classes) %>%
    purrr::map(.f = function(.x) {
      .se_mat_x <- .se_mat[, classes == .x]
      .cor <- bigcor(x = .se_mat_x, fun = "cor", size = floor(ncol(.se_mat_x) / 2), method = type, verbose = FALSE)
      .corr <- abs(.cor[1:nrow(.cor), 1:ncol(.cor)])
      .rem_samp <- colnames(.se_mat_x)[rowMeans(.corr) > coef]
      assert_that(length(.rem_samp) > 0, msg = "One group has 0 sample remain.")
      message(glue::glue("Notice: {length(.rem_samp)} samples remained in group {.x}."))
      .rem_samp
    }) %>%
    purrr::reduce(.f = c) ->
  .remain_samples

  message(glue::glue("Notice: {.total_number_sample - length(.remain_samples)} filtered. {length(.remain_samples)} samples remained."))
  .se[, .remain_samples]
}

fn_remove_unwanted_variables <- function(.se, .vars = c("oc", "age", "lib.size", "library"), .th_pval = 0.05, .th_var_strong = 0.1) {
  .vars <- if (!"platinum" %in% .vars) c(.vars, "platinum") else .vars

  .mod <- model.matrix(~platinum, data = .se@colData)
  .mod0 <- model.matrix(~1, data = .se@colData)
  .svobj <- sva(dat = assay(.se), mod = .mod, mod0 = .mod0, n.sv = 100)

  fn_parallel_start(n_cores = 50)
  .matrix_w_corr_vars <- foreach(
    i = 1:ncol(.svobj$sv),
    .combine = rbind,
    .packages = c("magrittr")
  ) %dopar% {
    .vars %>%
      purrr::map_dbl(
        .f = function(.x) {
          if (.x %in% c("platinum", "oc", "library")) {
            aov(formula = .svobj$sv[, i] ~ .se@colData[, .x]) %>%
              broom::tidy() %>%
              dplyr::slice(1) %>%
              dplyr::pull(p.value)
          } else {
            cor.test(
              x = .svobj$sv[, i],
              y = .se@colData[, .x],
              method = "pearson"
            ) %>%
              broom::tidy() %>%
              dplyr::slice(1) %>%
              dplyr::pull(p.value)
          }
        }
      ) -> .vars_pval
    names(.vars_pval) <- .vars

    if (.vars_pval["platinum"] > .th_pval) {
      .strong_var <- names(sort(x = .vars_pval))[1]
    } else {
      .strong_var <- "platinum"
    }
    if (.vars_pval[.strong_var] < .th_var_strong) {
      .verdict <- names(.vars_pval[.strong_var])
    } else {
      .verdict <- NA
    }

    c(.vars_pval, verdict = .verdict)
  }
  fn_parallel_stop()
  rownames(.matrix_w_corr_vars) <- 1:ncol(.svobj$sv)

  .vars %>% 
    purrr::map(.f = function(.x) {
    unname(which(.matrix_w_corr_vars[, "verdict"] == .x))
  }) -> 
    confounding
  names(confounding) <- .vars
  confounding$na <- unname(which(is.na(.matrix_w_corr_vars[, "verdict"])))
  
  confounding$confounding <- setdiff(1:ncol(.svobj$sv), c(confounding$platinum))

  readr::write_rds(x = list(confounding = confounding,svobj = .svobj,mod = .mod, dbdat = assay(.se)), file = "data/rda/confounding-svobj.rds")
  
  .data_rm_be <- removeBatchEffect(assay(.se), design = .mod, covariates = .svobj$sv[, confounding$confounding])
  
  SummarizedExperiment(assays = .data_rm_be, colData = .se@colData[colnames(.data_rm_be), ])
}
# Normalize ---------------------------------------------------------------

total351.platinum.se.fs <- fn_filter_samples(.se = total351.platinum.se)
total351.platinum.se.norm <- fn_filter_genes_deseq_normalize(.se = total351.platinum.se.fs)
total351.platinum.se.norm.filter_sample <- fn_filter_samples_by_cor(.se = total351.platinum.se.norm)

total351.platinum.se.norm.filter_sample.rbe <- fn_remove_unwanted_variables(.se = total351.platinum.se.norm.filter_sample)

readr::write_rds(x = total351.platinum.se.norm.filter_sample.rbe, file = "data/rda/total351.platinum.se.norm.rds.gz", compress = "gz")


# Save image --------------------------------------------------------------

save.image(file = "data/rda/03-platinum-normalize.rda")
load(file = "data/rda/03-platinum-normalize.rda")


total351.platinum.se.fs$group = factor(total351.platinum.se.fs$platinum, levels = c("resistant", "sensitive" ), ordered = FALSE)

total351.platinum.se.fs_res <- results(DESeq(DESeqDataSet(total351.platinum.se.fs[rowMeans(assay(total351.platinum.se.fs)) >= 5, ], design = ~ group)))


total351.platinum.se.fs_res %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "ensg") %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::filter(baseMean > 5) %>% 
  dplyr::mutate(log2FC = ifelse(log2FoldChange > 4, 4, log2FoldChange)) %>%
  dplyr::mutate(log2FC = ifelse(log2FC < -4, -4, log2FC)) %>%
  # dplyr::mutate(log2FC = log2FoldChange) %>% 
  dplyr::mutate(FDR = -log10(padj)) %>% 
  dplyr::mutate(FDR = ifelse(FDR > 8, 8, FDR)) %>% 
  dplyr::mutate(
    color = dplyr::case_when(
      log2FC > log2(1.3) & padj < 0.05 ~ "red",
      log2FC < log2(1/1.3) & padj < 0.05 ~ "green",
      TRUE ~ "grey"
    )
  ) ->
  de

de %>% 
  dplyr::filter(color != "grey") %>% 
  dplyr::group_by(color) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  tibble::deframe() ->
  de_x

# de %>% dplyr::group_by(color) %>% dplyr::count()
de_x

de %>% 
  ggplot(aes(x = log2FC, y = FDR, color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values =  c("#6BAB62", "grey", "#9A84B2")) +
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    color = "black",
    linetype = 88,
    data = tibble::tibble(
      x1 = c(-Inf, log2(1.3), log2(1/1.3), log2(1.3)),
      y1 = -log10(0.05),
      x2 = c(log2(1/1.3), Inf, log2(1/1.3), log2(1.3) ),
      y2 = c(-log10(0.05), -log10(0.05), Inf, Inf))
  ) +
  # ggrepel::geom_text_repel(
  #   aes(label = GeneName),
  #   data = subset(.xd, color == "red") %>% 
  #     dplyr::arrange(-FDR, -abs(log2FC)) %>% 
  #     dplyr::slice(1:10),
  #   box.padding = 0.5,
  #   max.overlaps = Inf,
  #   # size = 6
  # ) +
  # ggrepel::geom_text_repel(
  #   aes(label = GeneName),
  #   data = subset(.xd, color == "green") %>% 
  #     dplyr::arrange(-FDR, -abs(log2FC)) %>% 
  #     dplyr::slice(1:10),
  #   box.padding = 0.5,
  #   max.overlaps = Inf,
  #   # size = 6
  # ) +
  scale_x_continuous(
    # limits = c(-4, 4),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    limits = c(
      0,
      ceiling(
        max(
          de %>% 
            dplyr::filter(!is.infinite(FDR)) %>% 
            dplyr::pull(FDR)
        ) / 10
      ) * 10
    )
  ) +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    axis.line.x.bottom = element_line(color = "black"),
    axis.line.y.left = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    
    legend.position = "none",
    
  ) +
  labs(
    x = "log2FC",
    y = "-log10(FDR)",
    title = glue::glue("Up (n={de_x[2]}), Down (n={de_x[1]}); |log2FC| > 1.3, FDR < 0.05; Pt response")
  )  ->
  p

de %>% 
  dplyr::filter(color != "grey") %>% 
  writexl::write_xlsx("data/review/deg.xlsx")


ggsave(
  filename = "de.pdf",
  plot = p,
  device = "pdf",
  path = "data/review",
  width = 6,
  height = 5
)


os.panel <- readr::read_rds(file='data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)


intersect(os.panel, de %>% dplyr::filter(color != "grey") %>% dplyr::pull(ensg))

de %>% 
  dplyr::filter(color != "grey") %>% 
  dplyr::filter(ensg %in% os.panel) 