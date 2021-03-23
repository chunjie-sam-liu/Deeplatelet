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
  .svobj <- sva(dat = assay(.se), mod = .mod, mod0 = .mod0, n.sv = 50)

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
  
  confounding$confounding <- setdiff(1:ncol(.svobj$sv), c(confounding$platinum, confounding$na))

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
