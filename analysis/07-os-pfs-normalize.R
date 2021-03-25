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

source(file = 'src/doparallel.R', local = TRUE)

# Load data ---------------------------------------------------------------

total416.os.se <- readr::read_rds(file = 'data/rda/total416.os.se.rds.gz')
total434.pfs.se <- readr::read_rds(file = 'data/rda/total434.pfs.se.rds.gz')

# Function ----------------------------------------------------------------

fn_filter_samples <- function(.se) {
  
  .lgl <- .se$X__mapped_reads > 5e6 & .se$mapping_rate >= 0.3
  .se[, .lgl]
}

fn_filter_genes_deseq_normalize <- function(.se, minCounts = 10, fSample = 0.9, thres = 3) {
  
  .keep_genes <- apply(X = assay(.se), MARGIN = 1, FUN = function(x) {
    sum(x >= minCounts) / length(x) > fSample
  })
  .se_filter_genes <- .se[.keep_genes, ]
  
  assay(.se_filter_genes) %>% 
    apply(MARGIN = 1, FUN = ineq::ineq, type = 'var') %>% 
    tibble::enframe() %>% 
    dplyr::filter(value < thres) %>% 
    dplyr::pull(name) -> 
    .keep_genes
  .se_filter_genes_ineq <- .se_filter_genes[.keep_genes, ]
  
  # DESeq vst transformation
  .vst <- varianceStabilizingTransformation(object = assay(.se_filter_genes_ineq), fitType = 'local')
  
  # Return SE
  SummarizedExperiment(assays = as.matrix(.vst), colData = as.data.frame(colData(.se_filter_genes_ineq)))
}
fn_filter_samples_by_cor <- function(.se, coef = 0.4, type = 'spearman') {
  # Filter the samples with lower correlation.
  assert_that(!missing(.se) && is(.se, "SummarizedExperiment"), msg = '.se must be provided with SummarizedExperiment.')
  
  .se_mat <- assay(.se)
  .total_number_sample <- ncol(.se_mat)
  
  .cor <- bigcor(x = .se_mat, fun = 'cor', size = floor(ncol(.se_mat) / 2), method = type, verbose = FALSE)
  .corr <- abs(.cor[1:nrow(.cor), 1:ncol(.cor)])
  .rem_samp <- colnames(.se_mat)[rowMeans(.corr) > coef]
  assert_that(length(.rem_samp) > 0, msg = 'One group has 0 sample remain.')
  message(glue::glue('Notice: {length(.rem_samp)} samples remained in group.'))
  
  .remain_samples <- .rem_samp
  
  message(glue::glue('Notice: {.total_number_sample - length(.remain_samples)} filtered. {length(.remain_samples)} samples remained.'))
  .se[, .remain_samples]
}

fn_remove_unwanted_variables <- function(.se, .vars = c('oc', 'age', 'lib.size', 'library'), .th_pval = 0.05, .th_var_strong = 0.1) {
  
  .vars <- if (! 'event' %in% .vars) c(.vars, 'event') else .vars
  # create model and null model
  .mod <-  model.matrix(~event, data = .se@colData)
  .mod0 <-  model.matrix(~1, data = .se@colData)
  .svobj <- sva(dat = assay(.se), mod = .mod, mod0 = .mod0, n.sv = 100)
  
  # confounding factors with surrogate variables
  # remove the factors that were identified as potential confounding variables from the dataset
  .matrix_w_corr_vars <- foreach(i = 1:ncol(.svobj$sv), .combine = rbind, .packages = c('magrittr')) %dopar% {
    .vars %>% purrr::map_dbl(.f = function(.x) {
      if (.x %in% c('event', 'oc', 'library')) {
        aov(formula = .svobj$sv[, i] ~ .se@colData[, .x]) %>%
          broom::tidy() %>% dplyr::slice(1) %>% dplyr::pull(p.value)
      } else {
        cor.test(
          x = .svobj$sv[, i],
          y = .se@colData[, .x],
          method = 'pearson'
        ) %>% broom::tidy() %>% dplyr::slice(1) %>% dplyr::pull(p.value)
      }
    }) -> .vars_pval
    names(.vars_pval) <- .vars
    
  # .strong_var <- names(sort(x = .vars_pval))[1]
    if (.vars_pval['event'] > .th_pval) {
      .strong_var <- names(sort(x = .vars_pval))[1]
    } else {
      .strong_var <- 'event'
    }
    
    if (.vars_pval[.strong_var] < .th_var_strong) {
      .verdict <- names(.vars_pval[.strong_var])
    } else {
      .verdict <- NA
    }
    
    c(.vars_pval, verdict = .verdict)
  }
  rownames(.matrix_w_corr_vars) <- 1:ncol(.svobj$sv)
  # remove the factors that were identified as potential confounding variables from the dataset
  .vars %>% purrr::map(.f = function(.x) {
    unname(which(.matrix_w_corr_vars[, 'verdict'] == .x))
  }) -> confounding
  names(confounding) <- .vars
  confounding$na <- unname(which(is.na(.matrix_w_corr_vars[, 'verdict'])))
  confounding$confounding <- setdiff(1:ncol(.svobj$sv), c(confounding$event))
  # confounding$confounding <- c(confounding$oc, confounding$age, confounding$lib.size, confounding$library)
  
  .data_rm_be <- removeBatchEffect(assay(.se), design = .mod, covariates = .svobj$sv[, confounding$confounding])
  SummarizedExperiment(assays = .data_rm_be, colData = .se@colData[colnames(.data_rm_be), ])
}

fn_se2df <- function(.se) {
  .expr <- assay(.se) %>% 
    t() %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var = 'barcode')
  .meta <- .se@colData %>%
    as.data.frame() %>%
    dplyr::select(barcode, oc, event, duration)
  
  .expr %>% 
    dplyr::left_join(.meta, by = 'barcode')
}

# Normalize ---------------------------------------------------------------
fn_parallel_start(n_cores = 50)

# os data normalization
total416.os.se.fs <- fn_filter_samples(.se = total416.os.se)
total416.os.se.norm <- fn_filter_genes_deseq_normalize(.se = total416.os.se.fs)
total416.os.se.norm.filter_sample <- fn_filter_samples_by_cor(.se = total416.os.se.norm)
total416.os.se.norm.filter_sample.rbe <- fn_remove_unwanted_variables(.se = total416.os.se.norm.filter_sample)
readr::write_rds(x = total416.os.se.norm.filter_sample.rbe, file = 'data/rda/total416.os.se.norm.rds.gz', compress = 'gz')

total416.os.se.norm.filter_sample.rbe.df <- fn_se2df(.se = total416.os.se.norm.filter_sample.rbe)
feather::write_feather(x = total416.os.se.norm.filter_sample.rbe.df, path = 'data/rda/total416.os.se.norm.feather')


# pfs data normalization
total434.pfs.se.fs <- fn_filter_samples(.se = total434.pfs.se)
total434.pfs.se.norm <- fn_filter_genes_deseq_normalize(.se = total434.pfs.se.fs)
total434.pfs.se.norm.filter_sample <- fn_filter_samples_by_cor(.se = total434.pfs.se.norm)
total434.pfs.se.norm.filter_sample.rbe <- fn_remove_unwanted_variables(.se = total434.pfs.se.norm.filter_sample)
readr::write_rds(x = total434.pfs.se.norm.filter_sample.rbe, file = 'data/rda/total434.pfs.se.norm.rds.gz', compress = 'gz')

total434.pfs.se.norm.filter_sample.rbe.df <- fn_se2df(.se = total434.pfs.se.norm.filter_sample.rbe)
feather::write_feather(x = total434.pfs.se.norm.filter_sample.rbe.df, path = 'data/rda/total434.pfs.se.norm.feather')


fn_parallel_stop()




# Save image --------------------------------------------------------------
save.image(file = 'data/rda/07-os-pfs-normalize.rda')

