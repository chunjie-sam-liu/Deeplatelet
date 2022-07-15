
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(survival)

# src ---------------------------------------------------------------------

source(file='src/doparallel.R', local = TRUE)

# Load data ---------------------------------------------------------------

total416.os.se <- readr::read_rds(file = 'data/rda/total416.os.se.norm.rds.gz')
total434.pfs.se <- readr::read_rds(file = 'data/rda/total434.pfs.se.norm.rds.gz')

# Function ----------------------------------------------------------------

fn_expr_duration_event <- function(.se, .type = 'os') {
  .cutoff <- ifelse(.type == 'os', 100, 60)
  .meta <- .se@colData %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    dplyr::select(barcode, oc, duration, event) %>% 
    dplyr::mutate(duration = ifelse(duration > .cutoff, .cutoff, duration))
  
  assay(.se) %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = 'barcode') %>% 
    tibble::as_tibble() %>% 
    dplyr::left_join(y = .meta, by = 'barcode') %>% 
    tidyr::gather(key = 'ensid', value = 'expr', -c('barcode', 'oc', 'duration', 'event')) %>% 
    dplyr::group_by(ensid) %>% 
    tidyr::nest() %>% 
    dplyr::ungroup()
}

fn_cox_model <- function(.x) {
  .cox <- tryCatch(
    expr = survival::coxph(survival::Surv(time = duration, event = event) ~ expr, data = .x),
    error = function(e) {NULL},
    warning = function(w) {NULL}
  ) 
  
  if(is.null(.cox)) {
    return(tibble::tibble(hazard_ratio = NA, ci_lower95 = NA, ci_upper95 = NA, coxp = NA))
  }
  
  .coxph <- summary(.cox)
  
  .coxph$conf.int %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    dplyr::select(hazard_ratio = 1, ci_lower95 = 3, ci_upper95 = 4) %>% 
    dplyr::mutate(coxp = .coxph$waldtest[3])
}

fn_coxr_each_gene <- function(.x) {
  .x_all <- fn_cox_model(.x) %>% dplyr::mutate(name = 'all')
  .x_oc521 <- fn_cox_model(.x %>% dplyr::filter(oc == 'OC521')) %>% dplyr::mutate(name = 'OC521')
  .x_oc44 <- fn_cox_model(.x %>% dplyr::filter(oc == 'OC44')) %>% dplyr::mutate(name = 'OC44')
  .x_oc79 <- fn_cox_model(.x %>% dplyr::filter(oc == 'OC79')) %>% dplyr::mutate(name = 'OC79')
  .x_oc172 <- fn_cox_model(.x %>% dplyr::filter(oc == 'OC172')) %>% dplyr::mutate(name = 'OC172')
  
  dplyr::bind_rows(.x_all, .x_oc521, .x_oc44, .x_oc79, .x_oc172)
}

fn_lasso <- function(.se, .hr, .type = 'os') {
  
  as.data.frame(.se@colData) %>% 
    dplyr::filter(duration > 0) %>% 
    # dplyr::filter(oc == 'OC521') %>%
    dplyr::pull(barcode) -> 
    .samples
  
  .d <- .se[.hr$ensid, .samples]
  
  .x <- t(assay(.d))
  .y <- as.data.frame(.d@colData) %>% 
    dplyr::select(time = duration, status = event) %>% 
    as.matrix()
  
  
  .fit <- glmnet::glmnet(
    x = .x,
    y = .y,
    family = 'cox'
  )
  
  # pdf(file = glue::glue('data/output/{.type}-multicox-lasso.pdf'))
  # plot(.fit, main = glue::glue('{.type}-multicox-lasso'))
  # dev.off()
  
  .cvfit <- glmnet::cv.glmnet(x = .x, y = .y, family = 'cox')
  
  # pdf(file = glue::glue('data/output/{.type}-multicox-lasso-lambda.pdf'))
  # plot(.cvfit, main=glue::glue('{.type}-multicox-lasso-lambda'))
  # dev.off()
  
  .coef.min <- coef(.cvfit, s = 'lambda.min')
  
  as.matrix(.coef.min) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'ensid') %>%
    dplyr::rename(coef = `1`) %>%
    dplyr::filter(abs(coef) > 0)
  # as.matrix(.coef.min) %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column(var = 'ensid') %>%
  #   dplyr::rename(coef = `1`)
}
# Analysis ----------------------------------------------------------------

total416.os.expr <- fn_expr_duration_event(.se = total416.os.se, .type = 'os')
total434.pfs.expr <- fn_expr_duration_event(.se = total434.pfs.se, .type = 'pfs')

cluster <- multidplyr::new_cluster(n = 20) 

multidplyr::cluster_library(cluster, 'magrittr')
multidplyr::cluster_assign(cluster, fn_cox_model = fn_cox_model, fn_coxr_each_gene = fn_coxr_each_gene)

# OS ----------------------------------------------------------------------


# Univariate cox ----------------------------------------------------------


total416.os.expr %>% 
  multidplyr::partition(cluster) %>% 
  dplyr::mutate(coxph = purrr::map(.x = data, .f = fn_coxr_each_gene)) %>% 
  dplyr::collect() ->
  total416.os.expr.coxph

readr::write_rds(x = total416.os.expr.coxph, file = 'data/rda/total416-os.expr.coxph.rds.gz', compress = 'gz')

total416.os.expr.coxph %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest(coxph) %>% 
  dplyr::filter(coxp < 0.05) %>% 
  dplyr::filter(name == 'OC521') ->
  total416.os.expr.coxph.hazard_ratio



# Multicox lasso ----------------------------------------------------------

total416.se.multicox <- fn_lasso(.se = total416.os.se, .hr = total416.os.expr.coxph.hazard_ratio, .type = 'OS')

total416.os.expr.coxph.hazard_ratio %>% 
  dplyr::filter(ensid %in% total416.se.multicox$ensid) %>%
  readr::write_rds(file = 'data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz', compress = 'gz')

# PFS ---------------------------------------------------------------------


total434.pfs.expr %>% 
  multidplyr::partition(cluster) %>% 
  dplyr::mutate(coxph = purrr::map(.x = data, .f = fn_coxr_each_gene)) %>% 
  dplyr::collect() ->
  total434.pfs.expr.coxph

readr::write_rds(x = total434.pfs.expr.coxph, file = 'data/rda/total434-pfs.expr.coxph.rds.gz', compress = 'gz')

total434.pfs.expr.coxph %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest(coxph) %>% 
  dplyr::filter(coxp < 0.05) %>% 
  dplyr::filter(name == 'OC521') ->
  total434.pfs.expr.coxph.hazard_ratio

readr::write_rds(x = total434.pfs.expr.coxph.hazard_ratio, file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz', compress = 'gz')


# Multicox lasso ----------------------------------------------------------

total434.se.multicox <- fn_lasso(.se = total434.pfs.se, .hr =total434.pfs.expr.coxph.hazard_ratio, .type = 'pfs')

total434.pfs.expr.coxph.hazard_ratio %>% 
  dplyr::filter(ensid %in% total434.se.multicox$ensid) %>%
  readr::write_rds(file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz', compress = 'gz')


# Intersection ------------------------------------------------------------

ensid_inter <- intersect(total416.os.expr.coxph.hazard_ratio$ensid, total434.pfs.expr.coxph.hazard_ratio$ensid)

total416.os.expr.coxph.hazard_ratio %>% 
  dplyr::filter(ensid %in% ensid_inter) ->
  total416.os.expr.coxph.hazard_ratio.inter

total416.se.multicox_inter <- fn_lasso(
  .se = total416.os.se, 
  .hr = total416.os.expr.coxph.hazard_ratio.inter, 
  .type = 'OS'
)


total434.pfs.expr.coxph.hazard_ratio %>% 
  dplyr::filter(ensid %in% ensid_inter) ->
  total434.pfs.expr.coxph.hazard_ratio.inter

total434.se.multicox_inter <- fn_lasso(
  .se = total434.pfs.se, 
  .hr = total434.pfs.expr.coxph.hazard_ratio.inter, 
  .type = 'pfs'
)


union(total416.se.multicox_inter$ensid, total434.se.multicox_inter$ensid) -> final_gene_set
# new write
total416.os.expr.coxph.hazard_ratio %>% 
  dplyr::filter(ensid %in% final_gene_set) %>%
  readr::write_rds(file = 'data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz', compress = 'gz')
total434.pfs.expr.coxph.hazard_ratio %>% 
  dplyr::filter(ensid %in% final_gene_set) %>%
  readr::write_rds(file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz', compress = 'gz')

# a -----------------------------------------------------------------------



ggvenn::ggvenn(
  data = list(
    a = total416.se.multicox_inter$ensid,
    b = total434.se.multicox_inter$ensid
  )
) ->p1


ggvenn::ggvenn(
  data = list(
    c = total416.se.multicox$ensid,
    d = total434.se.multicox$ensid
  )
) ->p2

ggvenn::ggvenn(
  data = list(
    a = total416.se.multicox_inter$ensid,
    b = total434.se.multicox_inter$ensid,
    c = total416.se.multicox$ensid,
    d = total434.se.multicox$ensid
  )
) ->p3
library(patchwork)

p3 | p1/p2
#
# Save image --------------------------------------------------------------

save.image(file = 'data/rda/09-os-pfs-univariate-cox-regression.rda')
load('data/rda/09-os-pfs-univariate-cox-regression.rda')
