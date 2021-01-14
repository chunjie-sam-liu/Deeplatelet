
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

fn_expr_duration_event <- function(.se) {
  .meta <- .se@colData %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    dplyr::select(barcode, oc, duration, event)
  
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
  .coxph <- tryCatch(
    expr = survival::coxph(survival::Surv(time = duration, event = event) ~ expr, data = .x),
    error = function(e) {1},
    warning = function(w) {1}
  ) %>% summary()
  
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
# Analysis ----------------------------------------------------------------

total416.os.expr <- fn_expr_duration_event(.se = total416.os.se)
total434.pfs.expr <- fn_expr_duration_event(.se = total434.pfs.se)

cluster <- multidplyr::new_cluster(n = 30) 
multidplyr::cluster_library(cluster, 'magrittr')
multidplyr::cluster_assign(cluster, fn_cox_model = fn_cox_model, fn_coxr_each_gene = fn_coxr_each_gene)

total416.os.expr %>% 
  multidplyr::partition(cluster) %>% 
  dplyr::mutate(coxph = purrr::map(.x = data, .f = fn_coxr_each_gene)) %>% 
  dplyr::collect() ->
  total416.os.expr.coxph

readr::write_rds(x = total416.os.expr.coxph, file = 'data/rda/total416-os.expr.coxph.rds.gz', compress = 'gz')

total434.pfs.expr %>% 
  multidplyr::partition(cluster) %>% 
  dplyr::mutate(coxph = purrr::map(.x = data, .f = fn_coxr_each_gene)) %>% 
  dplyr::collect() ->
  total434.pfs.expr.coxph

readr::write_rds(x = total434.pfs.expr.coxph, file = 'data/rda/total434-pfs.expr.coxph.rds.gz', compress = 'gz')


# Save image --------------------------------------------------------------

save.image(file = 'data/rda/0072-univariate-cox-regression.rda')
