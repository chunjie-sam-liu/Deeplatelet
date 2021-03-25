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
# 
# total416.os.se <- readr::read_rds(file = 'data/rda/total416.os.se.duration.norm.rds.gz')
# total434.pfs.se <- readr::read_rds(file = 'data/rda/total434.pfs.se.duration.norm.rds.gz')

total416.os.expr.coxph.hazard_ratio <- readr::read_rds(file='data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz')
total434.pfs.expr.coxph.hazard_ratio <- readr::read_rds(file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz')

# Function ----------------------------------------------------------------

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

# Transform to feather ----------------------------------------------------
total416.os.expr.coxph.hazard_ratio.df <- total416.os.se[total416.os.expr.coxph.hazard_ratio$ensid,] %>% fn_se2df()
total434.pfs.expr.coxph.hazard_ratio.df <- total434.pfs.se[total434.pfs.expr.coxph.hazard_ratio$ensid,] %>% fn_se2df()

feather::write_feather(x = total416.os.expr.coxph.hazard_ratio.df, path = 'data/rda/total416.os.se.norm.coxph.feather')
feather::write_feather(x = total434.pfs.expr.coxph.hazard_ratio.df, path = 'data/rda/total434.pfs.se.norm.coxph.feather')

# test features
total434.pfs.se[total416.os.expr.coxph.hazard_ratio$ensid,] %>% 
  fn_se2df() %>% 
  feather::write_feather(path = 'data/rda/total434.pfs.se.norm.coxph.test.feather')

# Save image --------------------------------------------------------------

save.image(file = 'data/rda/10-univariate-features-deep-learning.rda')
