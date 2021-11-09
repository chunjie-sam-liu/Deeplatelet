# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Nov  9 19:41:35 2021
# @DESCRIPTION: 18-enrichment.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# src ---------------------------------------------------------------------

source(file = "src/convertID.R")

# Load panel --------------------------------------------------------------

os.panel <- readr::read_rds(file='data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)
pfs.panel <- readr::read_rds(file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)


# Convert ID --------------------------------------------------------------

os.panel.df <- fn_convertId(ids = os.panel)
pfs.panel.df <- fn_convertId(ids = pfs.panel)

os.pfs.panel.df <- list(
  "OS gene panel" = os.panel.df,
  "PFS gene panel" = pfs.panel.df
)

writexl::write_xlsx(
  x = os.pfs.panel.df, 
  path = "data/newoutput/OS-PFS-Panel-Description.xlsx"
  )

