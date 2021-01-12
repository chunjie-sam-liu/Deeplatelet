
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
fn_get_duration_event <- function(.se) {
  .se@colData[,c('duration', 'event')]
}

# Analysis ----------------------------------------------------------------



