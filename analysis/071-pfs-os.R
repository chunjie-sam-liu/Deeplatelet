
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(survival)

# src ---------------------------------------------------------------------

source(file='src/doparallel.R', local = TRUE)

# Load data ---------------------------------------------------------------

total416.os.se <- readr::read_rds(file = 'data/rda/total416.os.se.rds.gz')
total434.pfs.se <- readr::read_rds(file = 'data/rda/total434.pfs.se.rds.gz')


# Function ----------------------------------------------------------------


# Analysis ----------------------------------------------------------------



