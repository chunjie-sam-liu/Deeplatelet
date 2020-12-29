
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

total351.platinum.se <- readr::read_rds(file = 'data/rda/total351.platinum.se.rds.gz')


# Function ----------------------------------------------------------------


# For filter
fn_filter_samples <- function(.se) {
  .m <- apply(X = assay(.se), MARGIN = 2, FUN = function(x) {sum(x) >= 8 * 2e6})
  .se[, .m]
}

# Normalize ---------------------------------------------------------------

fn_parallel_start(n_cores = 50)
# filter sample with total library size.
total351.platinum.se.fs <- fn_filter_samples(.se = total351.platinum.se)


# Save image --------------------------------------------------------------

save.image(file = 'data/rda/03-platinum-normalize.rda')
