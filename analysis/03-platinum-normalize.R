
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

total315.platinum.se <- readr::read_rds(file = 'data/rda/total315.platinum.se.rds.gz')


# Function ----------------------------------------------------------------


# Normalize ---------------------------------------------------------------





# Save image --------------------------------------------------------------

save.image(file = 'data/rda/03-platinum-normalize.rda')
