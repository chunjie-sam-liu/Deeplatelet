
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)

# Load data ---------------------------------------------------------------

metadata_platinum <- readr::read_rds(file = 'data/rda/metadata-platinum.rds.gz') %>% 
  dplyr::select(barcode, platinum = platinum_sensitivity)
total800.se <- readr::read_rds(file = 'data/rda/total800.se.rds.gz')
# src ---------------------------------------------------------------------

source(file = 'src/doparallel.R', local = TRUE)


# Merge data --------------------------------------------------------------

total800_meta <- total800.se@colData %>% 
  as.data.frame() %>% 
  tibble::as_tibble()

total800_meta %>% 
  dplyr::inner_join(metadata_platinum, by = 'barcode') %>% 
  dplyr::mutate(platinum = ordered(x = platinum, levels = c('sensitive', 'resistant'))) %>% 
  as.data.frame()->
  total315_meta_platinum

rownames(total315_meta_platinum) <- total315_meta_platinum$barcode

total315.platinum.se <- SummarizedExperiment::SummarizedExperiment(assays = assay(total800.se[,total315_meta_platinum$barcode]), colData = total315_meta_platinum)

readr::write_rds(x = total315.platinum.se, file = 'data/rda/total315.platinum.se.rds.gz', compress = 'gz')

# Save image --------------------------------------------------------------

save.image(file = 'data/rda/02-platinum-data-merge.rda')
