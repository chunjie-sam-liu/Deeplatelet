
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)

# Load data ---------------------------------------------------------------

metadata_os <- readr::read_rds(file = 'data/rda/metadata-os.rds.gz') %>% 
  dplyr::select(barcode, event = status, duration = os)
metadata_pfs <- readr::read_rds(file = 'data/rda/metadata-pfs.rds.gz') %>% 
  dplyr::select(barcode, event = palindromia, duration = pfs)
total800.se <- readr::read_rds(file = 'data/rda/total800.se.rds.gz')
# src ---------------------------------------------------------------------

source(file = 'src/doparallel.R', local = TRUE)

# Merge data --------------------------------------------------------------

total800_meta <- total800.se@colData %>% 
  as.data.frame() %>% 
  tibble::as_tibble()

total800_meta %>% 
  dplyr::inner_join(metadata_os, by = 'barcode') %>% 
  as.data.frame()->
  total416_meta_os
rownames(total416_meta_os) <- total416_meta_os$barcode

total416.os.se <- SummarizedExperiment::SummarizedExperiment(assays = assay(total800.se[,total416_meta_os$barcode]), colData = total416_meta_os)
readr::write_rds(x = total416.os.se, file = 'data/rda/total416.os.se.rds.gz', compress = 'gz')


# Pfs ---------------------------------------------------------------------


total800_meta %>% 
  dplyr::inner_join(metadata_pfs, by = 'barcode') %>% 
  as.data.frame() ->
  total434_meta_pfs
rownames(total434_meta_pfs) <- total434_meta_pfs$barcode

total434.pfs.se <- SummarizedExperiment::SummarizedExperiment(assays = assay(total800.se[,total434_meta_pfs$barcode]), colData = total434_meta_pfs)
readr::write_rds(x = total434.pfs.se, file = 'data/rda/total434.pfs.se.rds.gz', compress = 'gz')

# Save image --------------------------------------------------------------

save.image(file = 'data/rda/06-prepare-pfs-os.rda')
