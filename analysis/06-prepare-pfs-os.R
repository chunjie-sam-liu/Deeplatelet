
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)

# Load data ---------------------------------------------------------------

metadata_os <- readr::read_rds(file = 'data/rda/metadata-os.rds.gz') %>% 
  dplyr::rename(event = status, duration = os)
metadata_pfs <- readr::read_rds(file = 'data/rda/metadata-pfs.rds.gz') %>% 
  dplyr::rename(event = palindromia, duration = pfs)
total800.se <- readr::read_rds(file = 'data/rda/total800.se.rds.gz')
