
# Library -----------------------------------------------------------------

library(magrittr)
library(DESeq2)

# Load data ---------------------------------------------------------------



metadata_os <- readr::read_rds(file = 'data/rda/metadata-os.rds.gz') %>% 
  dplyr::select(barcode, event = status, duration = os)
metadata_pfs <- readr::read_rds(file = 'data/rda/metadata-pfs.rds.gz') %>% 
  dplyr::select(barcode, event = palindromia, duration = pfs)
wuhan.se <- readr::read_rds(file = 'data/rda/wuhan.se.rds.gz')

filepath_updated_0510_os <- "data/metadata/Riskgroup-os-5-10(1).xlsx"
filepath_updated_0510_pfs <- "data/metadata/Riskgroup-pfs-5-10(1).xlsx"

metadata_os_0510 <- readxl::read_excel(path = filepath_updated_0510_os) %>% 
  dplyr::select(barcode, event, duration)
metadata_pfs_0510 <- readxl::read_excel(path = filepath_updated_0510_pfs) %>% 
  dplyr::select(barcode, event, duration)

# src ---------------------------------------------------------------------

source(file = 'src/doparallel.R', local = TRUE)

# Merge data --------------------------------------------------------------

wuhan_meta <- wuhan.se@colData %>% 
  as.data.frame() %>% 
  tibble::as_tibble()

wuhan_meta %>% 
  # dplyr::inner_join(metadata_os, by = 'barcode') %>% 
  dplyr::inner_join(metadata_os_0510, by = "barcode") %>% 
  as.data.frame()->
  total416_meta_os
rownames(total416_meta_os) <- total416_meta_os$barcode

total416.os.se <- SummarizedExperiment::SummarizedExperiment(assays = assay(wuhan.se[,total416_meta_os$barcode]), colData = total416_meta_os)
readr::write_rds(x = total416.os.se, file = 'data/rda/total416.os.se.rds.gz', compress = 'gz')


# Pfs ---------------------------------------------------------------------


wuhan_meta %>% 
  # dplyr::inner_join(metadata_pfs, by = 'barcode') %>%
  dplyr::inner_join(metadata_pfs_0510, by = "barcode") %>% 
  as.data.frame() ->
  total434_meta_pfs
rownames(total434_meta_pfs) <- total434_meta_pfs$barcode

total434.pfs.se <- SummarizedExperiment::SummarizedExperiment(assays = assay(wuhan.se[,total434_meta_pfs$barcode]), colData = total434_meta_pfs)
readr::write_rds(x = total434.pfs.se, file = 'data/rda/total434.pfs.se.rds.gz', compress = 'gz')

# Save image --------------------------------------------------------------

save.image(file = 'data/rda/06-prepare-pfs-os.rda')
