
# Library ----------------------------------------------------------------

library(magrittr)
library(DESeq2)

# Load data ---------------------------------------------------------------


platinum.se <- readr::read_rds(file = "data/rda/total351.platinum.se.norm.rds.gz")
os.se <- readr::read_rds(file = "data/rda/total416.os.se.norm.rds.gz")
pfs.se <- readr::read_rds(file = "data/rda/total434.pfs.se.norm.rds.gz")


list(
  Platinum = platinum.se@colData %>% as.data.frame(),
  OS = os.se@colData %>% as.data.frame(),
  PFS = pfs.se@colData %>% as.data.frame()
) %>% 
  writexl::write_xlsx(path = "data/rda/platinum-os-pfs-used-data.xlsx")


# Save image --------------------------------------------------------------

save.image(file = "data/rda/10-export.rda")
