
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Load residual -----------------------------------------------------------

residual <- readxl::read_excel(path = "data/metadata/Residual.xlsx") %>% 
  dplyr::filter(Residual != "NA") %>% 
  dplyr::mutate(residual = ifelse(Residual == "R0", "R0", "non-R0")) %>% 
  dplyr::select(barcode, residual)

os_risk_group <- readxl::read_excel(path = "data/newoutput/Riskgroup-os.xlsx") %>% dplyr::left_join(residual, by = "barcode")
pfs_risk_group <- readxl::read_excel(path = "data/newoutput/Riskgroup-pfs.xlsx") %>% dplyr::left_join(residual, by = "barcode")



# Function ----------------------------------------------------------------


# Multi-Cox ---------------------------------------------------------------

# OS ----------------------------------------------------------------------
os_risk_group %>% 
  dplyr::select(barcode, oc, figo_stage, stage, CA125, age, histotype, platelet_count, riskscore, group, residual, event, duration) ->
  os_risk_group_s


# PFS ---------------------------------------------------------------------

pfs_risk_group %>% 
  dplyr::select(barcode, oc, figo_stage, stage, CA125, age, histotype, platelet_count, riskscore, group, residual, event, duration) ->
  pfs_risk_group_s


