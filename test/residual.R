# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Jul 27 08:02:00 2021
# @DESCRIPTION: residual.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

residual <- readxl::read_excel(path = "data/metadata/Residual.xlsx") %>% 
  dplyr::filter(Residual != "NA") %>% 
  dplyr::mutate(residual = ifelse(Residual %in% c("R0", "R1"), "R0", "non-R0")) %>% 
  dplyr::select(barcode, residual)


new_residual <- readxl::read_excel(path = "data/metadata/血小板-样本信息汇总0721V7.xlsx")


new_residual %>% 
  dplyr::select(barcode, res = 32) %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(residual = ifelse(res == "R0", "R0", "non-R0")) %>% 
  dplyr::select(-res) ->
  residual

writexl::write_xlsx(x = residual, path = "data/metadata/residual.xlsx")
