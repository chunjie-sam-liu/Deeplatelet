# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Jul 27 08:02:00 2021
# @DESCRIPTION: residual.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)


new_residual <- readxl::read_excel(path = "data/metadata/血小板-样本信息汇总0729不含临床信息.xlsx")


new_residual %>% 
  dplyr::select(barcode, res = 2, figo_stage, stage) %>% 
  dplyr::mutate(figo_stage = ifelse(figo_stage == "NA", NA, figo_stage)) %>% 
  dplyr::mutate(stage = ifelse(stage == "NA", NA, stage)) %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(residual = ifelse(res == "R0", "R0", "non-R0")) %>% 
  dplyr::select(-res) ->
  residual

writexl::write_xlsx(x = residual, path = "data/metadata/residual.xlsx")
