
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Path --------------------------------------------------------------------


filepath <- 'data/metadata/血小板预后2020年11月V1chijh.xlsx'


# Load data ---------------------------------------------------------------

metadata <- readxl::read_excel(path = filepath, sheet = 6, skip = 1, col_types = c('text', 'text', 'text', 'date', 'text', 'text', 'text', 'text', 'date', 'date', 'numeric', 'text', 'text', 'date', 'date', 'text', 'date', 'text', 'date', 'numeric'))



# Filter samples ----------------------------------------------------------


metadata %>% 
  dplyr::select(barcode, patient_ID, oc, time_on_diagnosis, end_of_chemotherapy, platinum_sensitivity, palindromia, palindromia2, time_of_palindromia, time_of_palindromia2, pfs, alive_type, alive_type2, time_of_die, time_of_followup, os) %>% 
  dplyr::select(barcode, patient_ID, oc, time_on_diagnosis, platinum_sensitivity, palindromia = palindromia2, pfs, alive_type, os) %>% 
  dplyr::filter(!is.na(time_on_diagnosis)) ->
  metadata

# save metadata
readr::write_rds(metadata, file = 'data/metadata/metadata.rds.gz', compress = 'gz')

# Platinum ----------------------------------------------------------------

metadata %>% 
  dplyr::select(barcode, patient_ID, oc, platinum_sensitivity) %>% 
  dplyr::filter(platinum_sensitivity != 'NA') %>% 
  dplyr::mutate(platinum_sensitivity = ifelse(platinum_sensitivity == 1, 'sensitive', 'resistant')) %>% 
  dplyr::mutate(platinum_sensitivity = factor(x = platinum_sensitivity, levels = c('sensitive', 'resistant'))) ->
  metadata_platinum

metadata_platinum %>% 
  dplyr::group_by(oc, platinum_sensitivity) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(oc = factor(x = oc, levels = c('OC521', 'OC44', 'OC79', 'OC172'))) %>% 
  ggplot(aes(x = oc, y = n, fill = platinum_sensitivity)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = n), vjust = -1, size = 6) +
  scale_fill_brewer(palette = 'Set1', name = 'Platinum') +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(limits = c('OC521', 'OC44', 'OC79', 'OC172'), labels = c('TC', 'DC', 'VC1', 'VC2')) +
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 1),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.position = 'bottom',
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  labs(
    y = 'Number of samples',
    title = 'Platinum sensitivity data distribution'
  ) ->
  metadata_platinum_plot
ggsave(
  filename = 'data/output/dist-platinum.pdf',
  plot = metadata_platinum_plot,
  device = 'pdf',
  width = 9,
  height = 8
)

# PFS ---------------------------------------------------------------------
metadata %>% 
  dplyr::select(barcode, patient_ID, oc, palindromia, pfs) %>% 
  dplyr::filter(!is.na(pfs)) %>% 
  dplyr::filter(palindromia != 'NA') %>% 
  dplyr::mutate(palindromia = ifelse(palindromia == '1', 1, 0)) ->
  metadata_pfs

# OS ----------------------------------------------------------------------
metadata %>% 
  dplyr::select(barcode, patient_ID, oc, status = alive_type, os) %>% 
  dplyr::filter(status != 'NA') %>% 
  dplyr::mutate(status = ifelse(status == 'alive', 0, 1)) ->
  metadata_os





