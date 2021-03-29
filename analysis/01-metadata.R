
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

# save rds
readr::write_rds(x = metadata_platinum, file = 'data/rda/metadata-platinum.rds.gz', compress = 'gz')

# plot
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
  dplyr::filter(pfs > 0) %>% 
  dplyr::mutate(palindromia = ifelse(palindromia == '1', 1, 0)) ->
  metadata_pfs

readr::write_rds(x = metadata_pfs, file = 'data/rda/metadata-pfs.rds.gz', compress = 'gz')


metadata_pfs %>% 
  dplyr::group_by(oc, palindromia) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(palindromia = factor(x = palindromia, levels = c(1, 0))) %>% 
  dplyr::mutate(oc = factor(x = oc, levels = c('OC521', 'OC44', 'OC79', 'OC172'))) %>% 
  ggplot(aes(x = oc, y = n, fill = palindromia)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = n), vjust = -1, size = 6) +
  scale_fill_brewer(palette = 'Set1', name = 'Palindromia') +
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
    title = 'Palindromia data distribution'
  ) ->
  metadata_pfs_dist_plot

ggsave(
  filename = 'data/output/dist-palindromia.pdf',
  plot = metadata_pfs_dist_plot,
  device = 'pdf',
  width = 9,
  height = 8
)

give.n <- function(x, .upper_limit = max(metadata_pfs$pfs)){
  .y = max(x) * c(1.2, 1.1)
  .label = c(glue::glue('count={length(x)}'), glue::glue('median={median(x)}'))
  print(.label)
  return(tibble::tibble(y = .y, label = .label)) 
}

metadata_pfs %>% 
  dplyr::mutate(palindromia = factor(x = palindromia, levels = c(1, 0))) %>% 
  dplyr::mutate(oc = factor(x = oc, levels = c('OC521', 'OC44', 'OC79', 'OC172'))) %>% 
  ggplot(aes(x = oc, y = pfs, color = palindromia)) +
  geom_boxplot() +
  scale_color_brewer(palette = 'Set1', name = 'Palindromia') +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(limits = c('OC521', 'OC44', 'OC79', 'OC172'), labels = c('TC', 'DC', 'VC1', 'VC2')) +
  stat_summary(
    fun.data = give.n,
    geom = 'text',
    hjust = 0.5,
    vjust = 0.9,
    position = position_dodge(width = 0.75),
    size = 5
  ) +
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 1),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = 'bottom',
    legend.key = element_rect(fill = NA)
  ) +
  labs(
    y = 'PFS (Month)',
    title = 'Progression free survival distribution'
  ) ->
  metadata_pfs_dist_pfs_plot

ggsave(
  filename = 'data/output/dist-pfs.pdf',
  plot = metadata_pfs_dist_pfs_plot,
  device = 'pdf',
  width = 9,
  height = 8
)
  

# OS ----------------------------------------------------------------------
metadata %>% 
  dplyr::select(barcode, patient_ID, oc, status = alive_type, os) %>% 
  dplyr::filter(status != 'NA') %>% 
  dplyr::filter(os > 0) %>% 
  dplyr::mutate(status = ifelse(status == 'alive', 0, 1)) ->
  metadata_os

readr::write_rds(x = metadata_os, file = 'data/rda/metadata-os.rds.gz', compress = 'gz')

give.n <- function(x, .upper_limit = max(metadata_os$os)){
  .y = max(x) * c(1.2, 1.1)
  .label = c(glue::glue('count={length(x)}'), glue::glue('median={median(x)}'))
  print(.label)
  return(tibble::tibble(y = .y, label = .label)) 
}

metadata_os %>% 
  dplyr::mutate(status = factor(x = status, levels = c(1, 0))) %>% 
  dplyr::mutate(oc = factor(x = oc, levels = c('OC521', 'OC44', 'OC79', 'OC172'))) %>% 
  ggplot(aes(x = oc, y = os, color = status)) +
  geom_boxplot() +
  scale_color_brewer(palette = 'Set1', name = 'Status') +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_discrete(limits = c('OC521', 'OC44', 'OC79', 'OC172'), labels = c('TC', 'DC', 'VC1', 'VC2')) +
  stat_summary(
    fun.data = give.n,
    geom = 'text',
    hjust = 0.5,
    vjust = 0.9,
    position = position_dodge(width = 0.5),
    size = 5
  ) +
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 1),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = 'bottom',
    legend.key = element_rect(fill = NA)
  ) +
  labs(
    y = 'OS (Month)',
    title = 'Overall survival distribution'
  ) ->
  metadata_os_dist_plot

ggsave(
  filename = 'data/output/dist-os.pdf',
  plot = metadata_os_dist_plot,
  device = 'pdf',
  width = 9,
  height = 8
)

# Save image --------------------------------------------------------------

save.image(file = 'data/rda/01-medata.rda')
load(file = 'data/rda/01-medata.rda')