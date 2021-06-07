
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(survival)

# Load data ---------------------------------------------------------------

metadata <- readxl::read_xlsx(path = "data/metadata/platelet-prognosis.xlsx", sheet = 1, skip = 1, col_names = F) %>%
  dplyr::select(1, 2, 3, 9, 10, 11, 12, 14, 17)

names(metadata) <- c("barcode", "hospital", "age", "pfs", "os", "pfs_status", "os_status", "platelet", "platelet_rate")


showtext::showtext_auto()

# PFS stat ----------------------------------------------------------------

metadata %>% 
  dplyr::mutate(hospital = factor(hospital, levels = c("河南省肿瘤医院","山东省肿瘤医院","同济医院","中山大学附属第一医院")))

metadata %>%
  dplyr::group_by(hospital, pfs_status) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = hospital, y = n, fill = pfs_status, label = n)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position =  position_dodge(width = 1), vjust = 0, hjust = 0.2, size = 6) +
  scale_fill_viridis_d(direction = -1, name = "复发") +
  scale_x_discrete(limits = (metadata$hospital %>% unique())[c(4, 2,1,3)]) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, color = 'black', size = 1),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.position = 'bottom',
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  labs(y = "Count") ->
  pfs_status_plot;pfs_status_plot

ggsave(
  filename = "PFS-status.pdf",
  plot = pfs_status_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 10,
  height = 7
)


# OS stat -----------------------------------------------------------------


metadata %>%
  dplyr::group_by(hospital, os_status) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = hospital, y = n, fill = os_status, label = n)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position =  position_dodge(width = 1), vjust = 0, hjust = 0.2, size = 6) +
  scale_fill_viridis_d(direction = -1, name = "死亡") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, color = 'black', size = 1),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.position = 'bottom',
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  labs(y = "Count") ->
  os_status_plot
ggsave(
  filename = "OS-status.pdf",
  plot = os_status_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 10,
  height = 7
)



# Platelet and rate ------------------------------------------------------

metadata %>%
  dplyr::mutate(
    pfs_status = as.numeric(plyr::revalue(x = pfs_status, replace = c("是" = 1, "否" = 0, "死亡" = 3, "进展" = 1))),
    platelet = as.numeric(platelet),
    platelet_rate = as.numeric(platelet_rate),
    os_status = as.numeric(plyr::revalue(x = os_status, replace = c("生存" = 0, "死亡" = 1)))
  ) %>%
  dplyr::mutate(pfs = pfs / 30, os = os / 30) ->
  metadata_clean

metadata_clean %>% 
  dplyr::mutate(hospital = plyr::revalue(x = hospital, replace = c(
    "中山大学附属第一医院" = "The First Affiliated Hospital of Sun Yat-sen University",
    "同济医院" = "Tongji Hospital",
    "山东省肿瘤医院" = "Shandong Tumour Hospital",
    "河南省肿瘤医院" = "Henan Tumour Hospital"
  ))) %>% 
  dplyr::mutate(pfs_status = as.character(pfs_status)) %>% 
  dplyr::mutate(pfs_status = plyr::revalue(x = pfs_status, replace = c("1" = "Yes", "0" = "No", "3" = "Dead"))) %>% 
  dplyr::mutate(os_status = as.character(os_status)) %>% 
  dplyr::mutate(os_status = plyr::revalue(x = os_status, replace = c("1" = "Dead", "0" = "Alive"))) %>% 
  dplyr::rename(Barcode = barcode, Hospital = hospital, Age = age, PFS = pfs, OS = os, `PFS status` = pfs_status, `OS status` = os_status, `Platelet count` = platelet, `Platelet rate` = platelet_rate) %>% 
  writexl::write_xlsx(path = "data/newoutput/biobank-platelet-count.xlsx")

# PFS platelet ------------------------------------------------------------

metadata_clean %>%
  dplyr::filter(pfs > 0) %>%
  dplyr::mutate(pfs = ifelse(pfs > 120, 120, pfs)) %>%
  dplyr::filter(pfs_status != 3)  ->
  metadata_clean_pfs

metadata_clean_pfs %>%
  dplyr::filter(!is.na(platelet)) %>%
  dplyr::mutate(PLT = as.factor(ifelse(platelet > 350, '>350', '<=350'))) %>%
  dplyr::mutate(age_group = as.factor(ifelse(age > 50, '>50', '<=50'))) ->
  metadata_clean_pfs_platelet

coxph(Surv(time = pfs, event = pfs_status) ~ PLT, data = metadata_clean_pfs_platelet) %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(
    term,
    hazard_ratio = estimate,
    low = conf.low,
    high = conf.high,
    pval = p.value
  )

survminer::ggsurvplot(
  fit = survfit(Surv(time = pfs, event = pfs_status) ~ PLT, data = metadata_clean_pfs_platelet),
  data = metadata_clean_pfs_platelet,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1')[c(2, 1)],
  break.time.by = 20,
  ggtheme = theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black', size=1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, color = 'black'),

    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
  ),

  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  risk.table.height = 0.3,
  risk.table.title = "",

  ncensor.plot = FALSE,
  surv.median.line = 'hv',

  legend = c(0.85, 0.85),
  legend.title = "Platelet count",
  legend.labs = c("<=350", ">350"),
  xlab = 'PFS (months)',
  ylab = 'Survival probability (PFS)'
) ->
  pfs_platelet_plot;pfs_platelet_plot
ggsave(
  filename ='PFS-platelet-count.pdf',
  plot = print(pfs_platelet_plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/newoutput',
  width = 6,
  height = 6
)

# PFS age -----------------------------------------------------------------

coxph(Surv(time = pfs, event = pfs_status) ~ age_group, data = metadata_clean_pfs_platelet) %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(
    term,
    hazard_ratio = estimate,
    low = conf.low,
    high = conf.high,
    pval = p.value
  )

survminer::ggsurvplot(
  fit = survfit(Surv(time = pfs, event = pfs_status) ~ age_group, data = metadata_clean_pfs_platelet),
  data = metadata_clean_pfs_platelet,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1')[c(2, 1)],
  break.time.by = 20,
  ggtheme = theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black', size=1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, color = 'black'),

    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
  ),

  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  risk.table.height = 0.3,
  risk.table.title = "",

  ncensor.plot = FALSE,
  surv.median.line = 'hv',

  legend = c(0.85, 0.85),
  legend.title = "Age",
  legend.labs = c("<=50", ">50"),
  xlab = 'PFS (months)',
  ylab = 'Survival probability (PFS)'
) ->
  pfs_age_plot;pfs_age_plot

ggsave(
  filename ='PFS-age.pdf',
  plot = print(pfs_age_plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/newoutput',
  width = 6,
  height = 6
)

# PFS platelet rate -------------------------------------------------------


metadata_clean_pfs %>%
  dplyr::filter(!is.na(platelet_rate)) %>%
  dplyr::mutate(PLT_rate = as.factor(ifelse(platelet_rate > 26, 'PLT_rate>25', 'PLT_rate<=25'))) ->
  metadata_clean_pfs_platelet_rate

coxph(Surv(time = pfs, event = pfs_status) ~ PLT_rate, data = metadata_clean_pfs_platelet_rate) %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(
    term,
    hazard_ratio = estimate,
    low = conf.low,
    high = conf.high,
    pval = p.value
  )

survminer::ggsurvplot(
  fit = survfit(Surv(time = pfs, event = pfs_status) ~ PLT_rate, data = metadata_clean_pfs_platelet_rate),
  data = metadata_clean_pfs_platelet_rate,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1')[c(2, 1)],
  break.time.by = 20,
  ggtheme = theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black', size=1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, color = 'black'),
    
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
  ),

  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  risk.table.height = 0.3,
  risk.table.title = "",
  
  ncensor.plot = FALSE,
  surv.median.line = 'hv',
  
  legend = c(0.85, 0.85),
  legend.title = "P-LCR",
  legend.labs = c("<=26", ">26"),
  xlab = 'PFS (months)',
  ylab = 'Survival probability (PFS)'
) ->
  pfs_platelet_rate_plot;pfs_platelet_rate_plot
ggsave(
  filename ='PFS-platelet-rate.pdf',
  plot = print(pfs_platelet_rate_plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/newoutput',
  width = 6,
  height = 6
)


# OS platelet -------------------------------------------------------------

metadata_clean %>%
  dplyr::filter(os > 0) %>%
  dplyr::mutate(os = ifelse(os > 120, 120, os)) ->
  metadata_clean_os


metadata_clean_os %>%
  dplyr::filter(!is.na(platelet)) %>%
  dplyr::mutate(PLT = as.factor(ifelse(platelet > 350, 'PLT>350', 'PLT<=350'))) %>%
  dplyr::mutate(age_group = as.factor(ifelse(age > 50, 'age>50', 'age<=50'))) ->
  metadata_clean_os_platelet

coxph(Surv(time = os, event = os_status) ~ PLT, data = metadata_clean_os_platelet) %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(
    term,
    hazard_ratio = estimate,
    low = conf.low,
    high = conf.high,
    pval = p.value
  )

survminer::ggsurvplot(
  fit = survfit(Surv(time = os, event = os_status) ~ PLT, data = metadata_clean_os_platelet),
  data = metadata_clean_os_platelet,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1')[c(2, 1)],
  break.time.by = 20,
  ggtheme = theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black', size=1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, color = 'black'),
    
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
  ),

  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  risk.table.height = 0.3,
  risk.table.title = "",
  
  ncensor.plot = FALSE,
  surv.median.line = 'hv',
  
  legend = c(0.85, 0.85),
  legend.title = "Platelet count",
  legend.labs = c("<=350", ">350"),
  xlab = 'OS (months)',
  ylab = 'Survival probability (OS)',
) ->
  os_platelet_plot;os_platelet_plot

ggsave(
  filename ='OS-platelet-count.pdf',
  plot = print(os_platelet_plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/newoutput',
  width = 6,
  height = 6
)

# OS age -----------------------------------------------------------------

coxph(Surv(time = os, event = os_status) ~ age_group, data = metadata_clean_os_platelet) %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(
    term,
    hazard_ratio = estimate,
    low = conf.low,
    high = conf.high,
    pval = p.value
  )

coxph(Surv(time = os, event = os_status) ~ age_group, data = metadata_clean_os_platelet) %>% 
  summary() -> 
  .coxph

.coxph$conf.int %>% 
  as.data.frame() %>% 
  tibble::as_tibble() %>% 
  dplyr::select(hazard_ratio = 1, ci_lower95 = 3, ci_upper95 = 4) %>% 
  dplyr::mutate(coxp = .coxph$waldtest[3])

survminer::ggsurvplot(
  fit = survfit(Surv(time = os, event = os_status) ~ age_group, data = metadata_clean_os_platelet),
  data = metadata_clean_os_platelet,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1')[c(2, 1)],
  break.time.by = 20,
  ggtheme = theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black', size=1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, color = 'black'),
    
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
  ),

  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  risk.table.height = 0.3,
  risk.table.title = "",
  
  ncensor.plot = FALSE,
  surv.median.line = 'hv',
  
  legend = c(0.85, 0.85),
  legend.title = "Age",
  legend.labs = c("<=50", ">50"),
  xlab = 'OS (months)',
  ylab = 'Survival probability (OS)',
) ->
  os_age_plot;os_age_plot

ggsave(
  filename ='OS-age.pdf',
  plot = print(os_age_plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/newoutput',
  width = 6,
  height = 6
)

# OS platelet rate -------------------------------------------------------


metadata_clean_os %>%
  dplyr::filter(!is.na(platelet_rate)) %>%
  dplyr::mutate(PLT_rate = as.factor(ifelse(platelet_rate > 26, 'PLT_rate>28', 'PLT_rate<=28'))) ->
  metadata_clean_os_platelet_rate

coxph(Surv(time = os, event = os_status) ~ PLT_rate, data = metadata_clean_os_platelet_rate) %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(
    term,
    hazard_ratio = estimate,
    low = conf.low,
    high = conf.high,
    pval = p.value
  )


survminer::ggsurvplot(
  fit = survfit(Surv(time = os, event = os_status) ~ PLT_rate, data = metadata_clean_os_platelet_rate),
  data = metadata_clean_os_platelet_rate,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1')[c(2, 1)],
  break.time.by = 20,
  ggtheme = theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black', size=1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, color = 'black'),
    
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
  ),

  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  risk.table.height = 0.3,
  risk.table.title = "",
  
  ncensor.plot = FALSE,
  surv.median.line = 'hv',
  
  legend = c(0.85, 0.85),
  legend.title = "P-LCR",
  legend.labs = c("<=26", ">26"),
  xlab = 'OS (months)',
  ylab = 'Survival probability (OS)',
) ->
  os_platelet_rate_plot;os_platelet_rate_plot
ggsave(
  filename ='OS-platelet-rate.pdf',
  plot = print(os_platelet_rate_plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/newoutput',
  width = 6,
  height = 6
)


# Save --------------------------------------------------------------------

save.image(file = "data/rda/00-plateletcount.rda")
load(file = "data/rda/00-plateletcount.rda")
