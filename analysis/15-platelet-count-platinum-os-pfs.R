
# Library -----------------------------------------------------------------


library(magrittr)
library(ggplot2)
library(survival)

# Load data ---------------------------------------------------------------

total351.platinum.se <- readr::read_rds(file = 'data/rda/total351.platinum.se.rds.gz')

metadata <- as.data.frame(total351.platinum.se@colData) %>% tibble::as_tibble()



# Function ----------------------------------------------------------------


# Analysis ----------------------------------------------------------------


# Platinum platelet count -------------------------------------------------


metadata %>% 
  dplyr::mutate(platelet_count = as.numeric(platelet_count)) %>% 
  dplyr::filter(!is.na(platelet_count)) %>%  
  ggpubr::ggboxplot(
    x = 'platinum',
    y = 'platelet_count',
    color = 'platinum',
    xlab = 'Platinum',
    ylab = 'Platelet count (in 10^9/L)',
    title = "Platinum sensitivity correlation with platelet count"
  ) +
  ggpubr::stat_compare_means(label.x = 1.1) +
  scale_color_manual(name = "Platinum", values = c("#B22222", "#00008B"), labels = c("Sensitive", "Resistant")) +
  scale_x_discrete(labels = c("Sensitive", "Resistant")) +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill = NULL, colour = "black", size = 1),
    axis.line = element_blank(),
    legend.position = c(0.5, 0.8),
    legend.direction = "horizontal",
    
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
  ) ->
  plot_platinum_platelet_count;plot_platinum_platelet_count
ggsave(
  filename = 'platinum-platelet-count-correlation.pdf',
  plot = plot_platinum_platelet_count,
  device = 'pdf',
  path = 'data/newoutput/',
  width = 3,
  height = 4
)


# platinum age ------------------------------------------------------------


metadata %>% 
  ggpubr::ggboxplot(
    x = 'platinum',
    y = 'age',
    color = 'platinum',
    xlab = 'Platinum',
    ylab = 'Age',
    title = 'Platinum sensitivity correlation with age'
  ) +
  ggpubr::stat_compare_means(label.x = 1.1) +
  scale_color_manual(name = "Platinum", values = c("#B22222", "#00008B"), labels = c("Sensitive", "Resistant")) +
  scale_x_discrete(labels = c("Sensitive", "Resistant")) +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill = NULL, colour = "black", size = 1),
    axis.line = element_blank(),
    legend.position = "None",
    
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
  )  ->
  plot_platinum_age;plot_platinum_age

ggsave(
  filename = 'platinum-age-correlation.pdf',
  plot = plot_platinum_age,
  device = 'pdf',
  path = 'data/newoutput/',
  width = 3,
  height = 4
)



# Platinum ca125 ----------------------------------------------------------
metadata$CA125 <- scale(log2(metadata$CA125))[, 1]
metadata %>% 
  dplyr::mutate(CA125 = CA125) %>% 
  ggpubr::ggboxplot(
    x = 'platinum',
    y = 'CA125',
    color = 'platinum',
    xlab = 'Platinum',
    ylab = 'Normalized CA125 level',
    title = 'Platinum sensitivity and correlation with CA125'
  ) +
  ggpubr::stat_compare_means(label.x = 1.1) +
  scale_color_manual(name = "Platinum", values = c("#B22222", "#00008B"), labels = c("Sensitive", "Resistant")) +
  scale_x_discrete(labels = c("Sensitive", "Resistant")) +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill = NULL, colour = "black", size = 1),
    axis.line = element_blank(),
    legend.position = "None",
    
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
  ) ->
  plot_platinum_ca125;plot_platinum_ca125
ggsave(
  filename = 'platinum-ca125-correlation.pdf',
  plot = plot_platinum_ca125,
  device = 'pdf',
  path = 'data/newoutput/',
  width = 3,
  height = 4
)

# Platelet count os  -----------------------------------------------
total416.os.se <- readr::read_rds(file = 'data/rda/total416.os.se.norm.rds.gz')

metadata.os <- as.data.frame(total416.os.se@colData) %>% tibble::as_tibble() %>% 
  dplyr::filter(duration > 0) %>%
  dplyr::mutate(duration = ifelse(duration > 100, 100, duration)) %>%
  dplyr::mutate(platelet_count = as.numeric(platelet_count)) %>% 
  dplyr::filter(!is.na(platelet_count)) %>% 
  dplyr::mutate(PLT = as.factor(ifelse(platelet_count > 350, 'PLT>350', 'PLT<=350'))) %>% 
  dplyr::mutate(ca125_group = as.factor(ifelse(CA125 > 1200, 'CA125>1200', 'CA125<1200'))) %>% 
  dplyr::mutate(age_group = as.factor(ifelse(age > 50, 'age>50', 'age<=50')))
  

coxph(formula = Surv(time = duration, event = event) ~ platelet_count, data = metadata.os)

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ PLT, data = metadata.os),
  data = metadata.os,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1'),
  break.time.by = 20,
  ggtheme = theme_bw(),
  
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  
  ncensor.plot = TRUE,
  surv.median.line = 'hv',
  
  legend = 'top',
  legend.title = 'Group',
  xlab = 'Time in months',
  ylab = 'Overall survival probability',
  title = 'Overall survival with platelet count'
) ->
  os_platelet_count

ggsave(
  filename ='os-platelet-count.pdf',
  plot = print(os_platelet_count, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)


coxph(formula = Surv(time = duration, event = event) ~ CA125, data = metadata.os)

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ ca125_group, data = metadata.os),
  data = metadata.os,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1'),
  break.time.by = 20,
  ggtheme = theme_bw(),
  
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  
  ncensor.plot = TRUE,
  surv.median.line = 'hv',
  
  legend = 'top',
  legend.title = 'Group',
  xlab = 'Time in months',
  ylab = 'Overall survival probability',
  title = 'Overall survival with CA125 level'
) ->
  os_ca125

ggsave(
  filename ='os-ca125.pdf',
  plot = print(os_ca125, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)


coxph(formula = Surv(time = duration, event = event) ~ age, data = metadata.os)

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ age_group, data = metadata.os),
  data = metadata.os,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1'),
  break.time.by = 20,
  ggtheme = theme_bw(),
  
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  
  ncensor.plot = TRUE,
  surv.median.line = 'hv',
  
  legend = 'top',
  legend.title = 'Group',
  xlab = 'Time in months',
  ylab = 'Overall survival probability',
  title = 'Overall survival with age'
) ->
  os_age

ggsave(
  filename ='os-age.pdf',
  plot = print(os_age, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)

# Platelet count pfs ------------------------------------------------------

total434.pfs.se <- readr::read_rds(file = 'data/rda/total434.pfs.se.norm.rds.gz')
metadata.pfs <- as.data.frame(total434.pfs.se@colData) %>% tibble::as_tibble() %>% 
  dplyr::filter(duration > 0) %>%
  dplyr::mutate(duration = ifelse(duration > 60, 60, duration)) %>%
  plyr::mutate(platelet_count = as.numeric(platelet_count)) %>% 
  dplyr::filter(!is.na(platelet_count)) %>% 
  dplyr::mutate(PLT = as.factor(ifelse(platelet_count > 350, 'PLT>350', 'PLT<=350'))) %>% 
  dplyr::mutate(ca125_group = as.factor(ifelse(CA125 > 1200, 'CA125>1200', 'CA125<=1200'))) %>% 
  dplyr::mutate(age_group = as.factor(ifelse(age > 50, 'age>50', 'age<=50')))

coxph(formula = Surv(time = duration, event = event) ~ platelet_count, data = metadata.pfs) 

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ PLT, data = metadata.pfs),
  data = metadata.pfs,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1'),
  break.time.by = 20,
  ggtheme = theme_bw(),
  
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  
  ncensor.plot = TRUE,
  surv.median.line = 'hv',
  
  legend = 'top',
  legend.title = 'Group',
  xlab = 'Time in months',
  ylab = 'Progression free survival probability',
  title = 'Progression free survival with platelet count'
) ->
  pfs_platelet_count

ggsave(
  filename ='pfs-platelet-count.pdf',
  plot = print(pfs_platelet_count, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)

coxph(formula = Surv(time = duration, event = event) ~ CA125, data = metadata.pfs) 

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ ca125_group, data = metadata.pfs),
  data = metadata.pfs,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1'),
  break.time.by = 20,
  ggtheme = theme_bw(),
  
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  
  ncensor.plot = TRUE,
  surv.median.line = 'hv',
  
  legend = 'top',
  legend.title = 'Group',
  xlab = 'Time in months',
  ylab = 'Progression free survival probability',
  title = 'Progression free survival with CA125 level'
) ->
  pfs_ca125

ggsave(
  filename ='pfs-ca125.pdf',
  plot = print(pfs_ca125, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)


coxph(formula = Surv(time = duration, event = event) ~ age, data = metadata.pfs)

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ age_group, data = metadata.pfs),
  data = metadata.pfs,
  pval = TRUE,
  pval.method = TRUE,
  palette = RColorBrewer::brewer.pal(n=4, name = 'Set1'),
  break.time.by = 20,
  ggtheme = theme_bw(),
  
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 6,
  
  ncensor.plot = TRUE,
  surv.median.line = 'hv',
  
  legend = 'top',
  legend.title = 'Group',
  xlab = 'Time in months',
  ylab = 'Overall survival probability',
  title = 'Overall survival with age'
) ->
  pfs_age

ggsave(
  filename ='pfs-age.pdf',
  plot = print(pfs_age, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)

# Save image --------------------------------------------------------------

save.image(file = 'data/rda/10-platelet-count-platinum-os-pfs.rda')
load(file = 'data/rda/10-platelet-count-platinum-os-pfs.rda')