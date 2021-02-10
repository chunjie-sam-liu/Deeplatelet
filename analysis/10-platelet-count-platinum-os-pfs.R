
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
  ggpubr::stat_compare_means() +
  ggthemes::theme_few() ->
  plot_platinum_platelet_count
ggsave(
  filename = 'platinum-platelet-count-correlation.pdf',
  plot = plot_platinum_platelet_count,
  device = 'pdf',
  path = 'data/output/',
  width = 6.5,
  height = 4.2
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
  ggpubr::stat_compare_means() +
  ggthemes::theme_few() ->
  plot_platinum_age

ggsave(
  filename = 'platinum-age-correlation.pdf',
  plot = plot_platinum_age,
  device = 'pdf',
  path = 'data/output/',
  width = 6.5,
  height = 4.2
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
  ggpubr::stat_compare_means() +
  ggthemes::theme_few() ->
  plot_platinum_ca125
ggsave(
  filename = 'platinum-ca125-correlation.pdf',
  plot = plot_platinum_ca125,
  device = 'pdf',
  path = 'data/output/',
  width = 6.5,
  height = 4.2
)

# Platelet count os  -----------------------------------------------
total416.os.se <- readr::read_rds(file = 'data/rda/total416.os.se.norm.rds.gz')

metadata.os <- as.data.frame(total416.os.se@colData) %>% tibble::as_tibble() %>% 
  dplyr::filter(duration > 0) %>%
  dplyr::mutate(duration = ifelse(duration > 100, 100, duration)) %>%
  dplyr::mutate(platelet_count = as.numeric(platelet_count)) %>% 
  dplyr::filter(!is.na(platelet_count)) %>% 
  dplyr::mutate(PLT = as.factor(ifelse(platelet_count > 350, 'PLT>350', 'PLT<=350'))) %>% 
  dplyr::mutate(ca125_group = as.factor(ifelse(CA125 > 1200, 'CA125>1200', 'CA125<1200')))
  

coxph(formula = Surv(time = duration, event = event) ~ platelet_count, data = metadata.os)

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ PLT, data = metadata.os),
  data = metadata.os,
  pval = TRUE
)


coxph(formula = Surv(time = duration, event = event) ~ CA125, data = metadata.os)

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ ca125_group, data = metadata.os),
  data = metadata.os,
  pval = TRUE
)

# Platelet count pfs ------------------------------------------------------

total434.pfs.se <- readr::read_rds(file = 'data/rda/total434.pfs.se.norm.rds.gz')
metadata.pfs <- as.data.frame(total434.pfs.se@colData) %>% tibble::as_tibble() %>% 
  # dplyr::filter(duration > 0) %>% 
  # dplyr::mutate(duration = ifelse(duration > 60, 60, duration)) %>% 
  plyr::mutate(platelet_count = as.numeric(platelet_count)) %>% 
  dplyr::filter(!is.na(platelet_count)) %>% 
  dplyr::mutate(PLT = as.factor(ifelse(platelet_count > 400, 'PLT>400', 'PLT<=400'))) %>% 
  dplyr::mutate(ca125_group = as.factor(ifelse(CA125 > 1200, 'CA125>1200', 'CA125<1200')))

coxph(formula = Surv(time = duration, event = event) ~ platelet_count, data = metadata.pfs) 

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ PLT, data = metadata.pfs),
  data = metadata.pfs,
  pval = TRUE
)

coxph(formula = Surv(time = duration, event = event) ~ CA125, data = metadata.pfs) 

survminer::ggsurvplot(
  fit = survfit(Surv(time = duration, event = event) ~ ca125_group, data = metadata.pfs),
  data = metadata.pfs,
  pval = TRUE
)

