
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(survival)

# src ---------------------------------------------------------------------

source(file='src/doparallel.R', local = TRUE)

# Load data ---------------------------------------------------------------

total416.os.se <- readr::read_rds(file = 'data/rda/total416.os.se.norm.rds.gz')
total434.pfs.se <- readr::read_rds(file = 'data/rda/total434.pfs.se.norm.rds.gz')


# Function ----------------------------------------------------------------
fn_get_duration_event <- function(.se) {
  .se@colData[,c('barcode', 'oc', 'duration', 'event')] %>% 
    as.data.frame() %>% 
    tibble::as_tibble()
}

fn_surv_plot <- function(.df, .ylab = 'Overall survival probability', .title = 'Overall survival data distribution') {
  .matchname <- c('OC521' = 'TC', 'OC44' = 'DC', 'OC79'='VC1', 'OC172' = 'CV2')
  .df %>%
    dplyr::mutate(group = plyr::revalue(oc, replace = .matchname)) %>% 
    dplyr::mutate(group = factor(x = group, levels = .matchname)) ->
    .df.group
  
  .fit <- survfit(Surv(time = duration, event = event) ~ group, data = .df.group)
  
  survminer::ggsurvplot(
    fit = .fit,
    data = .df.group,
    pval = TRUE,
    pval.mehtod = TRUE,
    palette = RColorBrewer::brewer.pal(n=4, name = 'Set1'),
    break.time.by = 20,
    ggtheme = theme_bw(),
    
    risk.table = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    risk.table.col = 'strata',
    risk.table.fontsize = 6,
    
    ncensor.plot = TRUE,
    surv.median.line = 'hv',
    
    legend = 'top',
    legend.labs = .matchname,
    legend.title = 'Group',
    xlab = 'Time in months',
    ylab = .ylab,
    title = .title
  )
}
# Analysis ----------------------------------------------------------------

# OS ----------------------------------------------------------------------


fn_get_duration_event(.se = total416.os.se) %>% 
  fn_surv_plot(.ylab = 'Overall survival probability', .title = 'Overall survival data distribution') -> os.plot

ggsave(
  filename ='os-plot.pdf',
  plot = print(os.plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)

# set cutoff for durations as 100 months.

total416.os.se.duration <- total416.os.se
total416.os.se.duration$duration <- ifelse(total416.os.se.duration$duration > 100, 100, total416.os.se.duration$duration)

fn_get_duration_event(.se = total416.os.se.duration) %>% 
  fn_surv_plot(.ylab = 'Overall survival probability', .title = 'Overall survival data distribution') -> os.duration.plot

ggsave(
  filename ='os-duration-plot.pdf',
  plot = print(os.duration.plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)


# PFS ---------------------------------------------------------------------



fn_get_duration_event(.se = total434.pfs.se) %>% 
  fn_surv_plot(.ylab = "Progression free survival probability", .title='Progression free survival data distribution') -> pfs.plot

ggsave(
  filename ='pfs-plot.pdf',
  plot = print(pfs.plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)

total434.pfs.se.duration <- total434.pfs.se
total434.pfs.se.duration$duration <- ifelse(total434.pfs.se.duration$duration > 60, 60, total434.pfs.se.duration$duration)

fn_get_duration_event(.se = total434.pfs.se.duration) %>% 
  fn_surv_plot(.ylab = "Progression free survival probability", .title='Progression free survival data distribution') -> pfs.duration.plot

ggsave(
  filename ='pfs-duration-plot.pdf',
  plot = print(pfs.duration.plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/output',
  width = 8,
  height = 9
)


# Save --------------------------------------------------------------------

save.image(file = 'data/rda/071-pfs-os.rda')
