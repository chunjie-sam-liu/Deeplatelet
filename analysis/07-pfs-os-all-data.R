
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
fn_get_duration_event <- function(.se, .type='os') {
  .cutoff <- ifelse(.type == 'os', 100, 60)
  .se@colData[,c('barcode', 'oc', 'duration', 'event')] %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(duration = ifelse(duration > .cutoff, .cutoff, duration))
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


fn_get_duration_event(.se = total416.os.se, .type = 'os') %>% 
  fn_surv_plot(.ylab = 'Overall survival probability', .title = 'Overall survival data distribution') -> os.plot

ggsave(
  filename ='OS-alldata-survival-plot.pdf',
  plot = print(os.plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/newoutput',
  width = 8,
  height = 9
)



# PFS ---------------------------------------------------------------------



fn_get_duration_event(.se = total434.pfs.se, .type = 'pfs') %>% 
  fn_surv_plot(.ylab = "Progression free survival probability", .title='Progression free survival data distribution') -> pfs.plot

ggsave(
  filename ='PFS-alldata-survival-plot.pdf',
  plot = print(pfs.plot, newpage = FALSE),
  device = 'pdf',
  path = 'data/newoutput',
  width = 8,
  height = 9
)


# Save --------------------------------------------------------------------

save.image(file = 'data/rda/07-pfs-os-all-data.rda')
