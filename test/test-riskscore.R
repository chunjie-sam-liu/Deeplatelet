
# Library -----------------------------------------------------------------

library(magrittr)
library(survival)
library(survminer)
library(survivalROC)
library(ggplot2)
library(doParallel)

# src ---------------------------------------------------------------------

source("src/utils.R")
source("src/doparallel.R")

# Function ----------------------------------------------------------------
fn_load <- function(.file) {
  .filepath = gsub('feather', 'tsv', glue::glue('data/rda/riskscore/{.file}'))
  # .d <- feather::read_feather(path = .filepath)
  .d <- readr::read_tsv(file = .filepath)
}

fn_bootstrap_cindex <- function(.d) {
  .d %>% 
    dplyr::sample_n(size = nrow(.d), replace = TRUE) -> 
    .dd
  .c <- survConcordance(Surv(durations, event) ~ riskscore, data = .dd)
  
  .c$concordance
}

fn_cindex_ci <- function(.d) {
  set.seed(1234)
  foreach(
    i = 1:500,
    .combine = rbind, 
    .packages = c('magrittr', 'survival', 'survivalROC'),
    .export = c('fn_bootstrap_cindex')
  ) %dopar% {
    fn_bootstrap_cindex(.d)
  } ->
    .boot
  .lower <- quantile(.boot[, 1], 0.025)
  .upper <- quantile(.boot[, 1], 0.975)
  glue::glue("{signif(.lower, 3)}-{signif(.upper, 3)}")
}

fn_surival_plot <- function(.d, .cohort = 'TC', .ylab = 'Overall survival probability') {
  
  .sd <- survdiff(formula = Surv(durations, event) ~ group, data = .d)
  .kmp <- 1 - pchisq(.sd$chisq, df = length(levels(as.factor(.d$group))) - 1)
  
  .cindex <- survConcordance(Surv(durations, event) ~ riskscore, data = .d)
  .cindex.ci <- fn_cindex_ci(.d)
  .cindex.label <- glue::glue("C-index: {signif(.cindex$concordance, digits = 3)} ({.cindex.ci})")
  
  .d %>% 
    dplyr::group_by(group) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(lab = glue::glue("{group} ({n})")) ->
    .legend.label
  
  ggsurvplot(
    fit = survfit(formula = Surv(durations, event) ~ group, data = .d),
    data = .d,
    pval = FALSE,
    pval.method = TRUE,
    pval.size = 7,
    palette = RColorBrewer::brewer.pal(n=3, name = 'Set1'),
    break.time.by = 20,
    ggtheme = theme(
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
    risk.table.col = 'strata',
    risk.table.fontsize = 5,
    risk.table.height = 0.25,
    
    ncensor.plot = FALSE,
    # surv.median.line = 'hv',
    
    legend = "top",
    legend.labs = .legend.label$group,
    legend.title = 'Risk',
    xlab = 'Time in months',
    ylab = .ylab
  ) ->
    .p
  
  .p$plot +
    theme(
      legend.position = c(0.235, 0.09),
      legend.direction = "horizontal"
    ) +
    annotate(geom = 'text', x = 0, y = 0.2, label = human_read_latex_pval(human_read(.kmp), .s = "Log-rank"), size = 5, hjust = 0, vjust = 0) +
    annotate(geom = 'text', x = 0, y = 0.12, label = .cindex.label, size = 5, hjust = 0, vjust = 0) +
    annotate(geom = 'text', x = 0, y = 0.29, label = .cohort, size = 5, hjust = 0, vjust = 0) ->
    .p$plot
  
  .p
}

fn_save_survival_plot <- function(.obj, .cohort, .type) {
  .filename <- glue::glue("test-final-{.type}-{.cohort}.pdf")
  ggsave(
    filename = .filename,
    plot = print(.obj, newpage = FALSE),
    device = 'pdf',
    path = "data/output",
    width = 6,
    height = 6.3
  )
}

fn_auc <- function(.d, .type) {
  .time <- ifelse(.type == 'os', 5 * 12, 3 * 12)
  survivalROC(
    Stime = .d$durations,
    status = .d$event,
    marker = .d$riskscore,
    predict.time = .time,
    method = 'NNE',
    span = 0.25 * nrow(.d)^(-0.2)
  )
}

fn_bootstrap_auc <- function(.d,.type) {
  .d %>% 
    dplyr::sample_n(size = nrow(.d), replace = TRUE) -> 
    .dd
  
  .sroc <- fn_auc(.dd, .type = .type)
  
  .sroc$AUC
  
}

fn_auc_ci <- function(.d, .type = .type) {
  set.seed(1234)
  foreach(
    i = 1:500, 
    .combine = rbind,
    .packages = c('magrittr', 'survival', 'survivalROC'),
    .export = c('fn_bootstrap_auc', 'fn_auc')
  ) %dopar% {
    fn_bootstrap_auc(.d, .type = .type)
  }->
    .boot
  .lower <- quantile(.boot[, 1], 0.025)
  .upper <- quantile(.boot[, 1], 0.975)
  glue::glue("{signif(.lower, 3)}-{signif(.upper, 3)}")
  
}

# Load pfs data -----------------------------------------------------------

pfs.train <- fn_load(.file = 'pfs.train.test.feather')
pfs.val <- fn_load(.file = 'pfs.val.test.feather')

pfs.merge <- pfs.train %>% 
  dplyr::bind_rows(pfs.val) %>% 
  dplyr::mutate(group = ifelse(riskscore > median(riskscore), 'High', 'Low'))

pfs.test1 <- fn_load(.file = 'pfs.test1.test.feather') %>% 
  dplyr::mutate(group = ifelse(riskscore > median(pfs.merge$riskscore), 'High', 'Low'))

pfs.test2 <- fn_load(.file = 'pfs.test2.test.feather') %>% 
  dplyr::mutate(group = ifelse(riskscore > median(pfs.merge$riskscore), 'High', 'Low'))


# pfs survival plot -------------------------------------------------------

fn_surival_plot(.d = pfs.merge, 'TC', .ylab = "Progression free survival probability") %>% 
  fn_save_survival_plot(.cohort = 'TC', .type = 'pfs')
fn_surival_plot(.d = pfs.test1, 'EV1', .ylab = "Progression free survival probability") %>%
  fn_save_survival_plot(.cohort = 'EV1', .type = 'pfs')
fn_surival_plot(.d = pfs.test2, 'EV2', .ylab = "Progression free survival probability") %>%
  fn_save_survival_plot(.cohort = 'EV2', .type = 'pfs')


# PFS ---------------------------------------------------------------------

pfs.merge.auc <- fn_auc(pfs.merge, .type = 'pfs')
pfs.merge.auc$ci <- fn_auc_ci(.d = pfs.merge, .type = 'pfs')
pfs.merge.auc$tb <- pfs.merge.auc[c('TP', 'FP')] %>% 
  as.data.frame() %>% 
  dplyr::mutate(cohort = 'TC')
pfs.test1.auc <- fn_auc(pfs.test1, .type = 'pfs')
pfs.test1.auc$ci <- fn_auc_ci(.d = pfs.test1, .type = 'pfs')
pfs.test1.auc$tb <- pfs.test1.auc[c('TP', 'FP')] %>% 
  as.data.frame() %>% 
  dplyr::mutate(cohort = 'EV1')
pfs.test2.auc <- fn_auc(pfs.test2, .type = 'pfs')
pfs.test2.auc$ci <- fn_auc_ci(.d = pfs.test2, .type = 'pfs')
pfs.test2.auc$tb <- pfs.test2.auc[c('TP', 'FP')] %>% 
  as.data.frame() %>% 
  dplyr::mutate(cohort = 'EV2')


pfs.legend.labels <- c(
  glue::glue("TC   {signif(pfs.merge.auc$AUC, 3)} ({pfs.merge.auc$ci})"),
  glue::glue("EV1 {signif(pfs.test1.auc$AUC, 3)} ({pfs.test1.auc$ci})"),
  glue::glue("EV2 {signif(pfs.test2.auc$AUC, 3)} ({pfs.test2.auc$ci})")
)

pfs.merge.auc$tb %>% 
  dplyr::bind_rows(pfs.test1.auc$tb) %>% 
  dplyr::bind_rows(pfs.test2.auc$tb) %>% 
  dplyr::mutate(cohort = factor(cohort, levels = c('TC', 'EV1', 'EV2'))) %>% 
  ggplot(aes(x = FP, y = TP, color = cohort)) +
  geom_path(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 11) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(
    name = 'AUC, 3 years',
    labels = pfs.legend.labels,
    values = RColorBrewer::brewer.pal(n=3, name = 'Set1')[c(2, 1, 3)]
  ) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    axis.line.x.bottom = element_line(color = 'black'),
    axis.line.y.left = element_line(color = 'black'),
    axis.ticks.length = unit(x = 0.2, units = 'cm'),
    axis.text = element_text(color = 'black', size = 14),
    axis.title = element_text(color = 'black', size = 16),
    
    legend.position = c(0.68, 0.2),
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(1.5, units = 'cm'),
    legend.spacing = unit(c(0,0,0,0), units = 'cm'),
    legend.title.align = 1,
    
    plot.margin = unit(c(1,1,0.5,0.5), units = 'cm'),
    plot.title = element_text(hjust = 0.5, size = 18)
  ) + 
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = "Progression free surival"
  ) ->
  pfs.auc.3years.plot;pfs.auc.3years.plot

ggsave(
  filename = 'data/output/test-final-pfs-3years-auc.pdf',
  plot = pfs.auc.3years.plot,
  device = 'pdf',
  width = 6,
  height = 5
)

# -------------------------------------------------------------------------

fn_parallel_stop()

# Save image --------------------------------------------------------------
save.image(file = 'data/rda/11-riskscore.rda')
