
# Library -----------------------------------------------------------------

library(magrittr)
library(survival)
library(survminer)
library(survivalROC)
library(ggplot2)

# Function ----------------------------------------------------------------
fn_load <- function(.file) {
  .filepath = glue::glue('data/rda/riskscore/{.file}')
  .d <- feather::read_feather(path = .filepath)
}

fn_surival_plot_function(.d, .ylab = 'Overall survival probability', .title = 'Overall survival data distribution') {
  
    ggsurvplot(
    fit = survfit(formula = Surv(durations, event) ~ group, data = .d),
    data = .d,
    pval = TRUE,
    pval.mehtod = TRUE,
    palette = RColorBrewer::brewer.pal(n=3, name = 'Set1'),
    break.time.by = 20,
    ggtheme = theme(
      
      panel.background = element_rect(fill = NA, colour = 'black'),
      
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16, face='bold', color = 'black')
      
    ),
    
    risk.table = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    risk.table.col = 'strata',
    risk.table.fontsize = 6,
    
    ncensor.plot = FALSE,
    # surv.median.line = 'hv',
    
    legend = 'top',
    # legend.labs = .legend,
    legend.title = 'Group',
    xlab = 'Time in months',
    ylab = .ylab,
    title = .title
  )
}

# Load data ---------------------------------------------------------------

os.train <- fn_load(.file = 'os.train.feather')
os.train$group <- ifelse(os.train$riskscore > median(os.train$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = os.train)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = os.train),
  data = os.train,
  pval = T,
  title = 'TC data high low risk'
)

os.val <- fn_load(.file = 'os.val.feather')
os.val$group <- ifelse(os.val$riskscore > median(os.val$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = os.val)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = os.val),
  data = os.val,
  pval = T,
  title = 'TC data high low risk'
)

os.merge <- dplyr::bind_rows(os.train, os.val)
os.merge$group <- ifelse(os.merge$riskscore > median(os.train$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = os.merge)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = os.merge),
  data = os.merge,
  pval = T,
  title = 'TC data high low risk'
)


os.test1 <- fn_load(.file = 'os.test1.feather')
os.test1$group <- ifelse(os.test1$riskscore > median(os.train$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = os.test1)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = os.test1),
  data = os.test1,
  pval = T,
  title = 'TC data high low risk'
)

os.test2 <- fn_load(.file = 'os.test2.feather')
os.test2$group <- ifelse(os.test2$riskscore > median(os.train$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = os.test2)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = os.test2),
  data = os.test2,
  pval = T,
  title = 'TC data high low risk'
)


pfs.train <- fn_load(.file = 'pfs.train.feather')
pfs.train$group <- ifelse(pfs.train$riskscore > median(pfs.train$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = pfs.train)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = pfs.train),
  data = pfs.train,
  pval = T,
  title = 'TC data high low risk'
)

pfs.val <- fn_load(.file = 'pfs.val.feather')
pfs.val$group <- ifelse(pfs.val$riskscore > median(pfs.val$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = pfs.val)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = pfs.val),
  data = pfs.val,
  pval = T,
  title = 'TC data high low risk'
)

pfs.merge <- dplyr::bind_rows(pfs.train, pfs.val)
pfs.merge$group <- ifelse(pfs.merge$riskscore > median(pfs.train$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = pfs.merge)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = pfs.merge),
  data = pfs.merge,
  pval = T,
  title = 'TC data high low risk'
)


pfs.test1 <- fn_load(.file = 'pfs.test1.feather')
pfs.test1$group <- ifelse(pfs.test1$riskscore > median(pfs.train$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = pfs.test1)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = pfs.test1),
  data = pfs.test1,
  pval = T,
  title = 'TC data high low risk'
)

pfs.test2 <- fn_load(.file = 'pfs.test2.feather')
pfs.test2$group <- ifelse(pfs.test2$riskscore > median(pfs.train$riskscore), 'high', 'low')
survConcordance(Surv(durations, event) ~ riskscore, data = pfs.test2)
ggsurvplot(
  fit = survfit(formula = Surv(durations, event) ~ group, data = pfs.test2),
  data = pfs.test2,
  pval = T,
  title = 'TC data high low risk'
)

