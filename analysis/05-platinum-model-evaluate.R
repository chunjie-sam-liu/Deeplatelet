# Library -----------------------------------------------------------------

library(magrittr)
library(assertthat)
library(DESeq2)
library(mlr)
library(mRMRe)
library(caret)
library(doParallel)
library(assertthat)

# src ---------------------------------------------------------------------

source(file = 'src/doparallel.R', local = TRUE)
source(file = 'src/utils.R', local = TRUE)
source(file = 'src/model-performance.R', local = TRUE)

# Load data ---------------------------------------------------------------
total351.platinum.se.norm <- readr::read_rds(file = 'data/rda/total351.platinum.se.norm.rds.gz')

# Function ----------------------------------------------------------------
fn_se2task <- function(.se) {
  # Task
  .d <- cbind(
    t(assay(.se)),
    data.frame(platinum = factor(x = as.character(.se@colData[, 'platinum']), levels = c('sensitive', 'resistant')))
  )
  
  .task <- mlr::makeClassifTask(id = "Panel-task", data = .d, target = "platinum", positive = "sensitive")
  
  # Partition
  .dataset = list(train = 'OC521', eval = 'OC44', test172 = 'OC172', test79 = 'OC79')
  .samples.panel <- purrr::map(.x = .dataset, .f = function(.x, .se) {
    .barcode_ind <- which(.se$oc == .x)
    setNames(object = .barcode_ind, nm = .se@colData[.barcode_ind, 'barcode'])
  }, .se = .se)
  
  .samples.panel$merge <- c(.samples.panel$train, .samples.panel$eval)
  list(task = .task, samples = .samples.panel)
}
fn_plot_tune_path <- function(.tune_result, .task_id) {
  .task_name <- gsub(pattern = "-task", replacement = "", x = .task_id)
  .hped <- mlr::generateHyperParsEffectData(tune.result = .tune_result, trafo = TRUE)
  
  .plot_hpe <- mlr::plotHyperParsEffect(
    hyperpars.effect.data = .hped,
    x = "iteration", y = "acc.test.mean",
    plot.type = "line"
  ) +
    theme_bw() +
    theme() +
    guides(shape = guide_legend(title = "Learner Status"), colour = guide_legend(title = "Learner Status")) +
    labs(
      x = "Iteration",
      y = "Accuracy test mean",
      title = glue::glue("Random search iteration for training {.task_name}")
    )
  ggsave(
    filename = glue::glue("01-Tune-parameter-{.task_name}.pdf"),
    plot = .plot_hpe,
    device = "pdf",
    path = "data/newoutput",
    width = 8,
    height = 4
  )
  
  .plot_hpe
}

fn_tune_hyperparameters <- function(.list) {
  .task <- .list$task
  .samples <- .list$samples
  .task_id <- getTaskId(x = .task)
  
  .task_for_tunes <- mlr::subsetTask(task = .task, subset = .samples$merge)
  
  set.seed(123)
  .cv10d <- mlr::makeResampleDesc(method = "CV", iters = 3, stratify = TRUE)
  .cv10i <- mlr::makeResampleInstance(desc = .cv10d, task = .task_for_tunes)
  .learner <- mlr::makeLearner(
    cl = "classif.ksvm",
    id = 'svm-learner',
    predict.type = 'prob'
  )
  .hypterparameters <- ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("C", lower = -10, upper = 10, trafo = function(x) 10 ^ x),
    ParamHelpers::makeNumericParam("sigma", lower = -10, upper = 10, trafo = function(x) 10 ^ x)
  )
  .tune_algorithm <-  mlr::makeTuneControlRandom(
    same.resampling.instance = TRUE, maxit = 50L
  )
  mlr::configureMlr(
    show.info = FALSE,
    on.learner.error = 'warn',
    on.measure.not.applicable = 'warn'
  )
  parallelMap::parallelStart(mode = 'multicore', cpus = 50L)
  .tune_result <-  mlr::tuneParams(
    learner = .learner,
    task = .task_for_tunes,
    resampling = .cv10i,
    measures = list(mlr::acc, mlr::mmce, mlr::auc),
    par.set = .hypterparameters,
    control = .tune_algorithm
  )
  parallelMap::parallelStop()
  .learner <- mlr::setHyperPars(learner = .learner, par.vals = .tune_result$x)
  fn_plot_tune_path(.tune_result = .tune_result, .task_id = .task_id)
  
  mlr::train(
    learner = .learner,
    task = .task,
    subset = .samples$merge
  )

}

fn_pred_perf_metrics <- function(.x, .task, .samples, .model) {
  .pred <- predict(
    object = .model,
    task = .task,
    subset = .samples[[.x]]
  )
  
  .pred_perf <- mlr::generateThreshVsPerfData(obj = .pred, measures = list(mlr::fpr, mlr::tpr, mlr::auc), gridsize = length(.samples[[.x]])) %>% 
    .$data %>% 
    tibble::add_column(cohort = .x, .before = 1)
  
  .pred_metrics <- fn_roc_95ci(.pred) %>% 
    dplyr::mutate_if(
      .predicate = rlang::is_list,
      .funs = ~purrr::map(.x = ., .f = function(.x) {glue::glue('{sprintf("%.3f", .x[1])} ({sprintf("%.3f", .x[2])} - {sprintf("%.3f", .x[3])})')})) %>% 
    dplyr::mutate_if(
      .predicate = rlang::is_double,
      .funs = ~purrr::map(.x = ., .f = function(.x) {sprintf("%.3f", .x)})
    ) %>%
    dplyr::mutate_if(
      .predicate = rlang::is_list,
      .funs = unlist
    ) %>% 
    tibble::add_column(cohort = .x, .before = 1)
  # names(.pred_metrics)  <- c("cohort", 'Accuracy (95% CI)', 'SN (95% CI)', 'SP (95% CI)', 'PPV (95% CI)', 'NPV (95% CI)', 'Kappa', 'F1')
  names(.pred_metrics) <- c("cohort", "AUC (95% CI)", "Accuracy (95% CI)", "SN (95% CI)", "SP (95% CI)", "PPV (95% CI)", "NPV (95% CI)", "Kappa", "F1")
  
  list(
    perf = .pred_perf,
    metrics = .pred_metrics
  )
}

fn_performance <- function(.list, .model) {
  .task <- .list$task
  .samples <- .list$samples
  .model <- .model
  
  .model_tuned_pred <- purrr::map(
    names(.samples),
    fn_pred_perf_metrics,
    .task = .task,
    .samples = .samples,
    .model = .model
    )
  names(.model_tuned_pred) <- names(.samples)
  
  list(
    model_tuned_pred = .model_tuned_pred
  )
}


# Panel -------------------------------------------------------------------

panel <- readr::read_rds(file = 'data/rda/panel.rds.gz')
total351.task.list <- fn_se2task(total351.platinum.se.norm[panel, ])
readr::write_rds(x = total351.task.list, file = 'data/rda/00-total351.task.list.rds.gz', compress = 'gz')


# Modeling ----------------------------------------------------------------

model <- fn_tune_hyperparameters(.list = total351.task.list)

readr::write_rds(x = model, file = 'data/rda/panel.model.rds.gz', compress = 'gz')


# Performance -------------------------------------------------------------

perf <- fn_performance(.list = total351.task.list, .model = model)

perf$model_tuned_pred %>% 
  purrr::map('metrics') %>% 
  purrr::reduce(.f = dplyr::bind_rows) %>% 
  dplyr::filter(cohort %in% c("merge", "test79", "test172")) %>% 
  dplyr::slice(3, 2, 1) %>% 
  dplyr::mutate(cohort = plyr::revalue(x = cohort, replace = c("merge" = "TC  ", "test79" = "EV1", "test172" = "EV2"))) %>% 
  dplyr::mutate(label = glue::glue("{cohort} {`AUC (95% CI)`}")) ->
  model_tuned_pred_label

model_tuned_pred_label %>% 
  dplyr::select(-label) %>% 
  readr::write_tsv(file = 'data/newoutput/merge-platinum-metrics.tsv')

model_tuned_pred_label %>% 
  dplyr::select(-label) %>% 
  writexl::write_xlsx(path = 'data/newoutput/merge-platinum-metrics.xlsx')

perf$model_tuned_pred %>% 
  purrr::map('perf') %>% 
  purrr::reduce(.f = dplyr::bind_rows) %>% 
  dplyr::filter(cohort %in% c("merge", "test79", "test172")) %>% 
  dplyr::mutate(cohort = plyr::revalue(x = cohort, replace = c("merge" = "TC", "test79" = "EV1", "test172" = "EV2"))) %>% 
  dplyr::mutate(cohort = factor(cohort, levels = c('TC', 'EV1', 'EV2'))) %>% 
  ggplot(aes(x = fpr, y = tpr, color = cohort)) +
  geom_path(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 11) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(
    name = 'AUC',
    labels = model_tuned_pred_label$label,
    values = RColorBrewer::brewer.pal(n=3, name = 'Set1')[c(2, 1, 3)]
  )  +
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
    title = "Platinum sensitivity"
  ) ->
  platinum_sensitivity_plot;platinum_sensitivity_plot

ggsave(
  filename = 'data/newoutput/merge-platinum-sensitivity.pdf',
  plot = platinum_sensitivity_plot,
  device = 'pdf',
  width = 6,
  height = 5
)


# Save image --------------------------------------------------------------

save.image(file = 'data/rda/05-platinum-model-evaluate.rda')
