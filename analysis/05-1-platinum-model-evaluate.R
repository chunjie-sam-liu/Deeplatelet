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

fn_tune_hyperparameters <- function(.list, .panel) {
  .task <- .list$task
  .samples <- .list$samples
  .task_id <- getTaskId(x = .task)
  
  .task_for_tunes <- mlr::subsetTask(task = .task, subset = .samples$train)
  
  # fix the train and test samples with sample index
  set.seed(123)
  .cv10d <- mlr::makeResampleDesc(method = "CV", iters = 10, stratify = TRUE)
  .cv10i <- mlr::makeResampleInstance(desc = .cv10d, task = .task_for_tunes)
  
  # using svm learner with predict type 'prob'
  .learner <- mlr::makeLearner(cl = "classif.ksvm", id = 'svm-learner', predict.type = 'prob')
  
  # hypterparameters search space
  .hypterparameters <- ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("C", lower = -10, upper = 10, trafo = function(x) 10 ^ x),
    ParamHelpers::makeNumericParam("sigma", lower = -10, upper = 10, trafo = function(x) 10 ^ x)
  )
  
  # hypterparameter optimize algorithm Random tuning
  .tune_algorithm <-  mlr::makeTuneControlRandom(
    same.resampling.instance = TRUE, maxit = 500L
  )
  # configue
  mlr::configureMlr(
    show.info = FALSE,
    on.learner.error = 'warn',
    on.measure.not.applicable = 'warn'
  )
  
  # tuning
  .tune_result <-  mlr::tuneParams(
    learner = .learner, task = .task_for_tunes,
    resampling = .cv10i,
    measures = list(mlr::auc, mlr::mmce, mlr::acc),
    par.set = .hypterparameters,
    control = .tune_algorithm
  )
  
  .hped <- mlr::generateHyperParsEffectData(.tune_result, trafo = TRUE)
  .plot_hpe <- mlr::plotHyperParsEffect(.hped, x = "iteration", y = "auc.test.mean", plot.type = "line")
  .learner_tuned <- mlr::setHyperPars(learner = .learner, par.vals = .tune_result$x)
  .model_tuned <- mlr::train(learner = .learner_tuned, task = .task, subset = .samples$train)
  
  .model_old <- mlr::train(learner = mlr::setHyperPars(learner = .learner, par.vals = as.list(.panel$hyperparameter)), task = .task, subset = .samples$train)
  
  list(
    hped = .hped,
    model_tuned = .model_tuned,
    hyperparameter = .tune_result$x,
    model_old = .model_old
  )
}

fn_pred_perf_metrics <- function(.x, .model) {
  .pred <- predict(object = .model, task = .task, subset = .samples[[.x]])
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
  names(.pred_metrics)  <- c("cohort", 'Accuracy (95% CI)', 'SN (95% CI)', 'SP (95% CI)', 'PPV (95% CI)', 'NPV (95% CI)', 'Kappa', 'F1')
  
  list(
    perf = .pred_perf,
    metrics = .pred_metrics
  )
}
fn_performance <- function(.list, .model) {
  .task <- .list$task
  .samples <- .list$samples
  .model_tuned <- .model$model_tuned
  .model_old <- .model$model_old
  
  .model_tuned_pred <- purrr::map(names(.samples), fn_pred_perf_metrics, .model = .model_tuned)
  names(.model_tuned_pred) <- names(.samples)
  
  .model_old_pred <- purrr::map(names(.samples), fn_pred_perf_metrics, .model = .model_old)
  names(.model_old_pred) <- names(.samples)
  
  list(
    model_tuned_pred = .model_tuned_pred,
    model_old_pred = .model_old_pred
  )
}


# Panel -------------------------------------------------------------------

panel <- readr::read_rds(file = 'data/rda/00-selected-features.rds.gz')
total351.task.list <- fn_se2task(total351.platinum.se.norm[panel$panel, ])
readr::write_rds(x = total351.task.list, file = 'data/rda/00-total351.task.list.rds.gz', compress = 'gz')


# Modeling ----------------------------------------------------------------

parallelMap::parallelStart(mode = 'multicore', cpus = 100)

#model <- fn_tune_hyperparameters(.list = total351.task.list, .panel = panel)

parallelMap::parallelStop()
readr::write_rds(x = model, file = 'data/rda/00-model-test.rds.gz', compress = 'gz')



# Performance -------------------------------------------------------------

perf <- fn_performance(.list = total351.task.list, .model = model)

perf$model_tuned_pred %>% 
  purrr::map('metrics') %>% 
  purrr::reduce(.f = dplyr::bind_rows) %>% 
  dplyr::filter(cohort %in% c("merge", "test79", "test172")) %>% 
  dplyr::slice(3, 2, 1) %>% 
  dplyr::mutate(cohort = plyr::revalue(x = cohort, replace = c("merge" = "TC  ", "test79" = "EV1", "test172" = "EV2"))) %>% 
  dplyr::mutate(label = glue::glue("{cohort} {`Accuracy (95% CI)`}")) ->
  model_tuned_pred_label

model_tuned_pred_label %>% 
  dplyr::select(-label) %>% 
  readr::write_tsv(file = 'data/output/final-platinum-sensitivity.tsv')

model_tuned_pred_label %>% 
  dplyr::select(-label) %>% 
  writexl::write_xlsx(path = 'data/output/final-platinum-sensitivity.xlsx')

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
  platinum_sensitivity_plot

ggsave(
  filename = 'data/output/final-platinum-sensitivity.pdf',
  plot = platinum_sensitivity_plot,
  device = 'pdf',
  width = 6,
  height = 5
)




perf$model_old_pred %>% 
  purrr::map('perf') %>% 
  purrr::reduce(.f = dplyr::bind_rows)
perf$model_old_pred %>% 
  purrr::map('metrics') %>% 
  purrr::reduce(.f = dplyr::bind_rows)


# Save image --------------------------------------------------------------

save.image(file = 'data/rda/051-platinum-model-evaluate.rda')
