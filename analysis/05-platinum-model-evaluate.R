
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

fn_se2total_task <- function(.se, .feats, .dataset = list(train = 'OC521', eval = 'OC44', test172 = 'OC172', test79 = 'OC79')) {
  # SE to task
  
  # Convert all se to mlr task
  .task.panel <- fn_se2task(.se = .se[.feats, ], .id = 'Panel-task')
  
  # Samples for panel
  .samples.panel <- purrr::map(.x = .dataset, .f = fn_task_inds, .se = .se[.feats, ])
  names(.samples.panel) <- glue::glue('panel.{names(.samples.panel)}')
  
  # Filter samples with no ca125
  .se_ca125 <- .se[.feats, !is.na(.se@colData$CA125)]
  .task.ca125 <- fn_se2task_ca125(.se = .se_ca125, .id = 'CA125-task')
  .task.panel.ca125 <- fn_se2task_panel_ca125(.se = .se_ca125, .id = 'Panel-CA125-task')
  
  # Samples for ca125
  .samples.ca125 <- purrr::map(.x = .dataset, .f = fn_task_inds, .se = .se_ca125[.feats, ])
  names(.samples.ca125) <- glue::glue('ca125.{names(.samples.ca125)}')
  
  # Return data
  list(
    task = list(
      panel = .task.panel,
      ca125 = .task.ca125,
      panel.ca125 = .task.panel.ca125
    ),
    samples = list(
      panel = .samples.panel,
      ca125 = .samples.ca125,
      panel.ca125 = .samples.ca125
    )
  )
}

fn_plot_tune_path <- function(.tune_result, .task_id) {
  .task_name <- gsub(pattern = '-task', replacement = '', x = .task_id)
  # save the random search iteration
  .hped <-  mlr::generateHyperParsEffectData(tune.result = .tune_result, trafo = TRUE)
  .plot_hpe <- mlr::plotHyperParsEffect(
    hyperpars.effect.data = .hped,
    x = "iteration", y = "auc.test.mean",
    plot.type = "line"
  ) +
    theme_bw() +
    theme() +
    guides(shape = guide_legend(title = 'Learner Status'), colour = guide_legend(title = 'Learner Status')) +
    labs(
      x = 'Iteration',
      y = 'Accuracy test mean',
      title = glue::glue('Random search iteration for training {.task_name}')
    )
  # save ifs to plot
  ggsave(
    filename = glue::glue('02-Tune-parameter-{.task_name}.pdf'),
    plot = .plot_hpe,
    device = 'pdf',
    path = 'data/output',
    width = 8,
    height = 4
  )
}

fn_tune_hyperparameters <- function(.task, .samples) {
  .task_id <- getTaskId(x = .task)
  .task_for_tunes <- mlr::subsetTask(task = .task, subset = c(.samples[[1]], .samples[[3]]))
  .task_rownames <- rownames(mlr::getTaskData(task = .task_for_tunes))
  .train.inds <- match(x = names(.samples[[1]]), .task_rownames)
  .test.inds <- match(x = names(c(.samples[[3]])), .task_rownames)
  
  # fix the train and test samples with sample index.
  .holdout <-  mlr::makeFixedHoldoutInstance(
    train.inds = .train.inds, test.inds = .test.inds,
    size = mlr::getTaskSize(.task_for_tunes)
  )
  
  # using svm learner with predict type 'prob'
  .learner <- mlr::makeLearner(cl = 'classif.ksvm', id = 'svm-learner', predict.type = 'prob')
  
  # hypterparameters search space
  .hypterparameters <- ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("C", lower = -10, upper = 10, trafo = function(x) 10 ^ x),
    ParamHelpers::makeNumericParam("sigma", lower = -10, upper = 10, trafo = function(x) 10 ^ x)
  )
  
  # hypterparameter optimize algorithm Random tuning
  .tune_algorithm <-  mlr::makeTuneControlRandom(
    same.resampling.instance = TRUE, maxit = 5000L
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
    resampling = .holdout,
    measures = list(mlr::auc, mlr::mmce, mlr::acc),
    par.set = .hypterparameters,
    control = .tune_algorithm
  )
  
  # plot random search and test mean auc
  fn_plot_tune_path(.tune_result = .tune_result, .task_id = .task_id)
  
  # tuned learner
  .leaner_tuned <- mlr::setHyperPars(learner = .learner, par.vals = .tune_result$x)
  .model_tuned <-  mlr::train(learner = .leaner_tuned, task = .task, subset = .samples[[1]])
  
  # .d <- predict(object = .model_tuned, task = .total_task$task$ca125, subset = .subset$ca125)
  
  return(.model_tuned)
}

fn_performace_dataset <- function(.dataset, .models, .total_task) {
  # swith datasets
  .subset <- fn_dataset_switch(.dataset = .dataset, .total_task = .total_task)
  # predict dataset with three models
  .pred <- list(
    panel = predict(object = .models$panel, task = .total_task$task$panel, subset = .subset$panel),
    ca125 = predict(object = .models$ca125, task = .total_task$task$ca125, subset = .subset$ca125),
    panel.ca125 = predict(object = .models$panel.ca125, task = .total_task$task$panel.ca125, subset = .subset$panel.ca125)
  )
  # plot and save
  .plot_auroc <- fn_plot_auroc(.pred = .pred, .subset = .subset, .dataset = .dataset)
  # save the metrics ci95 tables
  .table_metrics <- fn_metrics(.pred = .pred, .dataset = .dataset)
  # return auroc plot and metrics table
  list(auroc = .plot_auroc, metrics = .table_metrics)
}


# Prepare task ------------------------------------------------------------

# panel <- feats
panel <- readr::read_rds(file = 'data/rda/00-selected-features.rds.gz')
total351.task.list <- fn_se2total_task(.se = total351.platinum.se.norm, .feats = panel$panel)
readr::write_rds(x = total351.task.list, file = 'data/rda/00-total351.task.list.rds.gz', compress = 'gz')

# Modeling ----------------------------------------------------------------
# use 50 number of cores
parallelMap::parallelStart(mode = 'multicore', cpus = 100)
models <- list(
  # model with panel
  panel = fn_tune_hyperparameters(
    .task = total351.task.list$task$panel,
    .samples = total351.task.list$samples$panel
  ),
  # model with ca125
  ca125 = fn_tune_hyperparameters(
    .task = total351.task.list$task$ca125,
    .samples = total351.task.list$samples$ca125
  ),
  # model with panel and ca125
  panel.ca125 = fn_tune_hyperparameters(
    .task = total351.task.list$task$panel.ca125,
    .samples = total351.task.list$samples$panel.ca125
  )
)
parallelMap::parallelStop()

readr::write_rds(x = models, file = 'data/rda/00-model.rds.gz', compress = 'gz')


# Performance -------------------------------------------------------------

dataset_oc <- list(OC521 = 'OC521', OC44 = 'OC44', OC172 = 'OC172', OC79 = 'OC79')
performace_dataset <- dataset_oc %>%
  purrr::map(.f = fn_performace_dataset, .models = models, .total_task = total351.task.list)

# Save image --------------------------------------------------------------
save.image(file = 'data/rda/05-platinum-model-evaluate.rda', ascii = FALSE, compress = TRUE)
