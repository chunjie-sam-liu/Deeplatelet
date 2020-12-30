# Library -----------------------------------------------------------------

library(magrittr)
library(assertthat)
library(doParallel)
library(DESeq2)
library(mlr)
library(mRMRe)
library(caret)
library(glmnet)

# src ---------------------------------------------------------------------

source(file = 'src/doparallel.R', local = TRUE)
source(file = 'src/utils.R', local = TRUE)
source(file = 'src/model-performance.R', local = TRUE)
# Load data ---------------------------------------------------------------

total351.platinum.se.norm <- readr::read_rds(file = 'data/rda/total351.platinum.se.norm.rds.gz')

# Functions ---------------------------------------------------------------

# Feature selection
fn_ifs_eval <- function(.x, .feats, .train.se, .eval.se) {
  .feats <- head(.feats, .x)
  .train.feats.se <- .train.se[.feats, ]
  .eval.feats.se <- .eval.se[.feats, ]
  
  # 2. build model
  .train.feats.model <- fn_caret_train_model(.se = .train.feats.se)
  .bestune <- .train.feats.model$bestTune
  
  # 3. ensembl predict
  .train.auc <- fn_ensembl_predict(.model = .train.feats.model, .se = .train.feats.se, .title = 'Train')
  .eval.auc <- fn_ensembl_predict(.model = .train.feats.model, .se = .eval.feats.se, .title = 'Eval')
  
  # 4. return
  tibble::tibble(
    num = .x,
    train = .train.auc$alist$auc,
    eval = .eval.auc$alist$auc,
    train_error = list(rownames(.train.auc$predict)[.train.auc$predict$correct == FALSE]),
    eval_error = list(rownames(.eval.auc$predict)[.eval.auc$predict$correct == FALSE]),
    best_tune = list(.bestune)
  )
}
fn_ifs_plot <- function(.ifs.res){
  .ifs.res %>%
    dplyr::select_if(.predicate = purrr::negate(rlang::is_list)) %>%
    tidyr::gather(key = 'dataset', value = 'auc', -num) %>%
    ggplot(aes(x = num, y = auc, color = dataset)) +
    geom_line() +
    scale_color_manual(
      name = 'Dataset',
      limits = c('train', 'eval'),
      labels = c('Training set', 'Evaluation Set'),
      values = ggthemes::gdocs_pal()(5)[c(2, 5, 3, 4, 1)],
    ) +
    theme_bw() +
    theme() +
    labs(
      x = 'Number of genes in panel',
      y = 'AUC',
      title = 'Increment Feature Selection'
    ) ->
    .ifs.res.plot
  .ifs.res.plot
}
fn_select_features <- function(.se, .train_sample, .eval_sample) {
  .train.se <- .se[, .se$oc == .train_sample]
  .eval.se <- .se[, .se$oc == .eval_sample]
  
  .df <- fn_se2data.frame(.se = .train.se)
  # 1. standard variation filter
  .feats.sd.median <- fn_sd_median(.df = .df)
  
  # 3. lasso expression L1-norm
  # Using lasso to select features based on big data.
  set.seed(1223)
  .lasso <- caret::train(
    platinum ~., data = .df[, c('platinum', .feats.sd.median)],
    method = "glmnet", family = 'binomial',
    trControl = trainControl("cv", number = 10),
    tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-3, 3, length = 100))
  )
  
  .lasso.dgcmatrix <- coef(.lasso$finalModel, .lasso$bestTune$lambda)
  as.matrix(.lasso.dgcmatrix) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'Gene_ID') %>%
    tibble::as_tibble() %>%
    dplyr::filter(abs(`1`) > 0, grepl(pattern = 'ENSG', x = Gene_ID)) %>%
    dplyr::arrange(-abs(`1`)) %>%
    dplyr::pull(Gene_ID) ->
    .feats.lasso
  readr::write_rds(x = .feats.lasso, file = 'data/rda/00-lasso-feature.rds.gz', compress = 'gz')
  
  # 4. mrmr
  
  .df.elastic.lasso <- .df[, c('platinum', .feats.lasso)]
  
  .df.elastic.lasso.mrmr.dd <- mRMRe::mRMR.data(data = .df.elastic.lasso)
  
  .df.elastic.lasso.mrmr.dd.dd.filter <- mRMRe::mRMR.ensemble(
    data = .df.elastic.lasso.mrmr.dd,
    target_indices = c(1),
    solution_count = 1,
    feature_count = length(.feats.lasso)
  )
  
  .feats.mrmr <- .df.elastic.lasso.mrmr.dd.dd.filter@feature_names[.df.elastic.lasso.mrmr.dd.dd.filter@filters$`1`[, 1]]
  
  .feats.mrmr_matrix <- cbind(.df.elastic.lasso.mrmr.dd.dd.filter@filters$`1`, .df.elastic.lasso.mrmr.dd.dd.filter@scores$`1`)
  
  rownames(.feats.mrmr_matrix) <- .feats.mrmr
  colnames(.feats.mrmr_matrix) <- c('order', 'score')
  
  .feats.mrmr_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'Gene_ID') %>%
    tibble::as_tibble() %>%
    dplyr::arrange(-score) %>%
    dplyr::pull(Gene_ID) ->
    .feats.mrmr_sort
  
  # 5. ifs for evaluation
  .ifs.res <- foreach(
    i = 10:length(.feats.mrmr_sort),
    .combine = rbind,
    .packages = c('magrittr', 'DESeq2', 'ggplot2'),
    .export = c('fn_ifs_eval', 'fn_caret_train_model', 'fn_ensembl_predict', 'fn_predict', 'fn_getROC_AUC', 'fn_alist', 'fn_plot_auc', 'fn_se2data.frame')
  ) %dopar% {
    fn_ifs_eval(.x = i, .feats = .feats.mrmr_sort, .train.se = .train.se, .eval.se = .eval.se)
  }
  # plot ifs
  .ifs.res.plot <- fn_ifs_plot(.ifs.res = .ifs.res)
  # save ifs to plot
  ggsave(
    filename = '01-Training-evaluation-on-selected-features.pdf',
    plot = .ifs.res.plot,
    device = 'pdf',
    path = 'data/output',
    width = 8,
    height = 4
  )
  # .optimal <- .ifs.res %>%
  #   dplyr::filter(eval > 0.95) %>%
  #   dplyr::mutate(sub = eval - train) %>%
  #   dplyr::arrange(-sub)
  .optimal <- .ifs.res %>%
    dplyr::arrange(-eval)
  
  .optimal_feats_num <- .optimal %>% dplyr::pull(num) %>% head(1)
  .optimal_parameter <- .optimal %>% dplyr::pull(best_tune) %>% head(1) %>% .[[1]]
  
  list(
    panel = head(.feats.mrmr_sort, .optimal_feats_num),
    hyperparameter = .optimal_parameter
  )
}



# Select features ---------------------------------------------------------

# start parallel
fn_parallel_start(n_cores = 50)
# feature selection
feats <- fn_select_features(.se = total351.platinum.se.norm, .train_sample = 'OC521', .eval_sample = 'OC44')

fn_parallel_stop()

readr::write_rds(x = feats, file = 'data/rda/00-selected-features.rds.gz', compress = 'gz')
# Save image --------------------------------------------------------------

save.image(file = 'data/rda/04-platinum-feature.rda', ascii = FALSE, compress = TRUE)
