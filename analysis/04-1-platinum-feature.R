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

fn_se2data.frame <- function(.se) {
  .df <- as.data.frame(t(assay(.se)))
  .target <- .se@colData[, 'platinum']
  .target <- ordered(x = .target, levels = c("sensitive", "resistant"))
  
  .df.t <- cbind(platinum = .target, .df)
  
  .df.t
}

fn_sd_median <- function(.df) {
  .sd <- apply(
    X = .df[, -1],
    MARGIN = 2,
    FUN = sd
  )
  .median <- apply(
    X = .df[, -1],
    MARGIN = 2,
    FUN = median
  )
  intersect(
    names(.sd[.sd > 0]),
    names(.median[.median > 0])
  ) -> .feats.sd.median
  .feats.sd.median
}


fn_select_features <- function(.se, .train_samples) {
  .train.se <- .se[, .se$oc == .train_samples]
  
  .df <- fn_se2data.frame(.se = .train.se)
  
  .feats.sd.median <- fn_sd_median(.df = .df)
  
  set.seed(123)
  .elastic <-  caret::train(
    platinum ~ .,
    data = .df[, c('platinum', .feats.sd.median)],
    method = "glmnet",
    family = "binomial",
    trControl = caret::trainControl(method = "cv", number = 10),
    tuneLength = 10
  )
  
  .elastic.dgcmatrix <- coef(.elastic$finalModel, .elastic$bestTune$lambda)
  
  as.matrix(.elastic.dgcmatrix) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "Gene_ID") %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(abs(`1`) > 0, grepl(pattern = 'ENSG', x = Gene_ID)) %>%
    dplyr::arrange(-abs(`1`)) %>% 
    dplyr::pull(Gene_ID) ->
    .feats.regr
  
  readr::write_rds(x = .feats.regr, file = "data/rda/00-regr-feature.rds", compress = "gz")
  
  .df.regr <- .df[, c("platinum", .feats.regr)]
  .df.regr.mrmr.dd <- mRMRe::mRMR.data(data = .df.regr)
  .df.regr.mrmr.filter <- mRMRe::mRMR.ensemble(
    data = .df.regr.mrmr.dd,
    target_indices = c(1),
    solution_count = 1,
    feature_count = length(.feats.regr)
  )
  .feats.mrmr <- .df.regr.mrmr.filter@feature_names[.df.regr.mrmr.filter@filters$`1`[, 1]]
  .feats.mrmr_matrix <- cbind(.df.regr.mrmr.filter@filters$`1`, .df.regr.mrmr.filter@scores$`1`)
  rownames(.feats.mrmr_matrix) <- .feats.mrmr
  colnames(.feats.mrmr_matrix) <- c("order", "score")
  .feats.mrmr_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'Gene_ID') %>%
    tibble::as_tibble() %>%
    dplyr::arrange(-abs(score)) %>%
    dplyr::pull(Gene_ID) ->
    .feats.mrmr_sort
  
  
  set.seed(123)
  .train_ind <- sample(ncol(.train.se), size = ncol(.train.se) * 0.7)
  .test_ind <- setdiff(seq_len(ncol(.train.se)), .train_ind)
  .train.feats.se <- .train.se[.feats.mrmr_sort, .train_ind]
  .test.feats.se <- .train.se[.feats.mrmr_sort, .test_ind]
  
  .ifs.res <- foreach(
    i = 50:length(.feats.mrmr_sort),
    .combine = rbind,
    .packages = c('magrittr', 'DESeq2', 'ggplot2'),
    .export = c('fn_ifs_eval', 'fn_caret_train_model', 'fn_ensembl_predict', 'fn_predict', 'fn_getROC_AUC', 'fn_alist', 'fn_plot_auc', 'fn_se2data.frame')
  ) %dopar% {
    fn_ifs_eval(.x = i, .feats = .feats.mrmr_sort, .train.se = .train.feats.se, .eval.se = .test.feats.se)
  }
  
  readr::write_rds(.ifs.res, file = "data/rda/platinum-ifs.rds", compress = 'gz')
  
  ifs.res.plot <- fn_ifs_plot(.ifs.res = .ifs.res)
  
  ggsave(
    filename = '01-Training-evaluation-on-selected-features.pdf',
    plot = ifs.res.plot,
    device = 'pdf',
    path = 'data/output',
    width = 8,
    height = 4
  )
  
  .optimal <- .ifs.res %>%
    dplyr::arrange(-eval, num)
  
  .optimal_feats_num <- .optimal %>% dplyr::pull(num) %>% head(1)
  .optimal_parameter <- .optimal %>% dplyr::pull(best_tune) %>% head(1) %>% .[[1]]
  
  list(
    panel = head(.feats.mrmr_sort, .optimal_feats_num),
    hyperparameter = .optimal_parameter
  )
}



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
fn_caret_train_model <- function(.se, method = 'svmRadial') {
  .d <- fn_se2data.frame(.se = .se)
  .model <- caret::train(
    platinum ~ .,
    data = .d,
    method = method,
    # preProcess = c("center", "scale", "YeoJohnson"),
    metric = 'ROC',
    maximize = TRUE,
    trControl = caret::trainControl(
      method = "repeatedcv", number = 10,
      classProbs = TRUE, summaryFunction = twoClassSummary
    ),
    tuneGrid = NULL,
    tuneLength = 10,
    verbose = FALSE
  )
  .model
}
fn_ensembl_predict <- function(.model, .se, .title = 'OC521') {
  .d.pred <- fn_predict(.model = .model, .se = .se)
  .d.pred.alist <- fn_alist(.c = .se$platinum, .pred = .d.pred)
  .d.p <- fn_plot_auc(.alist = .d.pred.alist, .t = .title)
  list(predict = .d.pred, alist = .d.pred.alist, plot = .d.p)
}
# Feature selection -------------------------------------------------------

# start parallel
fn_parallel_start(n_cores = 95)
# feature selection

feats <- fn_select_features(
  .se = total351.platinum.se.norm,
  .train_samples = 'OC521'
  )
fn_parallel_stop()

readr::write_rds(x = feats, file = 'data/rda/00-selected-features.rds.gz', compress = 'gz')

# Save image --------------------------------------------------------------

save.image(file = 'data/rda/04-1-platinum-feature.rda', ascii = FALSE, compress = TRUE)
