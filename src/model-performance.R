
# Train and predict
fn_resample <- function(.se, .id) {
  .task <- fn_se2task(.se = .se, .id = .id)
  ksvm.lrn <- makeLearner(cl = 'classif.ksvm', predict.type = 'prob')
  resamdesc.inner <- makeResampleDesc(method = 'CV', iters = 10, stratify = TRUE, predict = 'both')
  resample(
    learner = ksvm.lrn,
    task = .task,
    resampling = resamdesc.inner,
    measures = mlr::auc,
    models = TRUE,
    keep.pred = TRUE
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
fn_predict <- function(.model, .se) {
  .d <- fn_se2data.frame(.se = .se)
  .p <- predict(.model, .d, type = "prob")
  rownames(.p) <- rownames(.d)
  .p$correct <- ifelse(.p$sensitive > .p$resistant, 'sensitive', 'resistant') == as.character(.d$platinum)
  .p
}
fn_getROC_AUC <- function(probs, true_Y) {
  probsSort <- sort(probs, decreasing = TRUE, index.return = TRUE)
  val <- unlist(probsSort$x)
  idx <- unlist(probsSort$ix)

  roc_y <- true_Y[idx]
  stack_x <- cumsum(roc_y == 2) / sum(roc_y == 2)
  stack_y <- cumsum(roc_y == 1) / sum(roc_y == 1)

  auc <- sum((stack_x[2:length(roc_y)] - stack_x[1:length(roc_y) - 1]) * stack_y[2:length(roc_y)])
  return(list(stack_x = stack_x, stack_y = stack_y, auc = auc))
}
fn_alist <- function(.c, .pred) {
  .t <- as.numeric(.c)
  aList <- fn_getROC_AUC(probs = .pred[[1]], true_Y = .t)
}
fn_plot_auc <- function(.alist, .t) {
  tibble::tibble(
    x = .alist$stack_x,
    y = .alist$stack_y
  ) %>%
    ggplot(aes(x = x, y = y)) +
    geom_path() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = 0.05, linetype = 2) +
    theme_bw() +
    labs(
      x = 'False Positive Rate',
      y = 'True Positive Rate',
      title = glue::glue("{.t}, AUC = {round(.alist$auc,3)}")
    )
}
fn_ensembl_predict <- function(.model, .se, .title = 'OC521') {
  .d.pred <- fn_predict(.model = .model, .se = .se)
  .d.pred.alist <- fn_alist(.c = .se$platinum, .pred = .d.pred)
  .d.p <- fn_plot_auc(.alist = .d.pred.alist, .t = .title)
  list(predict = .d.pred, alist = .d.pred.alist, plot = .d.p)
}


# Performance -------------------------------------------------------------

fn_roc_95ci <- function(.p) {
  # acc - Accuracy
  # tpr - True positive rate (Sensitivity, Recall)
  # tnr - True negative rate (Specificity)
  # ppv - Positive predictive value (Precision)
  # npv - Negative predictive value
  # kappa - Cohen's kappa
  # f1 - F1 measure
  # mlr::performance(pred = .pred$panel, measures = list(mlr::acc, mlr::tpr, mlr::tnr, mlr::ppv, mlr::npv, mlr::kappa, mlr::f1))
  .rocm <- mlr::calculateROCMeasures(.p)
  .kappa_f1 <- mlr::performance(pred = .p, measures = list(mlr::kappa, mlr::f1))
  .rocm_epi <- epiR::epi.tests(dat = t(.rocm$confusion.matrix))
  tibble::tibble(
    acc = list(unlist(.rocm_epi$rval$diag.acc)),
    tpr = list(unlist(.rocm_epi$rval$se)),
    tnr = list(unlist(.rocm_epi$rval$sp)),
    ppv = list(unlist(.rocm_epi$rval$ppv)),
    npv = list(unlist(.rocm_epi$rval$npv)),
    kappa = .kappa_f1['kappa'],
    f1 = .kappa_f1['f1']
  )

}

fn_performance <- function(.model_tuned, .task, .samples) {
  purrr::map(.x = .samples, .f = function(.x) {
    predict(object = .model_tuned, task = .task, subset = .x) %>%
      mlr::performance(measures = mlr::auc)
  }) %>%
    tibble::enframe() %>%
    tidyr::unnest()
}

fn_ifs_plot <- function(.ifs.res){
  # readr::write_rds(x = .ifs.res, path = file.path(path_analysis_out, '01-ifs-progress.rds.gz'), compress = 'gz')
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


fn_task_inds <- function(.set, .se) {
  .barcode_ind <- which(.se$oc == .set)
  setNames(object = .barcode_ind, nm = .se@colData[.barcode_ind, 'barcode'])
}


# help function for performance
fn_dataset_switch <- function(.dataset, .total_task) {
  list(
    panel = switch (
      .dataset,
      'OC521' = .total_task$samples$panel[[1]],
      'OC44' = .total_task$samples$panel[[2]],
      'OC172' = .total_task$samples$panel[[3]],
      'OC79' = .total_task$samples$panel[[4]]
    ),
    ca125 = switch(
      .dataset,
      'OC521' = .total_task$samples$ca125[[1]],
      'OC44' = .total_task$samples$ca125[[4]],
      'OC172' = .total_task$samples$ca125[[3]],
      'OC79' = .total_task$samples$ca125[[4]]
    ),
    panel.ca125 = switch(
      .dataset,
      'OC521' = .total_task$samples$ca125[[1]],
      'OC44' = .total_task$samples$ca125[[4]],
      'OC172' = .total_task$samples$ca125[[3]],
      'OC79' = .total_task$samples$ca125[[4]]
    )
  )
}
# plot ROC curve
fn_merge_auroc_plot <- function(.merge, .dataset) {
  .merge$curve %>%
    ggplot(aes(x = fpr, y = tpr, color = mod)) +
    geom_path() +
    scale_color_manual(values = c("#006400", "#B22222", "#00008B")) +
    theme_bw() +
    labs(title = glue::glue("{.dataset}, AUROC Platinum sensitive Vs. resistant"), x = "False Positive Rate", y = "True Positive Rate") +
    theme(
      legend.position = "none"
    ) +
    annotate("text", x = 0.5, y = 0.45, label = glue::glue('AUC (CA125): {round(.merge$auc$CA125, 3)}'), color = "#006400", fontface = 2) +
    annotate("text", x = 0.5, y = 0.5, label = glue::glue('AUC (Panel): {round(.merge$auc$Panel, 3)}'), color = "#B22222", fontface = 2) +
    annotate("text", x = 0.5, y = 0.55, label = glue::glue('AUC (Panel+CA125): {round(.merge$auc$`Panel+CA125`, 3)}'), color = "#00008B", fontface = 2)
}
fn_plot_auroc <- function(.pred, .subset, .dataset) {
  .pred_perf <- purrr::map(
    .x = .pred,
    .f = mlr::generateThreshVsPerfData,
    measures = list(mlr::fpr, mlr::tpr, mlr::auc),
    gridsize = length(.subset$panel)
  )
  .merge <- list(
    curve = dplyr::bind_rows(
      .pred_perf$panel.ca125$data %>% dplyr::mutate(mod = 'Panel+CA125'),
      .pred_perf$panel$data %>% dplyr::mutate(mod = 'Panel'),
      .pred_perf$ca125$data %>% dplyr::mutate(mod = 'CA125')
    ),
    auc = list(
      `Panel+CA125` = .pred_perf$panel.ca125$data$auc[1],
      `Panel` = .pred_perf$panel$data$auc[1],
      `CA125` = .pred_perf$ca125$data$auc[1]
    )
  )
  .plot <- fn_merge_auroc_plot(.merge = .merge, .dataset = .dataset)
  # save plot
  ggsave(
    filename = glue::glue('03-AUROC-merge-{.dataset}.pdf'),
    plot = .plot, device = 'pdf',
    path = 'data/output',
    width = 8, height = 5
  )
  return(.plot)
}
# convert confusionmatrix to 95 CI
fn_roc_95ci <- function(.p) {
  # acc - Accuracy
  # tpr - True positive rate (Sensitivity, Recall)
  # tnr - True negative rate (Specificity)
  # ppv - Positive predictive value (Precision)
  # npv - Negative predictive value
  # kappa - Cohen's kappa
  # f1 - F1 measure
  # mlr::performance(pred = .pred$panel, measures = list(mlr::acc, mlr::tpr, mlr::tnr, mlr::ppv, mlr::npv, mlr::kappa, mlr::f1))
  .rocm <- mlr::calculateROCMeasures(.p)
  .kappa_f1 <- mlr::performance(pred = .p, measures = list(mlr::kappa, mlr::f1))
  .rocm_epi <- epiR::epi.tests(dat = t(.rocm$confusion.matrix))
  tibble::tibble(
    acc = list(unlist(.rocm_epi$rval$diag.acc)),
    tpr = list(unlist(.rocm_epi$rval$se)),
    tnr = list(unlist(.rocm_epi$rval$sp)),
    ppv = list(unlist(.rocm_epi$rval$ppv)),
    npv = list(unlist(.rocm_epi$rval$npv)),
    kappa = .kappa_f1['kappa'],
    f1 = .kappa_f1['f1']
  )

}
fn_save_metrics <- function(.pred_metrics, .dataset) {
  .table_metrics <- .pred_metrics %>%
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
    )
  colnames(.table_metrics) <- c('Modality', 'Accuracy (95% CI)', 'SN (95% CI)', 'SP (95% CI)', 'PPV (95% CI)', 'NPV (95% CI)', 'Kappa', 'F1')
  # save metrics table to tsv and excel
  readr::write_tsv(
    x = .table_metrics,
    path = file.path('data/output', glue::glue('04-Metrics-CI95-{.dataset}.tsv'))
  )
  writexl::write_xlsx(
    x = .table_metrics,
    path = file.path('data/output', glue::glue('04-Metrics-CI95-{.dataset}.xlsx'))
  )
  # return
  .table_metrics
}
fn_metrics <- function(.pred, .dataset) {
  # get metrics ci95 with epiR fn_roc_95ci
  .pred_metrics <- purrr::map(.x = .pred, .f = fn_roc_95ci) %>%
    tibble::enframe() %>%
    tidyr::unnest(cols = value) %>%
    dplyr::mutate(name = c('Panel', 'CA125', 'Panel + CA125'))
  # write the table to disk
  .table_metrics <- fn_save_metrics(.pred_metrics = .pred_metrics, .dataset = .dataset)
  # return the tables
  .table_metrics
}
# performance
