
# Library -----------------------------------------------------------------

library(magrittr)
library(survival)
library(survminer)
library(survivalROC)

# src ---------------------------------------------------------------------

source(file='src/doparallel.R', local = TRUE)


# Load data ---------------------------------------------------------------
total416.df <- feather::read_feather(path = 'data/rda/total416.os.se.norm.coxph.feather')

# Function ----------------------------------------------------------------
fn_filter_data <- function(.oc, .data = total416.df) {
  .data %>%
    dplyr::filter(oc == .oc) %>%
    dplyr::select(-c(barcode, oc))
}

fn_as_formula <- function(.factors) {
  as.formula(glue::glue('Surv(duration, event) ~ {paste0(.factors, collapse = "+")}'))
}

fn_coxph <- function(.data) {
  .colnames <- colnames(.data)
  .factors <- .colnames[1:(length(.colnames)-2)]
  .formula <- fn_as_formula(.factors = .factors)
  coxph(formula = .formula, data = .data )
}

fn_survivalROC_helper <- function(.x, .d) {
  survivalROC(
    Stime = .d$duration,
    status = .d$event,
    marker = .d$lp,
    predict.time = .x,
    method = 'NNE',
    span = 0.25 * nrow(.d)^(-0.2)
  )
}

# Filter data -------------------------------------------------------------

train <- fn_filter_data(.oc = 'OC521')
eval <- fn_filter_data(.oc = 'OC44')
test1 <- fn_filter_data(.oc = 'OC79')
test2 <- fn_filter_data(.oc = 'OC172')


# Cox model ---------------------------------------------------------------
model <- fn_coxph(.data = train)
train_new <- train
train_new$lp <- predict(model, newdata = train)
train_new$risk <- ifelse(train_new$lp > median(train_new$lp), 'high', 'low')

ggsurvplot(
  fit = survfit(formula = Surv(duration, event) ~ risk, data = train_new),
  data = train_new,
  pval = T,
  title = 'TC data high low risk'
)

survConcordance(Surv(duration, event) ~ lp, data = train_new)

train_roc <- tibble::tibble(
  t = 12 * c(1, 2, 3, 4, 5, 6, 7, 8, 9)
) %>%
  dplyr::mutate(
    survivalROC = purrr::map(.x = t, .f = fn_survivalROC_helper, .d = train_new),
    auc = purrr::map_dbl(survivalROC, magrittr::extract2, 'AUC'),
    df_survivalROC = purrr::map(survivalROC, function(obj) {
    tibble::as_tibble(obj[c('cut.values', 'TP', 'FP')])
  })
) %>%
  dplyr::select(-survivalROC) %>%
  tidyr::unnest(df_survivalROC) %>%
  dplyr::arrange(t, FP, TP)

train_roc %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  geom_label(
    data = train_roc %>% dplyr::select(t, auc) %>% unique,
    mapping = aes(label =sprintf("%.3f", auc), x = 0.5, y = 0.5 )
  ) +
  facet_wrap(~t) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank()
  ) +
  labs(title = 'TC data AUC by the years')

eval_new <- eval
eval_new$lp <- predict(model, newdata = eval)
eval_new$risk <- ifelse(eval_new$lp > median(eval_new$lp), 'high', 'low')
ggsurvplot(
  fit = survfit(formula = Surv(duration, event) ~ risk, data = eval_new),
  data = eval_new,
  pval = T,
  title = 'DC data high low risk'
)

survConcordance(Surv(duration, event) ~ lp, data = eval_new)

survivalROC(
  Stime = eval_new$duration,
  status = eval_new$event,
  marker = eval_new$lp,
  predict.time = 60,
  method = 'NNE',
  span = 0.25 * nrow(eval_new)^(-0.2)
)

eval_roc <- tibble::tibble(
  t = 12 * c(1, 2, 3, 4, 5, 6, 7, 8, 9)
) %>%
  dplyr::mutate(
    survivalROC = purrr::map(.x = t, .f = fn_survivalROC_helper, .d = eval_new),
    auc = purrr::map_dbl(survivalROC, magrittr::extract2, 'AUC'),
    df_survivalROC = purrr::map(survivalROC, function(obj) {
      tibble::as_tibble(obj[c('cut.values', 'TP', 'FP')])
    })
  ) %>%
  dplyr::select(-survivalROC) %>%
  tidyr::unnest(df_survivalROC) %>%
  dplyr::arrange(t, FP, TP)

eval_roc %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  geom_label(
    data = eval_roc %>% dplyr::select(t, auc) %>% unique,
    mapping = aes(label =sprintf("%.3f", auc), x = 0.5, y = 0.5 )
  ) +
  facet_wrap(~t) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank()
  ) +
  labs(title = 'DC data AUC by the years')


test1_new <- test1
test1_new$lp <- predict(model, newdata = test1)
test1_new$risk <- ifelse(test1_new$lp > median(test1_new$lp), 'high', 'low')
ggsurvplot(
  fit = survfit(formula = Surv(duration, event) ~ risk, data = test1_new),
  data = test1_new,
  pval = T,
  title = 'VC1 data high low risk'
)

survConcordance(Surv(duration, event) ~ lp, data = test1_new)



test1_roc <- tibble::tibble(
  t = 12 * c(1, 2, 3, 4, 5, 6, 7, 8, 9)
) %>%
  dplyr::mutate(
    survivalROC = purrr::map(.x = t, .f = fn_survivalROC_helper, .d = test1_new),
    auc = purrr::map_dbl(survivalROC, magrittr::extract2, 'AUC'),
    df_survivalROC = purrr::map(survivalROC, function(obj) {
      tibble::as_tibble(obj[c('cut.values', 'TP', 'FP')])
    })
  ) %>%
  dplyr::select(-survivalROC) %>%
  tidyr::unnest(df_survivalROC) %>%
  dplyr::arrange(t, FP, TP)

test1_roc %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  geom_label(
    data = test1_roc %>% dplyr::select(t, auc) %>% unique,
    mapping = aes(label =sprintf("%.3f", auc), x = 0.5, y = 0.5 )
  ) +
  facet_wrap(~t) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank()
  ) +
  labs(title = 'VC1 data AUC by the years')



test2_new <- test2
test2_new$lp <- predict(model, newdata = test2)
test2_new$risk <- ifelse(test2_new$lp > median(train_new$lp), 'high', 'low')
ggsurvplot(
  fit = survfit(formula = Surv(duration, event) ~ risk, data = test2_new),
  data = test2_new,
  pval = T,
  title = 'VC2 data high low risk'
)

survConcordance(Surv(duration, event) ~ lp, data = test2_new)


test2_roc <- tibble::tibble(
  t = 12 * c(1, 2, 3, 4, 5, 6, 7, 8, 9)
) %>%
  dplyr::mutate(
    survivalROC = purrr::map(.x = t, .f = fn_survivalROC_helper, .d = test2_new),
    auc = purrr::map_dbl(survivalROC, magrittr::extract2, 'AUC'),
    df_survivalROC = purrr::map(survivalROC, function(obj) {
      tibble::as_tibble(obj[c('cut.values', 'TP', 'FP')])
    })
  ) %>%
  dplyr::select(-survivalROC) %>%
  tidyr::unnest(df_survivalROC) %>%
  dplyr::arrange(t, FP, TP)

test2_roc %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  geom_label(
    data = test2_roc %>% dplyr::select(t, auc) %>% unique,
    mapping = aes(label =sprintf("%.3f", auc), x = 0.5, y = 0.5 )
  ) +
  facet_wrap(~t) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank()
  ) +
  labs(title = 'VC2 data AUC by the years')


# Save image --------------------------------------------------------------
save.image(file = 'data/rda/074-coxph-modeling-os.rda')

