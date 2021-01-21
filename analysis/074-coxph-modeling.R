
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
  pval = T
)

survivalROC(
  Stime = train_new$duration,
  status = train_new$event,
  marker = train_new$lp,
  predict.time = 12,
  method = 'NNE',
  span = 0.25 * nrow(train_new)^(-0.2)
)

eval_new <- eval
eval_new$lp <- predict(model, newdata = eval)
eval_new$risk <- ifelse(eval_new$lp >-10.84, 'high', 'low')
ggsurvplot(
  fit = survfit(formula = Surv(duration, event) ~ risk, data = eval_new),
  data = eval_new,
  pval = T
)

survConcordance(Surv(duration, event) ~ lp, data = eval_new)

survivalROC(
  Stime = eval_new$duration,
  status = eval_new$event,
  marker = eval_new$lp,
  predict.time = 12,
  method = 'NNE',
  span = 0.25 * nrow(eval_new)^(-0.2)
)


test1_new <- test1
test1_new$lp <- predict(model, newdata = test1)
test1_new$risk <- ifelse(test1_new$lp > median(test1_new$lp), 'high', 'low')
ggsurvplot(
  fit = survfit(formula = Surv(duration, event) ~ risk, data = test1_new),
  data = test1_new,
  pval = T
)

survConcordance(Surv(duration, event) ~ lp, data = test1_new)

survivalROC(
  Stime = test1_new$duration,
  status = test1_new$event,
  marker = test1_new$lp,
  predict.time = 12,
  method = 'NNE',
  span = 0.25 * nrow(test1_new)^(-0.2)
)


test2_new <- test2
test2_new$lp <- predict(model, newdata = test2)
test2_new$risk <- ifelse(test2_new$lp > -10.84, 'high', 'low')
ggsurvplot(
  fit = survfit(formula = Surv(duration, event) ~ risk, data = test2_new),
  data = test2_new,
  pval = T
)

survConcordance(Surv(duration, event) ~ lp, data = test2_new)

survivalROC(
  Stime = test2_new$duration,
  status = test2_new$event,
  marker = test2_new$lp,
  predict.time = 30,
  method = 'NNE',
  span = 0.25 * nrow(test2_new)^(-0.2)
)
