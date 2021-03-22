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

fn_se2df <- function(.se) {
  .df <- as.data.frame(t(assay(.se)))
  .target <- .se$platinum
  .df$platinum <- .target
  .df
}

fn_sd_median <- function(.se) {
  .sd <- assay(.se) %>%
    apply(MARGIN = 1, FUN = sd)
  .md <- assay(.se) %>%
    apply(MARGIN = 1, FUN = median)
  
  names(.md[.md > 0])
}
fn_select_features <- function(.se) {
  .train.se <- .se[, .se$oc == "OC521"]
  
  .feats1 <- fn_sd_median(.se = .train.se)
  
  .df <- fn_se2df(.se = .train.se[.feats1, ])
  
  set.seed(123)

  fn_parallel_start(n_cores = 50)
  
  .glm <-  caret::train(
    platinum ~ .,
    data = .df,
    method = "glmnet",
    family = "binomial",
    trControl = caret::trainControl(method = "cv", number = 5),
    tuneGrid =  expand.grid(
      alpha = 1,
      lambda = 10^seq(-3, 3,length = 100)
    ),
    tuneLength = 5
  )
  
  fn_parallel_stop()
  
  coef(.glm$finalModel, .glm$bestTune$lambda) %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'ensid') %>%
    dplyr::filter(`1` != 0) %>%
    dplyr::filter(grepl(pattern = "ENSG", x = ensid)) %>%
    dplyr::pull(ensid)
  
}

# Feature selection -------------------------------------------------------

panel <- fn_select_features(.se = total351.platinum.se.norm)
readr::write_rds(x = panel, file = 'data/rda/panel.rds.gz', compress = 'gz')

# Save image --------------------------------------------------------------

save.image(file = 'data/rda/04-platinum-feature.rda', ascii = FALSE, compress = TRUE)
