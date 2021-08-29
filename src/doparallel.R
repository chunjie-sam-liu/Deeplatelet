library(assertthat)
# Start
fn_parallel_start <- function(n_cores = 10) {
  # https://github.com/ljdursi/beyond-single-core-R
  # detailed information for R parallel computing.
  n_detected_cores <- parallel::detectCores()

  assert_that(n_detected_cores > n_cores, msg = glue::glue("The n_core {n_cores} exceeds the total number of cores {n_detected_cores}."))
  print(glue::glue("Notice: Running process in parallel with {n_cores} cores."))

  global_cluster <<- parallel::makeForkCluster(nnodes = n_cores)

  doParallel::registerDoParallel(cl = global_cluster)
  print(glue::glue("Notice: {foreach::getDoParWorkers()} processes in background."))
}


# Stop
fn_parallel_stop <- function() {
  assert_that(not_empty(showConnections()), msg = "No cluster connections.")
  print(glue::glue("Notice: Stop cluster."))
  parallel::stopCluster(cl = global_cluster)
  foreach::registerDoSEQ()
  print(glue::glue("Notice: {foreach::getDoParWorkers()} workers in background."))
}
