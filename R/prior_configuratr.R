#' Prior configurator
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
prior_configuratr <- function(...) {
  config_params = list(
    init_mu = 1e-3,
    init_pi = 0.6,
    find_import = FALSE,
    n_iter = 1*10^3,
    move_alpha = TRUE,
    move_t_inf = TRUE,
    move_mu = TRUE,
    move_kappa = TRUE,
    move_eps = FALSE,
    move_lambda = FALSE,
    move_pi = TRUE,
    sample_every = 1,
    prior_mu = 3.06e-6,
    sd_mu = 8e-7,
    sd_pi = 0.05,
    prop_alpha_move = 0.1,
    prop_t_inf_move = 0.1,
    prior_pi = c(3, 2),
    max_kappa = 5,
    min_date = -28,
    pb = TRUE
  )
  args <- list(...)
  for (arg in names(args)) {
    if(exists(arg, where=config_params)) {
      config_params[[arg]] <- args[[arg]]
    }
  }
  return(config_params)
}
