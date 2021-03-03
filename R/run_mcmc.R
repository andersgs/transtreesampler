#' Sample trees

#' @param create_config_params
#' @param likelihood_funcs
#' @param data_params
#' @param seed
#' @param save
#'
#' @return
#' @export
#'
#' @examples
sample_trees <- function(create_config_params,
                         likelihood_funcs,
                         data_params,
                         seed=42,
                         save=FALSE) {

  # build priors
  priors <- do.call(outbreaker2::create_config, create_config_params)
  # priors <- outbreaker2::create_config(
  #   init_mu = 1e-3,
  #   init_pi = 0.6,
  #   find_import = FALSE,
  #   n_iter = 1*10^3,
  #   move_alpha = TRUE,
  #   move_t_inf = TRUE,
  #   move_mu = TRUE,
  #   move_kappa = TRUE,
  #   move_eps = FALSE,
  #   move_lambda = FALSE,
  #   move_pi = TRUE,
  #   sample_every = 1,
  #   sd_mu = sd,
  #   sd_pi = 0.05,
  #   prop_alpha_move = 0.1,
  #   prop_t_inf_move = 0.1,
  #   prior_mu = mu,
  #   prior_pi = c(3, 2),
  #   max_kappa = 5,
  #   min_date = -28,
  #   pb = TRUE
  # )


  model <- do.call(outbreaker2::custom_likelihoods, likelihood_funcs)
  #model <- outbreaker2::custom_likelihoods(genetic = custom_genetic_ll)

  data <- do.call(outbreaker2::outbreaker_data, data_params)
  # data <- outbreaker2::outbreaker_data(
  #   dna = aln %>% ape::as.DNAbin(.),
  #   #dates = onset_date,
  #   dates = sample_date,
  #   ids = sample_ids,
  #   w_dens = w_dens,
  #   f_dens = f_dens
  # )

  set.seed(seed)
  tictoc::tic()
  mcmc_out <- outbreaker2::outbreaker(
    data = data,
    config = priors,
    likelihoods = model
  )
  tictoc::toc()
  return(mcmc_out)
}
