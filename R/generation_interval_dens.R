
#' Generate distribution of generation times
#'
#' @description  We use a generation time distribution modeled as a discretized shifted gamma
#' distribution, with a mean of mu days, a standard deviation of sigma days, and shift
#' 1, using the function discr_si from the R package EpiEstim. We set the colonization
#' time distribution equal to the generation time distribution. This is the same
#' approach used in 'Bayesian Reconstruction of Disease Outbreaks by Combining
#' Epidemiologic and Genomic Data', Thibaut Jombart et al (original Outbreaker paper).
#'
#' @param mu
#' @param stdev
#' @param len
#'
#' @return
#' @export
#'
#' @examples
get_generation_internal_dens <- function(mu, stdev, len=20) {
  x <- seq(1, len)
  dens <- EpiEstim::discr_si(x, mu = mu, sigma = stdev)
  return(list(
    w_dens = dens,
    f_dens = dens
  ))
}
