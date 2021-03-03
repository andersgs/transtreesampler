#' The likelihood function factory
#'
#' @description This encapsulates all the data the likelihood function needs
#' to run into a single environment, and still fits
#' into a slightly modified outbreaker2 which allows
#' for i to be passed into the custom genetic
#' likelihood function.
#' this returns a function that can be consumed by outbreaker2
#'
#' @param ll_func
#' @param build_q_mat
#' @param aln
#' @param rate_param
#' @param base_freq
#'
#' @return
#' @export
#'
#' @examples
ll_factory <- function(ll_func,
                       build_q_mat,
                       aln,
                       rate_param,
                       base_freq,
                       cores=1,
                       strategy=future::sequential) {
  pat_counts <- transtreesampler::pw_pattern_count(aln, cores = cores)
  qmat <- build_q_mat(rate_param, base_freq)
  genetic_ll <- function(data, param, i = NULL) {
    return(ll_func(data,
                   param,
                   i,
                   qmat=qmat,
                   pat_counts = pat_counts))
  }
  return(genetic_ll)
}
