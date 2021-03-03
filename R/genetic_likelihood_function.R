#' Extend the genetic likelihood function for Outbreaker2
#'
#' @param rate_params
#' @param base_freq
#'
#' @description Build a transition matrix to be used
#' to calculate the prob matrix
#'
#' @export
#'
build_Q_mat <- function(rate_params, base_freq) {
  # rate parms needs to be in order:
  # a -> c
  # a -> g
  # c -> g
  # a -> t
  # a -> c
  # a -> t
  # base freq needs to be in order:
  # a, c, g, t
  n_states <- length(base_freq)
  qmat <- matrix(0, n_states, n_states)
  qmat[upper.tri(qmat)] <- rate_params
  qmat <- qmat * base_freq
  qmat <- qmat + t(qmat)
  diag(qmat) <- -1 * rowSums(qmat)
  if(!isSymmetric(qmat)) {
    stop("Matrix is not Symmetric, something is not quite right.")
  }
  stat <- ape::matexpo(qmat*10000)
  is_stattionary <- all.equal(stat, matrix(rep(0.25, n_states^2), nrow=4))
  if(!is_stattionary) {
    stop("Matrix is not stationary. Maybe something wrong with rate params. Are they in the right order?")
  }
  return(qmat)
}

#' Count lookup
#'
#' @param i
#' @param j
#' @param df
#'
#' @description A helper function to lookup counts during likelihood calculations.
#'
#' @return
#' @export
#'
#' @examples
count_lookup <- function(i, j, df=counts_df) {
  ix <- sort(i, j)
  return(counts_df %>%
           dplyr::filter(s1 == ix[1] & s2 == x[2]) %>%
           dplyr::pull(counts))
}

#' A generic custom genetic likelihood function
#'
#' @param data
#' @param param
#' @param i
#' @param qmat
#' @param pat_counts
#'
#' @return
#' @export
#'
#' @examples
genetic_ll <- function(data, param, i=NULL, qmat=NULL, pat_counts=NULL) {

  if(is.null(i)) {
    i <- seq_len(data$N)
  }

  i <- i[!is.na(param$alpha[i])]

  if(length(i) == 0) {
    return(0.0)
  }

  mu <- param$mu
  log_like <- 0.0
  for(this_i in i) {
    anc <- param$alpha[this_i]
    ix <- sort(c(this_i, anc))
    counts <- pat_counts %>%
      dplyr::filter(s1 == ix[1] & s2 == ix[2]) %>%
      dplyr::pull(counts) %>%
      magrittr::extract2(1)
    # question: are we underestimating the number of days here?
    # if we assume the MRCA of all the virions in a patient is
    # contained within the patient, and acquired during the infection
    # event bottleneck or shortly there after.
    # then, but counting only the days between the infection times
    # we are only accounting for the evolutionary time along the
    # lineage from infection. but the actual time between the
    # two observed sequences is going to be the t_inf to sample
    # time in the primary patient, and t_inf in p1 to t_inf in p2 plus
    # time to collection time in the secondary patient?
    # 2021-01-04: to address the issue above, I have changed the function below
    # to reflect that evolution time between the two sequences is likely to be
    # from the time of infection of the ancestor to the time of sampling in the
    # ancestor plus from time of infection in the ancestor to the time of
    # sampling in the infectee. --- This assumes that the infection in the
    # ancestor generates a bottleneck event that constrains the sequence to a
    # single variant
    days_between_cases = abs((data$dates[this_i] - param$t_inf[anc]) + (data$dates[anc] - param$t_inf[anc]))
    pmat <- ape::matexpo(qmat*param$mu*days_between_cases)
    log_p <- c(log(diag(pmat)), log(pmat[upper.tri(pmat)]))
    log_like <-  log_like + sum(log_p * counts)
  }

  return(log_like)

}

# run some tests
# uncomment to run
# fA = 0.300
# fC = 0.183
# fG = 0.196
# fT = 0.321
#
# base_f = c(fA, fC, fG, fT)
#
# ac = 0.245
# ag = 0.541
# at = 0.007
# cg = 0.273
# ct = 4.748
# gt = 1.000
#
# rate_p <- c(ac, ag, cg, at, ct, gt)
#
# (tmp <- build_Q_mat(rate_p, base_f))
