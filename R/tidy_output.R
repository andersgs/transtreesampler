#' Tidy outbreaker 2 MCMC output
#'
#' @param mcmc
#'
#' @return
#' @export
#'
#' @examples
tidy_outb2 <- function(mcmc) {
  res_lst <- as.list(mcmc)
  data_tbl <- tibble::as_tibble(res_lst) %>%
    tidyr::pivot_longer(cols=contains("_"),
                        names_to = c(".value", 'sample_index'),
                        names_pattern = "(.*)_(.*)")
  return(data_tbl)
}
