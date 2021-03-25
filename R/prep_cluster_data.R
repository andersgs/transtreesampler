#' Prep cluster data
#'
#' @param df
#' @param onset_date_col
#' @param sampling_date_col
#' @param max_samples
#'
#' @return
#' @export
#'
#' @examples
#'
prep_cluster_data <- function(df,
                              onset_date_col,
                              sampling_date_col,
                              max_samples=NULL) {
  df %<>%
    dplyr::mutate(sample_index = 1:nrow(.)) %>%
    dplyr::mutate(onset_date_days = as.integer(!!onset_date_col - min(!!sampling_date_col)))

  if (!is.null(max_samples)) {
    df %<>%
      dplyr::slice(1:max_samples)
  }
  return(df)
}
