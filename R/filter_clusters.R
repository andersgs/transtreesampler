#' Filter clusters
#'
#' @param df
#' @param cluster_id
#' @param cluster_col
#'
#' @return
#' @export
#'
#' @examples
filter_clusters <- function(df, cluster_id, cluster_col, min_size) {
  df %<>%
    dplyr::group_by(!!cluster_col) %>%
    # nest to make it easy to apply elements to each cluster
    tidyr::nest() %>%
    # calculate the number of samples per cluster
    dplyr::mutate(size = purrr::map_int(.x = data, ~nrow(.x))) %>%
    # arrange by size --- an aesthetic option
    dplyr::arrange(desc(size)) %>%
    # remove clusters that do not meet the minimum requirement
    dplyr::filter(size >= min_size)
  if (cluster_id != '') {
    df %<>%
      dplyr::filter(!!cluster_col == cluster_id)
  }
  if (nrow(df) == 0) {
    stop("No clusters meet the minimum requirements")
  }
  return(df)
}
