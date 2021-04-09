#' Title
#'
#' @param data
#' @param aln
#' @param threshold
#' @param onset_date_col
#' @param model
#' @param cluster_col
#' @param min_cluster_size
#'
#' @return
#' @export
#'
#' @examples
id_transmission_pairs <- function(data,
                                  aln,
                                  max_genetic_dist,
                                  min_day_diff,
                                  onset_date_col,
                                  case_id_col,
                                  model='TN93',
                                  max_day_diff = Inf,
                                  cluster_col = NULL,
                                  min_cluster_size=0) {
  if(!is.null(cluster_col)) {
    data %>%
      dplyr::group_by(!!cluster_col) %>%
      tidyr::nest() %>%
      dplyr::filter(purrr::map_lgl(data, .f=~nrow(.x)>=min_cluster_size)) %>%
      dplyr::mutate(genetic_dists = purrr::map(data,
                                               .f=~id_transmission_pairs(data = .x,
                                                                    aln = aln,
                                                                    max_genetic_dist=max_genetic_dist,
                                                                    min_day_diff = min_day_diff,
                                                                    max_day_diff = max_day_diff,
                                                                    onset_date_col = onset_date_col,
                                                                    case_id_col = case_id_col,
                                                                    model = model,
                                                                    cluster_col = NULL,
                                                                    min_cluster_size = 0)))

  } else {
    tmp <- data %>%
      dplyr::select(!!case_id_col, !!onset_date_col)

    data %>%
      dplyr::arrange(!!onset_date_col) %>%
      dplyr::pull(!!case_id_col) %>%
      (function(ids) {
        ape::dist.dna(aln[ids], model = model)
      }) %>%
      broom::tidy() %>%
      dplyr::rename(case_i = item1, case_j = item2) %>%
      dplyr::filter(distance <= max_genetic_dist) %>%
      dplyr::left_join(tmp, c(case_i = rlang::as_name(case_id_col))) %>%
      dplyr::rename(onset_i = !!onset_date_col) %>%
      dplyr::left_join(tmp, c(case_j = rlang::as_name(case_id_col))) %>%
      dplyr::rename(onset_j = !!onset_date_col) %>%
      dplyr::mutate(onset_diff = onset_j - onset_i) %>%
      dplyr::filter(onset_diff >= min_day_diff & onset_diff <= max_day_diff) %>%
      I()
  }
}
