#' Parse metadata
#'
#' @param metadata_file
#' @param data_dir
#' @param cluster_id
#' @param cluster_col
#' @param sample_col
#' @param aln
#'
#' @return
#' @export
#'
#' @examples
prep_data <- function(metadata_file,
                           data_dir,
                           cluster_id,
                           aln,
                           cluster_col=nat_comm_clusters,
                           sample_col=vic_id) {
  meta <- readr::read_csv(file.path(data_dir, metadata),
                          col_types = readr::cols()) %>%
    dplyr::filter({{ cluster_col }} == cluster_id) %>%
    dplyr::arrange({{ sample_col }})
  # pull out some relevant data points
  onset_date <- meta %>%
    dplyr::pull(onset_date) %>%
    as.Date()
  # sample date
  sample_date <- meta %>%
    dplyr::pull(diagnosis_date) %>%
    as.Date()
  # sample ids
  sample_ids <- meta %>%
    dplyr::pull({{ sample_col }})
  #load align
  aln <- phangorn::read.phyDat(file.path(data_dir, fasta), format = "fasta")
  aln <- aln[sample_ids]
  data <- list(
    meta = meta,
    onset_date = onset_date,
    dates = sample_date,
    ids = sample_ids,
    aln = aln
  )
  return(data)
}
