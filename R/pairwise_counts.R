#' Calculate frequency of changes for every sequence pair
#'
#' @param pair
#'
#' @description This counts the frequency of each pattern in a pair of aligned DNA
#' sequences.
#'
#' @export
#'
count_pat <- function(pair) {
  # given the indexes for the two sequences in the alignment
  # subset them out
  # print(c(s1, s2, str(aln)))
  # subs <- phangorn:::subset.phyDat(aln, subset = c(s1, s2))
  # # notify the user what pair is being worked on
  print(names(pair))
  # # to use the table function below to tally all the valid patterns, we
  # # first need to create a vector in which each entry has two characters
  # # the input is a single column of alignment.
  # # for instance, if the column of the alignment is "a", "a", then return "aa"
  # # to ensure that "c" "g" and "g", "c" are counted as the same pattern, we
  # # first lexically sort the column of the alignment.
  # #
  # # we do this because we are implementing time reversible models, so the direction
  # # of substitution doesn't matter, So, for instance, the prob of a -> c is the same as c -> a
  # # in the future, if we would wish to account for directionality, we could modify this
  # # function to distinguish between 'ac', and 'ca'.
  chars <- phangorn:::as.character.phyDat(pair)
  patterns <- apply(chars, 2, function(x) paste(sort(x), collapse = ''))
  # # we then transform the vector above into a factor with the appropriately valid
  # # levels --- these will exclude gaps, Ns, indels, and extended IUPAC codes.
  count <- table(factor(patterns,
                        levels = c('aa', 'cc', 'gg', 'tt',
                                   'ac', 'ag', 'at', 'cg', 'ct', 'gt')))
  return(as.integer(count))
}


#' Pairwise pattern count
#'
#' @param aln
#' @param cores
#'
#' @description This is a wrapper for count_pat to apply it to all possible pairs.
#' It returns a tibble with three columns:
#'  * s1 - index of the first sample in the comparison
#'  * s2 - index of the second sample in the comparison
#'  * counts - the counts of all valid observed nucleotide pairs (e.g., aa, ac, ag, etc..)
#'
#'  The mapping to pairs is parallelised with furrr
#'  @export
pw_pattern_count <- function(aln, cores=1, strategy=future::multisession) {
  future::plan(strategy, workers = cores)
  n_samps <-length(aln)
  df <- t(combn(n_samps, 2)) %>%
    tibble::as_tibble(.name_repair = function(x) paste("s", 1:length(x), sep="")) %>%
    dplyr::mutate(pairs = purrr::map2(.x = s1, .y = s2,
                                      ~phangorn:::subset.phyDat(aln, subset = c(.x, .y)))) %>%
    dplyr::mutate(counts = furrr::future_map(.x = pairs, ~transtreesampler::count_pat(.x)))
  return(df)
}
