// A NextFlow pipeline to generate Transmission Tree samples

nextflow.enable.dsl=2

Random rnd = new Random()

/*
conditions <- list(
  set1 = list(mu = 5.2, stdev = 1.72),
  set2 = list(mu = 3.95, stdev = 1.51),
  set3 = list(mu = 6.5, stdev = 1.9),
  set4 = list(mu = 5.2, stdev = 4.0)
)
*/

// conditions to test for w_dens and f_dens
conditions = [
  [id: "set1", mu: 5.2, stdev: 1.72],
  [id: "set2", mu: 3.95, stdev: 1.51],
  [id: "set3", mu: 6.5, stdev: 1.9],
  [id: "set4", mu: 5.2, stdev: 4.0]
  ]

process prep_data {
  publishDir "results/data/", mode: "copy"
  cpus 1

  input:
    path metadata
    path fasta

  output:
    path "*.Rdata", emit: rdata

  """
  #!/usr/bin/env Rscript
  # load needed libraries
  library(magrittr)
  # set column name mappers
  cluster_col <- quote($params.cluster_col)
  sample_id_col <- quote($params.sample_id_col)
  sampling_date_col <- quote($params.sampling_date_col)
  onset_date_col <- quote($params.onset_date_col)
  min_cluster_size <- $params.min_cluster_size

  # load data
  aln <- phangorn::read.phyDat("$fasta", format="fasta")
  mdt <- readr::read_csv("$metadata")

  # prepare the data for the likelihood and MCMC
  mdt_prep <- mdt %>%
  # group the data by the cluster ID column
  dplyr::group_by(!!cluster_col) %>%
  # nest to make it easy to apply elements to each cluster
  tidyr::nest() %>%
  # calculate the number of samples per cluster
  dplyr::mutate(size = purrr::map_int(.x = data, ~nrow(.x))) %>%
  # arrange by size --- an aesthetic option
  dplyr::arrange(desc(size)) %>%
  # remove clusters that do not meet the minimum requirement
  dplyr::filter(size >= min_cluster_size) %>%
  # add a column for sample index and add an onset_date_days column
  # to put the date of onset in the same unit as sampling days used by
  # outbreaker2
  dplyr::mutate(data = purrr::map(.x = data,
                                  .f = ~ .x %>%
                                    dplyr::mutate(sample_index = 1:nrow(.)) %>%
                                    #dplyr::slice(1:6) %>%
                                    dplyr::mutate(onset_date_days = as.integer(!!onset_date_col - min(!!sampling_date_col))))) %>%
  # get the sample ids
  dplyr::mutate(sample_ids = purrr::map(.x = data,
                                        ~dplyr::pull(.x,!!sample_id_col))) %>%
  # make sure sampling dates are a Date object
  dplyr::mutate(dates = purrr::map(.x = data,
                                   ~dplyr::pull(.x, !!sampling_date_col) %>%
                                   as.Date(.))) %>%
  # add the DNA column with the the relevant portion of the alignment for each
  # cluster
  dplyr::mutate(dna = purrr::map(.x = sample_ids, ~aln[.x]))

  # split the dataframe by
  splitTTR <- split(mdt_prep , mdt_prep[[cluster_col]])

  for (nm in names(splitTTR)) {
    save_name <- paste0("cluster_", nm, ".Rdata")
    saveRDS(splitTTR[[nm]], file=save_name)
  }
  """
}


process prep_likelihood {

  publishDir "results/likelihood/"
  cpus 8

  input:
    path rdata

  output:
    path "*_dna.Rdata", emit: rdata

  """
  #!/usr/bin/env Rscript
  # load needed libraries
  library(magrittr)
  # specify cores
  n_cores = $task.cpus

  # get parameter values
  base_f <- c($params.fA, $params.fC, $params.fG, $params.fT)
  rate_p <- c($params.ac, $params.ag, $params.cg, $params.at, $params.ct, $params.gt)

  # load the metadata and to the calculations
  ttr <- readRDS("$rdata") %>%
    dplyr::mutate(llfunc = purrr::map(.x = dna,
                                    ~transtreesampler::ll_factory(
    transtreesampler::genetic_ll,
    transtreesampler::build_Q_mat,
    .x,
    rate_p,
    base_f,
    cores = n_cores,
    strategy = future::multisession))) %>%
  dplyr::mutate(this_likelihood = purrr::map(.x = llfunc, ~list(genetic = .x)))
  new_name = stringr::str_replace(string = "$rdata", pattern=".Rdata", replacement="_dna.Rdata")
  saveRDS(ttr, new_name)
  """
}


process set_conditions {

  publishDir "results/conditions"

  input:
    path rdata
    each condition

  output:
    path "*_${condition.id}.Rdata", emit: rdata

  script:
  def condition_id = "$condition.id"
  def condition_mu = "$condition.mu"
  def condition_stdev = "$condition.stdev"
  """
  #!/usr/bin/env Rscript
  # load library
  library(magrittr)

  condition = tibble::tibble(condition="$condition_id",
                              mu=$condition_mu,
                              stdev=$condition_stdev) %>%
  dplyr::mutate(gen_intervals = purrr::map2(.x = mu,
                                            .y = stdev,
                                            .f = ~transtreesampler::get_generation_internal_dens(.x, .y)))

  ttr <- readRDS("$rdata") %>%
              tidyr::expand_grid(condition)
  new_name = stringr::str_replace(string = "$rdata", pattern=".Rdata", replacement="_${condition_id}.Rdata")
  saveRDS(ttr, new_name)
  """
}


process run_mcmc {

  publishDir "results/mcmc"
  cpus 1
  errorStrategy "retry"
  maxRetries 3

  input:
    path rdata

  output:
    path "*_mcmc.Rdata", emit: rdata

  script:
  def paranoid = params.paranoid ? "TRUE" : "FALSE"
  def seed = rnd.nextInt(10000)
  """
  #!/usr/bin/env Rscript
  # load libraries
  library("magrittr")
  # run mcmc
  mcmc <- readRDS("$rdata") %>%
    dplyr::mutate(this_data = purrr::pmap(.l = list(sample_ids, dna, dates, gen_intervals),
                                        .f = function(w, x, y, z) list(ids=w,
                                                                        dna=x %>% phangorn:::as.DNAbin.phyDat(.),
                                                                        dates=y,
                                                                        w_dens=z\$w_dens,
                                                                        f_dens=z\$f_dens
                                                                       ))) %>%
  dplyr::mutate(this_config = list(transtreesampler::prior_configuratr(n_iter = $params.n_iter,
                                                                       sample_every = $params.sample_every,
                                                                       prior_mu = $params.prior_mu,
                                                                       sd_mu = $params.sd_mu,
                                                                       paranoid =$paranoid))) %>%
  dplyr::mutate(ttr = purrr:::pmap(.l=list(this_config, this_likelihood, this_data),
                                  .f=transtreesampler::sample_trees, seed=$seed))
  new_name = stringr::str_replace(string = "$rdata", pattern=".Rdata", replacement="_mcmc.Rdata")
  saveRDS(mcmc, new_name)
  """
}


process tidy_output {

  publishDir "results/tidy_output"

  input:
    path rdata

  output:
    path "*_tidy.Rdata", emit: rdata

  script:
  """
  #!/usr/bin/env Rscript
  #load libraries
  library(magrittr)

  # set column name mappers
  cluster_col <- quote($params.cluster_col)
  sample_id_col <- quote($params.sample_id_col)
  sampling_date_col <- quote($params.sampling_date_col)
  onset_date_col <- quote($params.onset_date_col)

  # tidy the output
  output <- readRDS("$rdata") %>%
    dplyr::mutate(tidy_data = purrr::map(.x=ttr,
                                       .f=transtreesampler::tidy_outb2)) %>%
    dplyr::mutate(tidy_data = purrr::map2(.x = data, .y = tidy_data,
                                        .f = ~ .x %>%
                                          dplyr::left_join(.y, c('sample_index' = 'sample_index')))) %>%
  dplyr::select(nat_comm_clusters, condition, tidy_data) %>%
  tidyr::unnest(c(tidy_data)) %>%
  (function(df) {
    sub_df <- df %>%
      dplyr::select(sample_index, onset_date_days) %>%
      dplyr::distinct() %>%
      dplyr::rename(ancestral_onset_date_days = onset_date_days)
    df %>%
      dplyr::left_join(sub_df,  c('alpha' = 'sample_index'))
  }) %>%
  dplyr::rename_with(~ stringr::str_c("outb2_", .x), .cols = c(step, post, like, prior, mu, pi, eps, lambda, alpha, kappa, t_inf)) %>%
  dplyr::rename_with(~ stringr::str_c("meta_", .x), .cols = c(!!cluster_col, condition, sample_index, !!sample_id_col, !!onset_date_col, !!sampling_date_col, onset_date_days, ancestral_onset_date_days)) %>%
  dplyr::select(dplyr::starts_with("outb2"), dplyr::starts_with("meta"))

    new_name_rds = stringr::str_replace(string = "$rdata", pattern=".Rdata", replacement="_tidy.Rdata")
    saveRDS(output, new_name_rds)
  """
}


process gather_clusters {

  publishDir "results/cluster_data"

  input:
    tuple val(cluster_id), path(rdata)

  output:
    path "*.csv"

  """
  #!/usr/bin/env Rscript
  library(magrittr)

  rdata <- list.files(pattern="*.Rdata")
    dplyr::bind_rows(purrr::map(rdata, readRDS)) %>%
    dplyr::arrange(meta_condition, outb2_step) %>%
    dplyr::rename(indexOfCase = meta_sample_index,
                indexOfAncestor = outb2_alpha,
                numberOfIntermediaries = outb2_kappa,
                symptomOnsetTime = meta_onset_date_days,
                symptomOnsetTimeOfAncestor = meta_ancestral_onset_date_days,
                timeOfInfectionEvent = outb2_t_inf) %>%
    readr::write_csv("cluster_${cluster_id}.csv")

  """
}


workflow {
  prep_data(params.metadata, params.fasta)
  prep_likelihood(prep_data.out.rdata.flatten())
  set_conditions(prep_likelihood.out.rdata, conditions)
  run_mcmc(set_conditions.out.rdata)
  tidy_output(run_mcmc.out.rdata)
  gather_clusters(tidy_output.out.rdata.map(it -> tuple(it.baseName.split("_")[1], it)).groupTuple())
}
