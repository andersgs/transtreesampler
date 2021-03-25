# Transmission Tree Sampler

## Background

The goal is to generate posterior samples of transmission trees in order to 
estimate epidemiological parameters of interest.

*WARNING* --- the current implementation is heavily optimized for running SARS-CoV-2
analyses.

## Installing the package

To install this package, you will need `remotes` packages:

```R
install.packages("remotes")
```

Once you have successfully install the `remotes` package, you can run:

```R
remotes::install_github("andersgs/transtreesampler")
```

You can optionally specify the `Ncpus` option to speed things up a bit:

```R
remotes::install_github("andersgs/transtreesampler", Ncpus=8)
```

## Installing NextFlow

The easiest way to run an analysis with this package is to use `NextFlow`. To 
install `NextFlow`, I recommend you use a package manager like `Brew` or `Conda`:

```bash
brew install nextflow
```

```bash
conda install -c bioconda nextflow
```

## Your first run

Once you have the environment with `R` and `NextFlow` installed, as well as the
package, you can get do a test run with the following command (*note*: this 
assumes that you have clusters of samples with at least 15 cases, the default 
minimum number of cases for a cluster to be analyzed --- see below how to change
this). This will run a very short MCMC chain (50 iterations sampled every 1 iteration,
 this is not sufficient for any analyses, just to test that things are working):

```bash
nextflow run andersgs/transtreesampler --fasta /full/path/to/fasta.fa --metadata /full/path/to/metadata.csv
```

The metadata file must have the following columns and data, and must be in CSV
format (this will change in subsequent versions to be more generic --- but see below
if your column headers don't exactly match the headers in the table):

| vic_id | nat_comm_clusters | onset_date | diagnosis_date |
|:------:|:-----------------:|:----------:|:--------------:|
| seq1 | 1 | 2020-02-02 | 2020-02-01 |
| seq2 | 1 | 2020-02-03 | 2020-02-02 |
| seq3 | 2 | 2020-02-04 | 2020-02-03 |
| seq4 | 2 | 2020-02-05 | 2020-02-04 |


The multiFASTA file assumes the sequences are aligned, and have identifiers that 
match the identifiers in the `vic_id` column of the `metadata` file:

```bash
>seq1
ATCG
>seq2
ATCC
>seq3
TTCC
>seq4
TTGC
```

### My column headers don't match those specifed above

That is ok, you can change that by using the following options (change only those
that apply):

```bash
nextflow run andersgs/transtreesampler --fasta /full/path/to/fasta.fa \
                     --metadata /full/path/to/metadata.csv \
                     --cluster_col <my_cluster_id> \
                     --sample_id_col <my_sample_id> \
                     --sampling_date_col <my_sample_date_col> \
                     --onset_date_col <my_date_onset_col>
```

### I want to change the minimum size of clusters

To change the minimum size of clusters to be analyzed (default 15), you need to 
specify the `--min_cluster_size` flag:

```bash
nextflow run andersgs/transtreesampler --fasta /full/path/to/fasta.fa \
                     --metadata /full/path/to/metadata.csv \
                     --min_cluster_size 20
```


### I only want to run a single cluster

If you have metadata with several clusters, but you only want to run a single
one of them through tree sampling process, you can use the `--filter_cluster_id`
option to select the cluster of interest. In the code below, we restrict the 
analysis only to the cluster identified as `67`:

```bash
nextflow run andersgs/transtreesampler --fasta /full/path/to/fasta.fa \
                     --metadata /full/path/to/metadata.csv \
                     --filter_cluster_id 67
```

### I want to run fewer samples per cluster as a test run

Let say you want to run some quick tests, and don't want to run a full set of 
samples in your clusters. To restrict the number of samples included in the 
analysis, you can use the `--max_samples` option. Below, we restrict our 
analysis to the first 6 samples in each cluster:

**NOTE**: This option really should only be used for testing.

```bash
nextflow run andersgs/transtreesampler --fasta /full/path/to/fasta.fa \
                     --metadata /full/path/to/metadata.csv \
                     --max_samples 6
```

### Test run worked, now how do I run a proper length MCMC

To run a proper length MCMC, you need to specify two options that set the 
length of the MCMC and the thinning (how frequently samples are kept):

```bash
nextflow run andersgs/transtreesampler --fasta /full/path/to/fasta.fa \
                     --metadata /full/path/to/metadata.csv \
                     --n_iter 100000 \
                     --sample_every 100
```

## Outputs

The outputs of each step of the workflow will be stored in a `results` folder within your 
project directory. Each step of the workflow will store its outputs in a 
`sub-folder` of the `results` folder named after that sub-step:

* `data` sub-folder will store the output of the first step where the metadata
 and sequence data are merged, and clusters of size greater than or equal to the
 minimum size are identified. This first step will output as many files as there
 are clusters that meet the minimum size criterion.
 
 * `likelihood` sub-folder will store the output calculating
  occurrence of all the possible pairwise base permutation patterns (e.g., AA, CC, TT, GG, AC, 
  etc.) --- this tally is used in the likelihood function and calculating it once
  avoids having to calculate repeatedly throughout the MCMC. The calculation is
  added to the data from the previous step.
  
* `conditions` sub-folder adds to the data from the previous step different assumptions
  about the `w_dens` and `f_dens` distribution. At the moment, four conditions 
  are tested (all geared towards SARS-CoV-2). This leads to four files per 
  cluster.
  
* `mcmc` sub-folder stores the output for each MCMC run.

* `tidy_output` sub-folder stores the output of the MCMC data in a tidy format. 
  Note that here, all column names produced by `outbreaker2` are prefixed with
  `outb2_`, and all column names originated from the data are prefixed with 
  `meta_`

* `cluster_output` sub-folder groups data from different conditions from the same
  cluster into a single data.frame, and outputs it as a CSV file. Here, some
  column names are changed to match requirements for downstream analyses.
  
  
## Authors

* Kurnia Susvitasari
* Jessica Stockdale
* Ben Sobkowiak
* Anders Goncalves da Silva
* Paul Tupper
* Caroline Colijn





