# Get results in a list for 1 scenario

Get results in a list for 1 iteration, using ss3sim::get_results_iter

## Usage

``` r
SSMSE_summary_all(
  dir = getwd(),
  scenarios = NULL,
  run_parallel = FALSE,
  n_cores = NULL,
  overwrite = FALSE
)
```

## Arguments

- dir:

  Path to the directory containing the scenarios, either relative or
  absolute. Defaults to the working directory.

- scenarios:

  A character vector of scenarios in dir from which to extract
  summaries. If left as NULL, the summaries will be extracted from all
  folders in dir.

- run_parallel:

  Option to use parallel processing on iterations. Defaults to FALSE

- n_cores:

  how many cores to use if running in parallel defaults to n_cores
  available - 1 (also capped at one less than the number of cores
  available)

- overwrite:

  Allow existing files to be overwritten?

## Value

A list of 3 data frames called scalar, ts, and dq (for derived
quantities). These lists contain information for multiple model runs
(estimation models and operating models) for 1 iteration.Also writes 3
.csv files with the contents of this list of dataframes to dir and 3.csv
files with scenario specific results in each of the scenario foldurs..

## See also

[`get_results_all`](https://rdrr.io/pkg/ss3sim/man/get_results_all.html)
