# Get results in a list for 1 scenario

Get results in a list for 1 iteration, using ss3sim::get_results_iter

## Usage

``` r
SSMSE_summary_scen(dir = getwd())
```

## Arguments

- dir:

  Path to the directory for 1 scenario, either relative or absolute.
  Defaults to the working directory.

## Value

A list of 3 data frames called scalar, ts, and dq (for derived
quantities). These lists contain information for multiple model runs
(estimation models and operating models) for 1 iteration.Also writes 3
.csv files with the contents of this list of dataframes to dir.

## See also

[`get_results_scenario`](https://rdrr.io/pkg/ss3sim/man/get_results_scenario.html)
