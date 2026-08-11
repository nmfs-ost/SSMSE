# Get results in a list for 1 iteration

Get results in a list for 1 iteration, using ss3sim::get_results_iter

## Usage

``` r
SSMSE_summary_iter(dir)
```

## Arguments

- dir:

  Path to the directory for 1 iteration of 1 scenario.

## Value

A list of 3 data frames called scalar, timeseries, and derived (for
derived quantities). These lists contain information for multiple model
runs (estimation models and operating models) for 1 iteration.

## See also

[`get_results_iter`](https://rdrr.io/pkg/ss3sim/man/get_results_iter.html)
