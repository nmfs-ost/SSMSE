# Example Performance Metric: calculate the average SSB over a range of years for each iteration

Example performance metric that calculates the average Spawning Stock
Biomass SSB (units as in the simulations) over a range of years for each
iteration of each scenario in the SSMSE simulation run.

## Usage

``` r
get_SSB_avg(summary, min_yr, max_yr)
```

## Arguments

- summary:

  Summary returned from running
  [`SSMSE_summary_all()`](SSMSE_summary_all.md)

- min_yr:

  The first year to include in the average

- max_yr:

  The last year to include in the average

## Value

A tibble containing the average SSB by iteration and scenario.

## Examples

``` r
if (FALSE) { # \dontrun{
avg_ssb <- get_SSB_avg(run_SSMSE_summary, min_yr = 10, max_yr = 105)
avg_ssb
} # }
```
