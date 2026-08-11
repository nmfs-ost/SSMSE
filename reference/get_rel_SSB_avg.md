# Example Performance Metric: Calculate the avg relative SSB (SSB/SSB unfished) over a range of years for each iteration

Example performance metric that calculates the average Spawning Stock
Biomass SSB (units as in the simulations) relative to the unfished SSB
over a range of years for each iteration of each scenario in the SSMSE
simulation run.

## Usage

``` r
get_rel_SSB_avg(summary, min_yr, max_yr)
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

A tibble containing the relative avg SSB per year by iteration and
scenario.

## Examples

``` r
if (FALSE) { # \dontrun{
rel_avg_ssb <- get_rel_SSB_avg(run_SSMSE_summary, min_yr = 10, max_yr = 105)
rel_avg_ssb
} # }
```
