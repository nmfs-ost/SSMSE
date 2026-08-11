# Flag potential convergence issues in SS3 model runs

Does basic checks for convergance of estimation model runs from
[`run_SSMSE()`](run_SSMSE.md) simulations. This function 1) warns if
there are parameters on bounds; 2) warns if the SSB in the EM is 2x as
large or half as small as the OM. Note these warnings may not mean that
the models have not converged, but can flag potential issues that can be
investigated further

## Usage

``` r
check_convergence(summary, min_yr, max_yr)
```

## Arguments

- summary:

  Summary returned from running
  [`SSMSE_summary_all()`](SSMSE_summary_all.md)

- min_yr:

  The first year of SSB checked

- max_yr:

  The last year of SSB checked

## Value

A tibble containing the SSB values in the EM relative to the OM by model
run of each iteration of each scenario.

## Examples

``` r
if (FALSE) { # \dontrun{
check_convergance(SSMSE_summary, min_yr = 101, max_yr = 120)
} # }
```
