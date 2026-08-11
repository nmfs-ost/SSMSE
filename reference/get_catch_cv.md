# Example Performance Metric: Calculate the coefficient of variation of catch

Example performance metric that calculates the coefficient of variation
(CV) of catch over a range of years in a Stock Synthesis data file. This
function aggregates across fleets, so may not be appropriate for models
with multiple fleets.

## Usage

``` r
get_catch_cv(datfile, yrs)
```

## Arguments

- datfile:

  Path to the Stock Synthesis data file containing catch

- yrs:

  A vector containing a range of years. Years are as defined in the
  Stock Synthesis data file.

## Value

The catch coefficient of variation, a number.

## Examples

``` r
if (FALSE) { # \dontrun{
catch_cv <- get_catch_cv(datfile = "mod/dat.ss", yrs = c(20:50, 75:100))
catch_cv
} # }
```
