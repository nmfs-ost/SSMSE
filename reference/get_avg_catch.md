# Example Performance Metric: Calculate average catch over a range of years

Example performance metric that calculates average catch over a range of
years in a Stock Synthesis data file. This function aggregates across
fleets, so may not be appropriate for models with multiple fleets.

## Usage

``` r
get_avg_catch(datfile, yrs)
```

## Arguments

- datfile:

  Path to the Stock Synthesis data file containing catch

- yrs:

  A vector containing a range of years. Years are as defined in the
  Stock Synthesis data file.

## Value

The average catch, a number.

## Examples

``` r
if (FALSE) { # \dontrun{
avg_catch <- function(datfile = "ss3model/dat.ss", yrs = 30:75) {
  avg_catch
}
} # }
```
