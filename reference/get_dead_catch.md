# Get dead catch from the timeseries Report.sso table

Get dead catch from the timeseries Report.sso table

## Usage

``` r
get_dead_catch(timeseries, units_of_catch)
```

## Arguments

- timeseries:

  The timeseries table from
  [`r4ss::SS_output()`](https://r4ss.github.io/r4ss/reference/SS_output.html).

- units_of_catch:

  From datalist, the catch units. A named list where the names are the
  fleets (to provide an extra check)

## Value

a data frame with retained catch by year, Era, Seas, fleet, and units
(long format)
