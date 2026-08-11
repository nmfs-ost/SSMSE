# Get the Fishing mortality from the timeseries Report.sso table

Get the Fishing mortality from the timeseries Report.sso table

## Usage

``` r
get_F(timeseries, fleetnames, fleetnames_all)
```

## Arguments

- timeseries:

  The timeseries table from
  [`r4ss::SS_output()`](https://r4ss.github.io/r4ss/reference/SS_output.html).

- fleetnames:

  A vector of fleet names, in the order they appear in the SS3 model.

- fleetnames_all:

  A vector of ALL fleet names that are in the model in the order that
  they are specified in the model. This vector helps the function know
  which order the fleets appear in the model.

## Value

a list containing: F_df, a long dataframe with F by year, Era, Seas, and
fleet; F_rate, a data frame with F for the time frame of the model only
by year, Seas, and fleet, ordered as the ss3.par file expects; init_F, a
named vector of initial F values by Season and fleet, ordered (and
named) as SS expects; and F_rate_fcast, a dataframe of forecasted F by
year, Seas, and fleet, ordered as SS would expect in F_rate.
