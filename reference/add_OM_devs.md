# Add in future parameter values

Add in future parameter values

## Usage

``` r
add_OM_devs(ctl, dat, parlist, timeseries, future_om_dat)
```

## Arguments

- ctl:

  A control file as read in by
  [`r4ss::SS_readctl`](https://r4ss.github.io/r4ss/reference/SS_readctl.html).

- dat:

  A data file as read in by
  [`r4ss::SS_readdat`](https://r4ss.github.io/r4ss/reference/SS_readdat.html).

- parlist:

  A parameter file as read in by
  [`r4ss::SS_readpar_3.30`](https://r4ss.github.io/r4ss/reference/SS_readpar_3.30.html)

- timeseries:

  The timeseries table from
  [`r4ss::SS_output()`](https://r4ss.github.io/r4ss/reference/SS_output.html).

- future_om_dat:

  A data frame with random sample data for future parameter

## Author

Nathan Vaughan
