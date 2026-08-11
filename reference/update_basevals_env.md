# Update a sequence of base parameter annual values to account for a time varying environmental effects

Update a sequence of base parameter annual values to account for a time
varying environmental effects

## Usage

``` r
update_basevals_env(
  base_vals,
  base_years,
  temp_env,
  current_par,
  timeseries,
  temp_ctl,
  dat,
  base_range,
  base_bounds,
  parlist
)
```

## Arguments

- base_vals:

  A vector of base parameter values that will be updated to include the
  impact of a time varying environmental effects

- base_years:

  A vector of years for which the base values are needed

- temp_env:

  The time varying parameter lines for the environmental effects on the
  base parameter

- current_par:

  The index of the current parameter being updated

- timeseries:

  The timeseries table from
  [`r4ss::SS_output()`](https://r4ss.github.io/r4ss/reference/SS_output.html).

- temp_ctl:

  A subset of the control file representing the parameter section of
  interest (i.e. MG, SR, Q, or Selectivity)

- dat:

  A datafile as read in by r4ss::SS_readdat

- base_range:

  the difference between the base parameters max and min bounds

- base_bounds:

  The min and max bounds of the base parameter

- parlist:

  A parameter file as read in by
  [`r4ss::SS_readpar_3.30`](https://r4ss.github.io/r4ss/reference/SS_readpar_3.30.html)

## Value

A modified parameter series that incorporates the appropriate time
varying environmental effects.

## Author

Nathan Vaughan
