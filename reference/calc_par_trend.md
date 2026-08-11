# Calculate the parameter trend

Calculate the parameter trend

## Usage

``` r
calc_par_trend(
  val_info,
  val_line = c("mean", "sd", "cv", "ar_1_phi"),
  ref_parm_value,
  vals_df,
  parname,
  parlist,
  ctl,
  par_section,
  dat
)
```

## Arguments

- val_info:

  The line in the input df containing info about the parameter.

- val_line:

  Which line in val info to use.

- ref_parm_value:

  This is the historic parameter that the end trend value. Can be NA if
  the there is no line in val_info for the given parameter

- vals_df:

  The dataframe of the parameter values by year. Use to get start val
  and last year

- parname:

  Name of the parameter with devs from the SS3 model. will reference, if
  using a relative method.

- parlist:

  A parameter file as read in by
  [`r4ss::SS_readpar_3.30`](https://r4ss.github.io/r4ss/reference/SS_readpar_3.30.html)

- ctl:

  A control file as read in by
  [`r4ss::SS_readctl`](https://r4ss.github.io/r4ss/reference/SS_readctl.html).

- par_section:

  Which parameter section should this variable be in?

- dat:

  A data file as read in by
  [`r4ss::SS_readdat`](https://r4ss.github.io/r4ss/reference/SS_readdat.html).

## Value

A vector of values with length ncol(vals_df), the number of future
years.

## Author

Kathryn Doering
