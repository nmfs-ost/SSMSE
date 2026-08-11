# Change the years in the forecast file

This is both to increment years forward and/or to change absolute years
to relative years.

## Usage

``` r
change_yrs_fcast(
  fore,
  make_yrs_rel = TRUE,
  nyrs_increment = NULL,
  nyrs_fore = NULL,
  mod_styr,
  mod_endyr
)
```

## Arguments

- fore:

  A forecasting file read into R using r4ss::SS_readforecast()

- make_yrs_rel:

  Should the absolute years in the forecast file be changed to relative
  years? Defaults to TRUE.

- nyrs_increment:

  The number of years to increment forecasting period years. If NULL
  (the default value), will not be incremented.

- nyrs_fore:

  The number of years of forecasting to do. If NULL, do not change the
  number of forecasting years already specified in `fore`

- mod_styr:

  The first year of the model

- mod_endyr:

  The last year of the model `fore` assumes when read in. Note that the
  assumed model year will be different for the output if nyrs_increment
  is not NULL.

## Value

A forecasting file as an R list object

## Author

Kathryn Doering
