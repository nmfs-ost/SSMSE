# Check structure of forecast is suitable to use in the EM

Check structure of forecast is suitable to use in the EM

## Usage

``` r
check_EM_forecast(fore, n_flts_catch = NULL)
```

## Arguments

- fore:

  A forecast list read in using r4ss::SS_readforecast

- n_flts_catch:

  The number of fleets with catch. If NULL, this function will skip a
  check requiring this input.

## Value

Function mainly used for side effects, but returns TRUE invisibly if no
errors created.

## Author

Kathryn Doering
