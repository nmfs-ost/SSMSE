# Get the EM catch data frame

Get the data frame of catch for the next iterations when using a Stock
Synthesis Estimation model from the Report.sso file.

## Usage

``` r
get_EM_catch_df(EM_dir, dat)
```

## Arguments

- EM_dir:

  Path to the EM files

- dat:

  A SS3 datfile read into R using
  [`r4ss::SS_readdat()`](https://r4ss.github.io/r4ss/reference/SS_readdat.html)

## Value

A data frame of future catch

## Author

Kathryn Doering
