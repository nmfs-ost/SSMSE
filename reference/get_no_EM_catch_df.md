# Get the data frame of catch for the next iterations when not using an estimation model.

Get the data frame of catch for the next iterations when not using an
estimation model.

## Usage

``` r
get_no_EM_catch_df(OM_dir, yrs, MS = "last_yr_catch")
```

## Arguments

- OM_dir:

  The full path to the OM directory.

- yrs:

  A vector of years

- MS:

  The management strategy to use. Current options are: `"last_yr_catch"`
  which uses the previous year's catch; `"no_catch"` which uses 0 catch;
  `"EM"` which uses an stock synthesis model as the estimation method
  and the management strategy as defined in the forecast file of the
  stock synthesis estimation method; `"Interim"` to modify catch based
  on survey predictions between assessments. Users can also specify
  their own management strategies as a function. For example, if the
  function is called "my_ms" then the user should specify MS = "my_ms"
  and specify the path to the file containing the function in
  `custom_MS_source`.

## Value

A dataframe of future catch.

## Author

Kathryn Doering
