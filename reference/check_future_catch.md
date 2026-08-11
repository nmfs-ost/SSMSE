# Check future catch smaller than the last year's population size.

Note that it could still be possible to take out too much catch from the
population, so this may not catch all instances of too much catch

## Usage

``` r
check_future_catch(catch, OM_dir, catch_units = "bio", datfile = NULL)
```

## Arguments

- catch:

  A dataframe of catch values and its associated information to add to
  the OM. The column names are the same as in an SS3 data file (e.g.,
  year, season, fleet, catch, catch_se). length of the number of years
  (only works when catch is for 1 fleet)

- OM_dir:

  The full path to the OM directory.

- catch_units:

  What units is the catch in? "bio" for biomass or "num" for numbers?
  Defaults to "bio".

- datfile:

  The optional name (as a character string) of the datafile, presumed to
  exist in `OM_dir`. Defaults to NULL, and if is NULL, the function will
  get the datfile name from the starter.ss file in `OM_dir`.

## Author

Kathryn Doering
