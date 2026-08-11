# Extend the OM forward using next years' catch

Add the EM defined catch values for the next years.

## Usage

``` r
update_OM(
  OM_dir,
  catch = NULL,
  harvest_rate = NULL,
  catch_basis = NULL,
  F_limit = NULL,
  EM_pars = NULL,
  write_dat = TRUE,
  impl_error = NULL,
  verbose = FALSE,
  seed = NULL,
  n_F_search_loops = 20,
  tolerance_F_search = 0.001
)
```

## Arguments

- OM_dir:

  The full path to the OM directory.

- catch:

  A dataframe of catch values and its associated information to add to
  the OM. The column names are the same as in an SS3 data file (e.g.,
  year, season, fleet, catch, catch_se). Must input either a catch
  and/or a harvest rate data frame. If both are input the catch will
  override harvest rate as the management unit but harvest rate will be
  used as a starting guess for search.

- harvest_rate:

  A dataframe of harvest rate (F) values and associated information to
  add to the OM. The column names are as in an SS3 datafile. If harvest
  rate is input without a corresponding catch the OM will assume effort
  based management an use harvest rate directly with implementation
  error added.

- catch_basis:

  data frame with columns year, seas, fleet, basis that specifies if
  catch should reference retained biomass (1) or dead biomass (2). Any
  year/season/fleet not listed will assume a value of 1 referencing
  retained biomass. Entering -99 for any of year, season, or fleet will
  apply the basis across all values of that variable (i.e., a single row
  with -99, -99, -99, 1 would implement retained biomass for all cases)

- F_limit:

  data frame with columns year, fleet, season, limit that specifies a
  maximum F allowed in the OM or a negative value to specify a multiple
  of the historic maximum F. Any year/season/fleet not listed will
  assume a value of 1.5. Entering -99 for any of year, season, or fleet
  will apply the limit across all values of that variable (i.e., a
  single row with -99, -99, -99, -2 would implement a cap of twice the
  historic maximum F for all cases)

- EM_pars:

  a dataframe of parameter value updates to modify OM

- write_dat:

  Should the datafile be overwritten? Defaults to TRUE.

- impl_error:

  The implementation error

- verbose:

  Want verbose output? Defaults to FALSE.

- seed:

  A random seed so that reproducible results are possible.

- n_F_search_loops:

  Number of times to try to find an F that achieves the catches input in
  the OM. Defaults to 20.

- tolerance_F_search:

  How far apart the input catch and achieved catch can be in tried to
  find an F that achieves the catch input in the OM. Defaults to 0.001.

## Value

A new dat list object (format as created by r4ss::SS_readdat) that has
been extended forward as if read in by r4ss function SS_readdat

## Author

Kathryn Doering & Nathan Vaughan
