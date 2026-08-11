# Create the OM

This function manipulates the OM as needed so that it can be used as an
operating model.

## Usage

``` r
create_OM(
  OM_out_dir,
  overwrite = TRUE,
  writedat = TRUE,
  verbose = FALSE,
  nyrs = NULL,
  nyrs_assess = NULL,
  nscen = 1,
  scen_name = NULL,
  niter = 1,
  future_om_dat = NULL,
  verify_OM = TRUE,
  sample_struct_hist = NULL,
  sample_struct = NULL,
  seed = NULL
)
```

## Arguments

- OM_out_dir:

  The full path to the directory in which the OM is run.

- overwrite:

  Allow existing files to be overwritten?

- writedat:

  Should a new datafile be written?

- verbose:

  Want verbose output? Defaults to FALSE.

- nyrs:

  Number of years beyond the years included in the OM to run the MSE. A
  single integer value.

- nyrs_assess:

  The number of years between assessments. This is used to structure the
  forecast file for use in the OM.

- nscen:

  The scenario number

- scen_name:

  The scenario name

- niter:

  the iteration number

- future_om_dat:

  An optional data_frame including changes that should be made after the
  end year of the input model. Including parameter variations,
  recruitment deviations, and implementation errors.

- verify_OM:

  Should the model be run without estimation and some basic checks done
  to verify that the OM can run? Defaults to TRUE.

- sample_struct_hist:

  An optional list including which years should be sampled for the
  historical period for the data generated from the OM. If this is left
  as NULL, then the same sampling scheme will be used as in the OM's
  data file. If it is not NULL, then each year.

- sample_struct:

  A optional list including which years, seasons, and fleets should be
  added from the OM into the EM for different types of data. If NULL,
  the data structure will try to be inferred from the pattern found for
  each of the datatypes within the EM datafiles. Include this structure
  for the number of years to extend the model out. Note that the data
  should be specified using the list component names and column names as
  in would be used in
  [`r4ss::SS_readdat()`](https://r4ss.github.io/r4ss/reference/SS_readdat.html).
  The `run_SSMSE_iter` function examples give an example of what this
  structure should be. Running the function
  [`create_sample_struct()`](create_sample_struct.md) will also produce
  a sample_struct object in the correct form. Can be NULL only when MS
  is not EM.

- seed:

  A random seed so that reproducible results are possible.

## Value

A modified datafile

## Author

Kathryn Doering & Nathan Vaughan
