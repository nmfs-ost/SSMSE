# Add new data to an existing EM dataset

This should be used for the feedback loops when an EM is used.

## Usage

``` r
add_new_dat(
  OM_dat,
  EM_datfile,
  sample_struct,
  EM_dir,
  nyrs_assess,
  do_checks = TRUE,
  new_datfile_name = NULL,
  verbose = FALSE
)
```

## Arguments

- OM_dat:

  An valid SS3 data file read in using r4ss. In particular, this should
  be sampled data.

- EM_datfile:

  Datafile name run in previous iterations with the EM. Assumed to exist
  in EM_dir.

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

- EM_dir:

  Absolute or relative path to the Estimation model directory.

- nyrs_assess:

  The number of years between assessments. E.g., if an assessment is
  conducted every 3 years, put 3 here. A single integer value.

- do_checks:

  Should checks on the data be performed? Defaults to TRUE.

- new_datfile_name:

  An optional name of a file to write the new datafile to. If NULL, a
  new datafile will not be written.

- verbose:

  Want verbose output? Defaults to FALSE.

## Value

A new SS3 datafile containing the data in EM_datfile with new data from
OM_dat appended

## Author

Kathryn Doering
