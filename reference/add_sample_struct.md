# Add in years of sampling data needed

Add in years of sampling data needed

## Usage

``` r
add_sample_struct(sample_struct, dat)
```

## Arguments

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

- dat:

  A datafile as read in by r4ss::SS_readdat
