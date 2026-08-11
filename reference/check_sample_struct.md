# Check sample_struct_list

Check that list object sample_struct_list has the expected form,
including the correct names, correct column names (as in r4ss), and that
all values in the dataframes are integer or numeric. This does not check
for if numeric or interger values make sense given the model used.

## Usage

``` r
check_sample_struct(
  sample_struct,
  valid_names = list(catch = c("Yr", "Seas", "FltSvy", "SE"), CPUE = c("Yr", "Seas",
    "FltSvy", "SE"), lencomp = c("Yr", "Seas", "FltSvy", "Sex", "Part", "Nsamp"), agecomp
    = c("Yr", "Seas", "FltSvy", "Sex", "Part", "Ageerr", "Lbin_lo", "Lbin_hi", "Nsamp"),
    meanbodywt = c("Yr", "Seas", "FltSvy", "Part", "Type", "Std_in"), MeanSize_at_Age_obs
    = c("Yr", "Seas", "FltSvy", "Sex", "Part", "Ageerr", "N_"))
)
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

- valid_names:

  The list to compare sample_struct to.

## Author

Kathryn Doering
