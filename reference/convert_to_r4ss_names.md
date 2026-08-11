# Convert user input to r4ss data names

Convert user input to r4ss data names

## Usage

``` r
convert_to_r4ss_names(
  sample_struct,
  convert_key = data.frame(df_name = c(rep("catch", 4), rep("CPUE", 4), rep("lencomp",
    6), rep("agecomp", 9), rep("meanbodywt", 6), rep("MeanSize_at_Age_obs", 7)),
    r4ss_name = c("year", "seas", "fleet", "catch_se", "year", "month", "index",
    "se_log", "year", "month", "fleet", "sex", "part", "Nsamp", "year", "month", "fleet",
    "sex", "part", "ageerr", "Lbin_lo", "Lbin_hi", "Nsamp", "year", "month", "fleet",
    "part", "type", "stderr", "year", "month", "fleet", "sex", "part", "ageerr", "N_"),
    sample_struct_name = c("Yr", 
     "Seas", "FltSvy", "SE", "Yr", "Seas", "FltSvy",
    "SE", "Yr", "Seas", "FltSvy", "Sex", "Part", "Nsamp", "Yr", "Seas", "FltSvy", "Sex",
    "Part", "Ageerr", "Lbin_lo", "Lbin_hi", "Nsamp", "Yr", "Seas", "FltSvy", "Part",
    "Type", "Std_in", "Yr", "Seas", "FltSvy", "Sex", "Part", "Ageerr", "N_"),
    stringsAsFactors = FALSE)
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

- convert_key:

  Data frame defining how r4ss names relate to the sample_struct names.
  For now, a 1:1 relationship is assumed.
