# Get the sampling scheme in a data file.

Determine what the default sampling scheme is for a given data file.
Produces a list object with the sampling scheme, which can be modified,
if desired.

## Usage

``` r
get_init_samp_scheme(
  dat,
  dat_types = c("CPUE", "lencomp", "agecomp", "meanbodywt", "MeanSize_at_Age_obs")
)
```

## Arguments

- dat:

  An SS3 data file

- dat_types:

  Types of data to include
