# Remove the historical sampling structure

Remove the historical sampling structure

## Usage

``` r
rm_sample_struct_hist(sample_struct_hist, dat)
```

## Arguments

- sample_struct_hist:

  An optional list including which years should be sampled for the
  historical period for the data generated from the OM. If this is left
  as NULL, then the same sampling scheme will be used as in the OM's
  data file. If it is not NULL, then each year.

- dat:

  The data file, as read in using r4ss
