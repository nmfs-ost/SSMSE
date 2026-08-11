# Change the OM data to match the format of the original EM data

This does the technical part of changing the EM data. Note this may be
unnecessary

## Usage

``` r
get_EM_dat(OM_dat, EM_dat, do_checks = TRUE)
```

## Arguments

- OM_dat:

  An SS3 data file read in by as a list read in using r4ss from the
  operating model

- EM_dat:

  An SS3 data file read in by as a list read in using r4ss from the
  estimation model

- do_checks:

  Should checks on the data be performed? Defaults to TRUE.

## Value

A data list in the same format that can be read/written by r4ss that has
index, lcomps, and age comps from OM_dat, but with the same structure as
EM_dat.

## Author

Kathryn Doering
