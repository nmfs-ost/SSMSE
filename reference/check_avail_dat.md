# check all index years/fleets in EM available in OM. (but not vice versa) a general function that can be used

check all index years/fleets in EM available in OM. (but not vice versa)
a general function that can be used

## Usage

``` r
check_avail_dat(
  EM_dat,
  OM_dat,
  list_item = "CPUE",
  colnames = c("year", "month", "index")
)
```

## Arguments

- EM_dat:

  An SS3 data file read in using r4ss for an EM

- OM_dat:

  An SS3 data file read in using r4ss for an OM

- list_item:

  A component in both EM_dat and OM_dat to check values for. This should
  be a single string value.

- colnames:

  The column names of data to append together.

## Author

Kathryn Doering
