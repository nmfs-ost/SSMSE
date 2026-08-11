# Add the deviation changes from the list obj to an existing df

Add the deviation changes from the list obj to an existing df

## Usage

``` r
add_dev_changes(fut_list, scen, iter, parlist, dat, vals_df, nyrs, ctl)
```

## Arguments

- fut_list:

  A single change input

- scen:

  The scenario name

- iter:

  The iteration name

- parlist:

  A parameter file as read in by
  [`r4ss::SS_readpar_3.30`](https://r4ss.github.io/r4ss/reference/SS_readpar_3.30.html)

- dat:

  A data file as read in by
  [`r4ss::SS_readdat`](https://r4ss.github.io/r4ss/reference/SS_readdat.html).

- vals_df:

  The dataframe with future om values

- nyrs:

  The number of years to extend the model forward

- ctl:

  A control file as read in by
  [`r4ss::SS_readctl`](https://r4ss.github.io/r4ss/reference/SS_readctl.html).

## Value

A modified version of vals_df with the new changes applied.

## Author

Kathryn Doering
