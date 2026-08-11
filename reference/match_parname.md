# Match parameter name to parameter names in the par file

Match parameter name to parameter names in the par file

## Usage

``` r
match_parname(list_pars, parlist)
```

## Arguments

- list_pars:

  the parameter names to find

- parlist:

  A parameter file as read in by
  [`r4ss::SS_readpar_3.30`](https://r4ss.github.io/r4ss/reference/SS_readpar_3.30.html)

## Value

A dataframe containing the parameter name and which object it is in in
the par object.

## Author

Kathryn Doering
