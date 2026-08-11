# Calculate uncertainty and biases in historic composition data

Calculate uncertainty and biases in historic composition data

## Usage

``` r
Sim_comp(
  Comp_uncert,
  data_exp,
  bins,
  years = NULL,
  seasons = NULL,
  fleets = NULL,
  sexes = NULL
)
```

## Arguments

- Comp_uncert:

  A list object representing the output from the calc_comp_var function

- data_exp:

  A vector representing the expected composition values for which to
  draw a random observation dataset

- bins:

  A vector object including the composition bins

- years:

  A vector of the years to simulate data for. default is all years in
  data_exp if NULL.

- seasons:

  A vector of the seasons to simulate data for. default is all seasons
  in data_exp if NULL.

- fleets:

  A vector of the fleets to simulate data for. default is all fleets in
  data_exp if NULL.

- sexes:

  A vector of the sexes to simulate data for. default is all sexes in
  data_exp if NULL.

## Value

A list object with uncertainty and bias characteristics to inform data
simulation.

## Author

Nathan R. Vaughan
