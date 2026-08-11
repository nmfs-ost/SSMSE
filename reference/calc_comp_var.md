# Calculate uncertainty and biases in historic composition data

Calculate uncertainty and biases in historic composition data

## Usage

``` r
calc_comp_var(
  data_obs,
  data_exp,
  bins,
  fleets = NULL,
  years = NULL,
  seasons = NULL,
  merge_sexes = TRUE,
  sexes = NULL,
  merge_seasons = TRUE,
  merge_fleets = FALSE
)
```

## Arguments

- data_obs:

  A data frame of observed composition data extracted from SS3 .dat file

- data_exp:

  A data frame of the expected composition data as estimated by an SS3
  assessment model

- bins:

  A vector object including the composition bins

- fleets:

  A vector of the fleet numbers to analyze composition uncertainty for
  (Default is all fleets if NULL)

- years:

  A vector of the years to include when calculating composition
  uncertainty (Default is all years if NULL)

- seasons:

  A vector of the seasons to include when calculating composition
  uncertainty (Default is all years if NULL)

- merge_sexes:

  TRUE/FALSE should sexes be merged to calculate variance and biases
  (Defaults to TRUE)

- sexes:

  A vector of the sexes to analyze composition uncertainty for (Default
  is all sexes if NULL)

- merge_seasons:

  TRUE/FALSE should seasons be merged to calculate variance and biases
  (Defaults to TRUE)

- merge_fleets:

  TRUE/FALSE should fleets be merged to calculate variance and biases
  (Defaults to FALSE)

## Value

A list object with uncertainty and bias characteristics to inform data
simulation.

## Author

Nathan R. Vaughan
