# Plot comp data, expected values, and sampled data for 1 scenario

Creates a plot that can be used to see how sampling lines up with data
and expected values for the index of abundance

## Usage

``` r
plot_comp_sampling(dir = getwd(), comp_type = c("agecomp", "lencomp"))
```

## Arguments

- dir:

  Path to the directory containing 1 scenario. Defaults to the current
  working directory.

- comp_type:

  Type of composition data, age or length. Defaults to age.

## Value

A list containing 2 components: 1) the ggplot object and 2) the
dataframe used to make the ggplot object

## Author

Kathryn Doering
