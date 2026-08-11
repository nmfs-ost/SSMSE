# Set the initial global, scenario, and iteration seeds

Set the initial global, scenario, and iteration seeds

## Usage

``` r
set_MSE_seeds(seed = NULL, iter_vec)
```

## Arguments

- seed:

  Input a fixed seed to replicate previous simulation runs. Seed can be
  a single value for a global seed, n_scenarios+1 length vector for
  scenario specific and a global seed, n_iterations+n_scenarios+1 length
  vector for iteration scenario and global seeds. Can also be a list
  object with a single value under `seed[["global"]]`, a vector under
  `seed[["scenario"]]`, and a multiple vectors for iteration specific
  seeds under `seed[["iter"]][[1:n_scenarios]]`.

- iter_vec:

  The number of iterations per scenario. A vector of integers in the
  same order as scen_name_vec.

## Value

A list of length 3 with 1) the global seed value; 2) the scenario seed
values; and 3) the iteration seed values.

## Examples

``` r
seeds <- set_MSE_seeds(seed = seq(10, 80, by = 10), iter_vec = c(2, 3))
```
