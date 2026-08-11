# Sample vals from normal random, lognormal random, or modified AR-1 process.

Sample vals from normal random, lognormal random, or modified AR-1
process.

## Usage

``` r
sample_vals(mean, sd, ar_1_phi = 0, ndevs, dist = c("normal", "lognormal"))
```

## Arguments

- mean:

  A single value or vector of mean parameters

- sd:

  A single value or vector of sd parameter

- ar_1_phi:

  The phi (coefficient) value for an ar 1 process. Should be between -1
  and 1. 0 means an AR 1 process will NOT be used. 1 indicates a random
  walk. A single value or vector.

- ndevs:

  The number of sampled values to expect

- dist:

  The distribution to sample from.

## Author

Kathryn Doering
