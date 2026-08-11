# Update a sequence of base parameter annual values to account for a time varying deviation effects

Update a sequence of base parameter annual values to account for a time
varying deviation effects

## Usage

``` r
update_basevals_dev(
  base_vals,
  temp_dev,
  dev_seq,
  current_par,
  temp_ctl,
  base_range,
  base_bounds
)
```

## Arguments

- base_vals:

  A vector of base parameter values that will be updated to include the
  impact of a time varying deviations

- temp_dev:

  The time varying parameter lines for the deviations on the base
  parameter

- dev_seq:

  A vector of the parameter deviations to be applied to the base values

- current_par:

  The index of the current parameter being updated

- temp_ctl:

  A subset of the control file representing the parameter section of
  interest (i.e. MG, SR, Q, or Selectivity)

- base_range:

  the difference between the base parameters max and min bounds

- base_bounds:

  The min and max bounds of the base parameter

## Value

A modified parameter series that incorporates the appropriate
deviations.

## Author

Nathan Vaughan
