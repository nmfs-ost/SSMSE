# Develop different operating models

This is a utility to help a user create new operating models starting
from the same model. For now, it is only possible to adjust 1 parameter
value

## Usage

``` r
develop_OMs(
  OM_name = NULL,
  OM_in_dir = NULL,
  out_dir = getwd(),
  par_name,
  par_vals,
  refit_OMs = TRUE,
  hess = FALSE
)
```

## Arguments

- OM_name:

  Name of the operating model (OM).

- OM_in_dir:

  Relative or absolute path to the operating model, if using a model
  outside of the SSMSE package. Should be a string.

- out_dir:

  Path where the new models will be written. Defaults to the current
  working directory.

- par_name:

  Name of the parameter to modify

- par_vals:

  Vector of parameter values to modify in the OM. Assume these will be
  fixed so phase will be set as negative.

- refit_OMs:

  Should the models be refit to data? Defaults to TRUE

- hess:

  Should the hessian be estimated if reffiting the OMs? defaults to
  FALSE
