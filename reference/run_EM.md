# Run the estimation model

Runs the estimation model and performs checks if desired.

## Usage

``` r
run_EM(
  EM_dir,
  hess = FALSE,
  check_converged = TRUE,
  set_use_par = FALSE,
  verbose = FALSE
)
```

## Arguments

- EM_dir:

  Absolute or relative path to the estimation model directory

- hess:

  Get the hessian during model run? Defaults to FALSE. Not estimating
  the hessian will speed up the run, but no estimates of error will be
  generated.

- check_converged:

  Perform checks to see if the model converged? Defaults to TRUE.

- set_use_par:

  Should input values be read from the .par file? If TRUE, will change
  setting in the starter file; otherwise, will use the setting already
  in the starter file, which may or may not read from the .par file.

- verbose:

  Want verbose output? Defaults to FALSE.

## Author

Kathryn Doering
