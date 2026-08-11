# Initial run of the OM

This function is used to initialize the OM and get either expected
values or bootstrap.

## Usage

``` r
run_OM(
  OM_dir,
  boot = TRUE,
  nboot = 1,
  init_run = FALSE,
  verbose = FALSE,
  debug_par_run = TRUE,
  sample_catch = FALSE,
  seed = NULL
)
```

## Arguments

- OM_dir:

  The full path to the OM directory.

- boot:

  Return the bootstrap dataset? If TRUE, function returns the number
  bootstrapped dataset specified in `nboot`. If FALSE, it returns the
  expected values.

- nboot:

  The number bootstrapped data set. This value is only used if
  `boot = TRUE`. Note that this numbering does NOT correspond with the
  numbering in section of r4ss::SS_readdat. E.g., specifying section = 3
  in SS_readdat is equivalent to specifying nboot = 1.

- init_run:

  Is this the initial iteration of the OM? Defaults to FALSE.

- verbose:

  Want verbose output? Defaults to FALSE.

- debug_par_run:

  If set to TRUE, and the run fails, a new folder called error_check
  will be created, and the model will be run from control start values
  instead of ss3.par. The 2 par files are then compared to help debug
  the issue with the model run. Defaults to TRUE.

- sample_catch:

  Should catch be sampled or fixed at the OM values? Defaults to FALSE.

- seed:

  A random seed so that reproducible results are possible.

## Author

Kathryn Doering
