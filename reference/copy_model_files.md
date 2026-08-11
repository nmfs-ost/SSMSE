# Copy OM and EM model files

Copy OM and EM model files from input to output location.

## Usage

``` r
copy_model_files(
  OM_in_dir = NULL,
  OM_out_dir = NULL,
  EM_in_dir = NULL,
  EM_out_dir = NULL,
  verbose = FALSE
)
```

## Arguments

- OM_in_dir:

  Relative or absolute path to the operating model, if using a model
  outside of the SSMSE package. Should be a string.

- OM_out_dir:

  The full path to the directory in which the OM is run.

- EM_in_dir:

  Relative or absolute path to the estimation model,

- EM_out_dir:

  Relative or absolute path to the estimation model, if using a model
  outside of the SSMSE package.

- verbose:

  Want verbose output? Defaults to FALSE.

## Value

TRUE, if copying is successful
