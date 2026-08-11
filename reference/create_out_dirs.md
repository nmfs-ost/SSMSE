# create the OM directory

Create an OM directory within the out_dir specified (named by the value
of niter)

## Usage

``` r
create_out_dirs(
  out_dir,
  niter,
  OM_name,
  OM_in_dir,
  EM_name = NULL,
  EM_in_dir = NULL
)
```

## Arguments

- out_dir:

  The directory to which to write output. IF NULL, will default to the
  working directory.

- niter:

  The number iteration

- OM_name:

  Name of the operating model (OM).

- OM_in_dir:

  Relative or absolute path to the operating model, if using a model
  outside of the SSMSE package. Should be a string.

- EM_name:

  Name of the EM model

- EM_in_dir:

  Relative or absolute path to the estimation model,

## Value

A list with 2 named components each of length 1 characters. The
components are: OM_dir, where OM will be run, and OM_in_dir, where the
model files will be copied from.
