# clean the initial model files

clean the initial model files

## Usage

``` r
clean_init_mod_files(
  OM_out_dir,
  EM_out_dir = NULL,
  MS = "EM",
  overwrite = FALSE
)
```

## Arguments

- OM_out_dir:

  The full path to the directory in which the OM is run.

- EM_out_dir:

  Relative or absolute path to the estimation model, if using a model
  outside of the SSMSE package.

- MS:

  The management strategy to use. Current options are: `"last_yr_catch"`
  which uses the previous year's catch; `"no_catch"` which uses 0 catch;
  `"EM"` which uses an stock synthesis model as the estimation method
  and the management strategy as defined in the forecast file of the
  stock synthesis estimation method; `"Interim"` to modify catch based
  on survey predictions between assessments. Users can also specify
  their own management strategies as a function. For example, if the
  function is called "my_ms" then the user should specify MS = "my_ms"
  and specify the path to the file containing the function in
  `custom_MS_source`.

- overwrite:

  Allow existing files to be overwritten?
