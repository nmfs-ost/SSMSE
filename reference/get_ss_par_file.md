# Locate the SS3 parameter file used by a model run.

This helper prefers the current `ss3.par` filename when present and
falls back to the legacy `ss.par` name for older model directories.

## Usage

``` r
get_ss_par_file(dir)
```

## Arguments

- dir:

  Directory containing the SS3 model files.

## Value

A character string giving the full path to the parameter file.
