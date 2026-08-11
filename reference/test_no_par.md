# Change a model from running with par to running without par

The intention of this function is to help troubleshooting issues with
the par file. It is intended mostly to help troubleshooting while
developing the SSMSE package, but may also be helpful with runtime
testing.

## Usage

``` r
test_no_par(orig_mod_dir, new_mod_dir)
```

## Arguments

- orig_mod_dir:

  The original model directory

- new_mod_dir:

  The new model directory (folder need not exist)
