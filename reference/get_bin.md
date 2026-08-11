# Get SS3 binary/executable location in package

Get the binary/executable location in the package SSMSE. This function
is from [ss3sim](https://github.com/ss3sim/ss3sim).

## Usage

``` r
get_bin(bin_name = "ss3")
```

## Arguments

- bin_name:

  Name of SS3 binary, defaults to "ss3"

## Value

The path to an SS3 binary. If using the GitHub version of the package,
this will be an internal binary. Otherwise, this function will search
for a version of the binary in your path. See the ss3sim vignette.

## Examples

``` r
if (FALSE) { # \dontrun{
get_bin()
} # }
```
