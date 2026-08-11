# Change dataset from OM into format for EM

Change dataset from OM into format for EM

## Usage

``` r
change_dat(OM_datfile, EM_datfile, EM_dir, do_checks = TRUE, verbose = FALSE)
```

## Arguments

- OM_datfile:

  Filename of the datfile produced by the OM within the EM_dir.

- EM_datfile:

  Filename of the datfile from the original EM within the EM_dir.

- EM_dir:

  Absolute or relative path to the Estimation model directory.

- do_checks:

  Should checks on the data be performed? Defaults to TRUE.

- verbose:

  Want verbose output? Defaults to FALSE.

## Value

the new EM data file. Side effect is saving over the OM_dat file in
EM_dir.

## Author

Kathryn Doering

## Examples

``` r
if (FALSE) { # \dontrun{
# TODO: Add example
} # }
```
