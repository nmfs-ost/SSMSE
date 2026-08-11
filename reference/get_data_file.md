# Locate the SS3 data file(s) associated with a model run.

The helper recognizes both legacy and current output layouts, including
`data.ss_new`, `data_echo.ss_new`, `data_expval.ss`, and
`data_boot_001.ss`.

## Usage

``` r
get_data_file(dir, file_type = c("all", "input", "expected", "bootstrap"))
```

## Arguments

- dir:

  Directory containing the SS3 model files.

- file_type:

  Which type of output file to look for. Supported values are `all`,
  `input`, `expected`, and `bootstrap`.

## Value

The first matching data file name in the requested priority order, or
`NULL` if no matching file is present.
