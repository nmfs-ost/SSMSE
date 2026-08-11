# Create the sample_struct list

Create a sampling structure list using the pattern in a data file and a
year range. NAs are added if no pattern is found (and rm_NAs = FALSE).
The types of structure that are added to this list (given their presence
in the dat file) with their names as called in the list object in
parentheses are: catch (catch), relative indices (CPUE), length
composition (lencomp), age composition (agecomp), mean body weight
(meanbodywt), and mean size at age (MeanSize_at_Age_obs). Details for
creating the sample structure list are available in the [sampling
options section of the SSMSE user
manual](https://nmfs-ost.github.io/SSMSE/manual/SSMSE.html#sampling-options).

## Usage

``` r
create_sample_struct(dat, nyrs, rm_NAs = FALSE)
```

## Arguments

- dat:

  An r4ss list object read in using r4ss::SS_readdat() or the path
  (relative or absolute) to an SS3 data file to read in.

- nyrs:

  Number of years beyond the years included in the dat file to run the
  MSE. A single integer value.

- rm_NAs:

  Should all NAs be removed from dataframes? Defaults to FALSE.

## Value

A sample_struct list object, where each list element is a dataframe
containing sampling values. If there were no data for the type, NA is
returned for the element.

## Author

Kathryn Doering

## Examples

``` r
OM_path <- system.file("extdata", "models", "cod", "ss3.dat", package = "SSMSE")
# note there is a warning for lencomp because it does not have a consistent pattern
sample_struct <- create_sample_struct(OM_path, nyrs = 20)
#> Warning: Pattern not found for lencomp: FltSvy 1, Seas 1. Returning NA for year in this dataframe.
print(sample_struct)
#> $catch
#>     Yr Seas FltSvy    SE
#> 1  101    1      1 0.005
#> 2  102    1      1 0.005
#> 3  103    1      1 0.005
#> 4  104    1      1 0.005
#> 5  105    1      1 0.005
#> 6  106    1      1 0.005
#> 7  107    1      1 0.005
#> 8  108    1      1 0.005
#> 9  109    1      1 0.005
#> 10 110    1      1 0.005
#> 11 111    1      1 0.005
#> 12 112    1      1 0.005
#> 13 113    1      1 0.005
#> 14 114    1      1 0.005
#> 15 115    1      1 0.005
#> 16 116    1      1 0.005
#> 17 117    1      1 0.005
#> 18 118    1      1 0.005
#> 19 119    1      1 0.005
#> 20 120    1      1 0.005
#> 
#> $CPUE
#>    Yr Seas FltSvy  SE
#> 1 105    7      2 0.2
#> 2 110    7      2 0.2
#> 3 115    7      2 0.2
#> 4 120    7      2 0.2
#> 
#> $lencomp
#>   Yr Seas FltSvy Sex Part Nsamp
#> 1 NA    1      1   0    0   125
#> 
#> $agecomp
#>    Yr Seas FltSvy Sex Part Ageerr Lbin_lo Lbin_hi Nsamp
#> 1 105    1      2   0    0      1      -1      -1   500
#> 2 110    1      2   0    0      1      -1      -1   500
#> 3 115    1      2   0    0      1      -1      -1   500
#> 4 120    1      2   0    0      1      -1      -1   500
#> 
#> $meanbodywt
#> [1] NA
#> 
#> $MeanSize_at_Age_obs
#> [1] NA
#> 
```
