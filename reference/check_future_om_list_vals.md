# Check structure of a future OM list against the scen_list and standardize output

Checks that a future OM list is valid when compared with the scen_list
inputs

## Usage

``` r
check_future_om_list_vals(future_om_list, scen_list)
```

## Arguments

- future_om_list:

  An optional list of lists including changes that should be made after
  the end year of the input model. Each first level list element
  outlines 1 change to be made to the operating model. To see an
  example, try running
  [`create_future_om_list`](create_future_om_list.md). Defaults to NULL,
  which implies that the model will be extended forward in time assuming
  the original model structure.

- scen_list:

  The list object of scenarios specifying inputs created by
  [`SSMSE::create_scen_list`](create_scen_list.md).

## Value

The future_om_list with implicit arguments made explicit
