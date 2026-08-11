# Create the devs dataframe for a scenario and iteration from user input

This function parses user inputs to convert it into a dataframe of
deviations.

## Usage

``` r
convert_future_om_list_to_devs_df(
  future_om_list,
  scen_name,
  niter,
  om_mod_path,
  nyrs,
  global_seed = 123
)
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

- scen_name:

  The scenario name

- niter:

  The iteration number

- om_mod_path:

  Path to the OM files. Used to reference parameter names.

- nyrs:

  The total number of years that the model will be extended forward.

- global_seed:

  A global seed to set, then pull new seeds from. Defaults to 123.

## Value

A list including 3 dataframes and one list: devs_df, the additive
deviations relative to the base values; base_df, the base values of the
parameter with deviations; abs_df, the absolute future values by year
(first col) and parameter (parameters in different cols). Also includes
a modified version of the future_om_list which includes the seed applied
to each list component (note that this is not the ultimate seed used for
sampling, as additional) seeds are generated from this seed based on the
scenario, iteration, and option for randomness (replicate across
scenarios or randomize across scenarios). Note that no OM files are
modified or created as part of this function (i.e., it does not have
side effects).

## Author

Kathryn Doering
