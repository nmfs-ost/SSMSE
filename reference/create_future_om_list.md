# Helper function to create future om list objects

The future_om_list objects specify changes to make in the future to the
OM as the OM is extended forward in time. In particular, this function
helps users create these objects. For now, just returns examples based
on cod model that comes with SSMSE. To learn more about the options
available for future_om_list object, see the [structure of
future_om_list section of the SSMSE user
manual](https://nmfs-ost.github.io/SSMSE/manual/SSMSE.html#structurefutureom).

## Usage

``` r
create_future_om_list(
  example_type = c("model_change", "custom"),
  list_length = 1
)
```

## Arguments

- example_type:

  Type of example future_om_list object to create. Options are
  "model_change" or "custom". Defaults to "model_change".

- list_length:

  The length of the example list to create. Defaults to 1. For now, just
  replicates the same list.

## Examples

``` r
example_future_om_list <-
  create_future_om_list(example_type = "custom", list_length = 2)
```
