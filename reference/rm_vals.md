# remove vals in 2 list components with the same name

From 2 list components with the same name, remove vals that aren't in
the compare object

## Usage

``` r
rm_vals(return_obj, compare_obj, name_in_obj, colnames)
```

## Arguments

- return_obj:

  the object (containing list component of name in obj) that will be
  modified. Only combinations of the columns found in compare object
  will be retained

- compare_obj:

  the object (containing list component of name_in_obj) that return_obj
  will be compared to

- name_in_obj:

  the name of the list elements to use; the same name must be in
  return_obj and compare_obj. This list element must be a data frame
  with the same column names

- colnames:

  The column names within the name_in_obj list components to compare.

## Value

`return_obj[[name_in_obj]]`, modified to only include elements present
in `compare_obj[[name_in_obj]]`.

## Author

Kathryn Doering
