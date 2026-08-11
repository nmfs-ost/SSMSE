# Package index

## All functions

- [`EM()`](EM.md) : Use EM as the management strategy option.
- [`Interim()`](Interim.md) : Interim assessment management strategy
- [`SSMSE-package`](SSMSE.md) [`SSMSE`](SSMSE.md) : SSMSE: A package for
  Management Strategy Evaluation (MSE) using Stock Synthesis (SS3)
- [`SSMSE_summary_all()`](SSMSE_summary_all.md) : Get results in a list
  for 1 scenario
- [`SSMSE_summary_iter()`](SSMSE_summary_iter.md) : Get results in a
  list for 1 iteration
- [`SSMSE_summary_scen()`](SSMSE_summary_scen.md) : Get results in a
  list for 1 scenario
- [`Sim_comp()`](Sim_comp.md) : Calculate uncertainty and biases in
  historic composition data
- [`add_OM_devs()`](add_OM_devs.md) : Add in future parameter values
- [`add_dev_changes()`](add_dev_changes.md) : Add the deviation changes
  from the list obj to an existing df
- [`add_new_dat()`](add_new_dat.md) : Add new data to an existing EM
  dataset
- [`add_sample_struct()`](add_sample_struct.md) : Add in years of
  sampling data needed
- [`calc_comp_var()`](calc_comp_var.md) : Calculate uncertainty and
  biases in historic composition data
- [`calc_par_trend()`](calc_par_trend.md) : Calculate the parameter
  trend
- [`change_dat()`](change_dat.md) : Change dataset from OM into format
  for EM
- [`change_yrs_fcast()`](change_yrs_fcast.md) : Change the years in the
  forecast file
- [`check_EM_forecast()`](check_EM_forecast.md) : Check structure of
  forecast is suitable to use in the EM
- [`check_OM_dat()`](check_OM_dat.md) : check that an OM data set has at
  least the same data as an estimation model
- [`check_avail_dat()`](check_avail_dat.md) : check all index
  years/fleets in EM available in OM. (but not vice versa) a general
  function that can be used
- [`check_catch_df()`](check_catch_df.md) : Check the catch dataframe
- [`check_convergence()`](check_convergence.md) : Flag potential
  convergence issues in SS3 model runs
- [`check_dir()`](check_dir.md) : Check that the directory for an OM is
  valid
- [`check_future_catch()`](check_future_catch.md) : Check future catch
  smaller than the last year's population size.
- [`check_future_om_list_str()`](check_future_om_list_str.md) : Check
  the general structure of a future OM list and standardize values
- [`check_future_om_list_vals()`](check_future_om_list_vals.md) : Check
  structure of a future OM list against the scen_list and standardize
  output
- [`check_sample_struct()`](check_sample_struct.md) : Check
  sample_struct_list
- [`check_scen_list()`](check_scen_list.md) : Check structure of the
  object scen_list
- [`clean_init_mod_files()`](clean_init_mod_files.md) : clean the
  initial model files
- [`combine_cols()`](combine_cols.md) : function that creates a combined
  column to the list_item of interest
- [`convert_future_om_list_to_devs_df()`](convert_future_om_list_to_devs_df.md)
  : Create the devs dataframe for a scenario and iteration from user
  input
- [`convert_to_r4ss_names()`](convert_to_r4ss_names.md) : Convert user
  input to r4ss data names
- [`copy_model_files()`](copy_model_files.md) : Copy OM and EM model
  files
- [`create_OM()`](create_OM.md) : Create the OM
- [`create_future_om_list()`](create_future_om_list.md) : Helper
  function to create future om list objects
- [`create_out_dirs()`](create_out_dirs.md) : create the OM directory
- [`create_sample_struct()`](create_sample_struct.md) : Create the
  sample_struct list
- [`create_scen_list()`](create_scen_list.md) : Create scen_list object
  to use in run_SSMSE function.
- [`develop_OMs()`](develop_OMs.md) : Develop different operating models
- [`get_EM_catch_df()`](get_EM_catch_df.md) : Get the EM catch data
  frame
- [`get_EM_dat()`](get_EM_dat.md) : Change the OM data to match the
  format of the original EM data
- [`get_F()`](get_F.md) : Get the Fishing mortality from the timeseries
  Report.sso table
- [`get_SSB_avg()`](get_SSB_avg.md) : Example Performance Metric:
  calculate the average SSB over a range of years for each iteration
- [`get_avg_catch()`](get_avg_catch.md) : Example Performance Metric:
  Calculate average catch over a range of years
- [`get_bin()`](get_bin.md) : Get SS3 binary/executable location in
  package
- [`get_catch_cv()`](get_catch_cv.md) : Example Performance Metric:
  Calculate the coefficient of variation of catch
- [`get_catch_sd()`](get_catch_sd.md) : Example Performance Metric:
  Calculate Standard Deviation of Catch
- [`get_data_file()`](get_data_file.md) : Locate the SS3 data file(s)
  associated with a model run.
- [`get_dead_catch()`](get_dead_catch.md) : Get dead catch from the
  timeseries Report.sso table
- [`get_full_sample_struct()`](get_full_sample_struct.md) : Get the full
  sample structure from user input
- [`get_impl_error_matrix()`](get_impl_error_matrix.md) : Put
  implementation error of 0 into a matrix
- [`get_init_samp_scheme()`](get_init_samp_scheme.md) : Get the sampling
  scheme in a data file.
- [`get_input_value()`](get_input_value.md) : return a value from a data
  frame
- [`get_no_EM_catch_df()`](get_no_EM_catch_df.md) : Get the data frame
  of catch for the next iterations when not using an estimation model.
- [`get_performance_metrics()`](get_performance_metrics.md) : get basic
  data to calculate performance metrics
- [`get_rel_SSB_avg()`](get_rel_SSB_avg.md) : Example Performance
  Metric: Calculate the avg relative SSB (SSB/SSB unfished) over a range
  of years for each iteration
- [`get_retained_catch()`](get_retained_catch.md) : Get retained catch
  from the timeseries Report.sso table
- [`get_ss_par_file()`](get_ss_par_file.md) : Locate the SS3 parameter
  file used by a model run.
- [`get_total_catch()`](get_total_catch.md) : Example Performance
  Metric: Calculate total catch over a range of years
- [`last_yr_catch()`](last_yr_catch.md) : Last year catch used in the
  future for management strategy
- [`locate_in_dirs()`](locate_in_dirs.md) : Locate the OM model files
- [`match_parname()`](match_parname.md) : Match parameter name to
  parameter names in the par file
- [`no_catch()`](no_catch.md) : No Catch in the future management
  strategy
- [`parse_MS()`](parse_MS.md) : Parse management strategy options
- [`plot_comp_sampling()`](plot_comp_sampling.md) : Plot comp data,
  expected values, and sampled data for 1 scenario
- [`plot_index_sampling()`](plot_index_sampling.md) : Plot index data,
  expected values, and sampled data for 1 scenario
- [`r4ss_obj_err()`](r4ss_obj_err.md) : Error if object is not an r4ss
  object
- [`rm_sample_struct_hist()`](rm_sample_struct_hist.md) : Remove the
  historical sampling structure
- [`rm_vals()`](rm_vals.md) : remove vals in 2 list components with the
  same name
- [`run_EM()`](run_EM.md) : Run the estimation model
- [`run_OM()`](run_OM.md) : Initial run of the OM
- [`run_SSMSE()`](run_SSMSE.md) : run an MSE using SS3 OMs
- [`run_SSMSE_iter()`](run_SSMSE_iter.md) : Run one iteration of an MSE
  using SS3 OM
- [`run_SSMSE_scen()`](run_SSMSE_scen.md) : Run an MSE scenario using
  SS3 OM
- [`sample_vals()`](sample_vals.md) : Sample vals from normal random,
  lognormal random, or modified AR-1 process.
- [`set_MSE_seeds()`](set_MSE_seeds.md) : Set the initial global,
  scenario, and iteration seeds
- [`test_no_par()`](test_no_par.md) : Change a model from running with
  par to running without par
- [`update_OM()`](update_OM.md) : Extend the OM forward using next
  years' catch
- [`update_basevals_blocks()`](update_basevals_blocks.md) : Update a
  sequence of base parameter annual values to account for a time varying
  block effects
- [`update_basevals_dev()`](update_basevals_dev.md) : Update a sequence
  of base parameter annual values to account for a time varying
  deviation effects
- [`update_basevals_env()`](update_basevals_env.md) : Update a sequence
  of base parameter annual values to account for a time varying
  environmental effects
- [`update_ss3_version()`](update_ss3_version.md) : function for package
  developers to update SS3 input files
