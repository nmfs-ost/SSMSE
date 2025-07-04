#' run an MSE using SS3 OMs
#'
#' High level function to run a management strategy evaluation using Stock
#' Synthesis as the Operating model(s). For more examples and information on how
#' to use SSMSE, see the [SSMSE user manual](https://nmfs-ost.github.io/SSMSE/manual/).
#' @template run_SSMSE_scen_list
#' @template seed_input
#' @template future_om_list
#' @template run_SSMSE_params_all
#' @template parallel
#' @template verbose
#' @template F_search_loops
#' @export
#' @author Kathryn Doering & Nathan Vaughan
#' @examples
#' \dontrun{
#' my_dir <- file.path(tempdir(), "ex-run_SSMSE")
#' dir.create(my_dir)
#' # For the EM, use the specified data structure
#' my_sample_struct_list <- list(
#'   NULL,
#'   list(
#'     catch = data.frame(
#'       Yr = 101:106,
#'       Seas = 1,
#'       FltSvy = 1,
#'       SE = 0.05
#'     ),
#'     CPUE = data.frame(
#'       Yr = c(102, 105),
#'       Seas = 7,
#'       FltSvy = 2,
#'       SE = 0.01
#'     ),
#'     lencomp = data.frame(
#'       Yr = c(102, 105), Seas = 1,
#'       FltSvy = 1, Sex = 0,
#'       Part = 0, Nsamp = 100
#'     ),
#'   )
#' )
#' # Use the default parameter values, except for the once specified.
#' # Note that the scen_list, either specified or internally created in the
#' # function is returned.
#' input_list <- run_SSMSE(
#'   scen_name_vec = c("scen_1", "scen_2"),
#'   out_dir_scen_vec = my_dir,
#'   iter_vec = c(2, 2),
#'   OM_name_vec = c("cod", "cod"),
#'   OM_in_dir_vec = NULL,
#'   EM_name_vec = c(NA, "cod"),
#'   EM_in_dir_vec = NULL,
#'   MS_vec = c("no_catch", "EM"),
#'   use_SS_boot_vec = TRUE,
#'   nyrs_vec = 6,
#'   nyrs_assess_vec = 3,
#'   scope = c("2", "1", "3"),
#'   rec_dev_pattern = c(
#'     "none", "rand", "AutoCorr_rand",
#'     "AutoCorr_Spec", "vector"
#'   ),
#'   rec_dev_pars = NULL,
#'   impl_error_pattern = c("none", "rand", "user"),
#'   impl_error_pars = NULL,
#'   verbose = FALSE,
#'   seed = NULL,
#'   sample_struct_list = my_sample_struct_list
#' )
#' unlink(my_dir, recursive = TRUE)
#' }
run_SSMSE <- function(scen_name_vec,
                      out_dir_scen_vec = NULL,
                      iter_vec,
                      OM_name_vec = NULL,
                      OM_in_dir_vec = NULL,
                      EM_name_vec = NULL,
                      EM_in_dir_vec = NULL,
                      run_EM_last_yr = FALSE,
                      MS_vec = c("EM", "no_catch", "Interim"),
                      custom_MS_source = NULL,
                      use_SS_boot_vec = TRUE,
                      nyrs_vec,
                      nyrs_assess_vec,
                      sample_struct_list = NULL,
                      future_om_list = NULL,
                      sample_struct_hist_list = NULL,
                      sample_catch_vec = FALSE,
                      interim_struct_list = NULL,
                      verbose = FALSE,
                      seed = NULL,
                      n_F_search_loops = 20,
                      tolerance_F_search = 0.001,
                      run_parallel = FALSE,
                      n_cores = NULL) {
  if (!is.null(custom_MS_source)) {
    source(custom_MS_source)
  }
  # input checks
  if (!all(MS_vec %in% c("EM", "no_catch", "Interim"))) {
    invalid_MS <- MS_vec[unlist(lapply(MS_vec, function(x) !exists(x)))]
    invalid_MS <- invalid_MS[!invalid_MS %in% c("EM", "no_catch", "Interim")]
    if (length(invalid_MS) > 0) {
      stop(
        "Invalid management strategies in MS_vec: ",
        paste0(invalid_MS, collapse = ", "), ". Please check your function",
        " names and make sure they are all available in the global ",
        "environment, if not built into SSMSE."
      )
    }
  }

  # if(length(sample_catch_vec) == 1) {
  #   sample_catch_vec <-
  #     rep(sample_catch_vec, length.out = length(scen_name_vec))
  # }


  # make sure the output directories exist
  result <- lapply(out_dir_scen_vec, function(x) if (!dir.exists(x)) dir.create(x, showWarnings = FALSE))

  # check and add implicit inputs to the future_om_list
  future_om_list <- check_future_om_list_str(future_om_list = future_om_list)
  # Note that all input checks are done in the check_scen_list function.
  # construct scen_list from other parameters.
  scen_list <- create_scen_list(
    scen_name_vec = scen_name_vec,
    out_dir_scen_vec = out_dir_scen_vec,
    iter_vec = iter_vec,
    OM_name_vec = OM_name_vec,
    OM_in_dir_vec = OM_in_dir_vec,
    EM_name_vec = EM_name_vec,
    EM_in_dir_vec = EM_in_dir_vec,
    MS_vec = MS_vec,
    use_SS_boot_vec = use_SS_boot_vec,
    nyrs_vec = nyrs_vec,
    nyrs_assess_vec = nyrs_assess_vec,
    sample_struct_list = sample_struct_list,
    sample_struct_hist_list = sample_struct_hist_list,
    sample_catch_vec = sample_catch_vec,
    interim_struct_list = interim_struct_list
  )
  # check list and change if need to duplicate values.
  scen_list <- check_scen_list(scen_list, verbose = verbose)
  # check future OM list with the scen_list
  future_om_list <- check_future_om_list_vals(
    future_om_list = future_om_list,
    scen_list = scen_list
  )

  # First reset the R random seed
  set.seed(seed = NULL)
  # Now set the global, scenario, and iteration seeds that will be used as needed
  seed <- set_MSE_seeds(
    seed = seed,
    iter_vec = unlist(lapply(scen_list, function(scen) scen["iter"]))
  )

  # make sure values are the correct length
  nyrs_vec <- unlist(lapply(scen_list, function(scen) scen["nyrs"]))
  nyrs_assess_vec <- unlist(lapply(scen_list, function(scen) scen["nyrs_assess"]))
  iter_vec <- unlist(lapply(scen_list, function(scen) scen["iter"]))

  for (i in seq_along(scen_list)) {
    scen_seed <- vector(mode = "list", length = 3)
    names(scen_seed) <- c("global", "scenario", "iter")
    scen_seed[["global"]] <- seed[["global"]]
    scen_seed[["scenario"]] <- seed[["scenario"]][i]
    scen_seed[["iter"]] <- seed[["iter"]][[i]]
    scen_list[[i]][["scen_seed"]] <- scen_seed
  }

  if (run_parallel) {
    if (!is.null(n_cores)) {
      n_cores <- min(max(n_cores, 1), (parallel::detectCores() - 1))
      cl <- parallel::makeCluster((n_cores), setup_strategy = "sequential")
      doParallel::registerDoParallel(cl, cores = (n_cores))
    } else {
      cl <- parallel::makeCluster((parallel::detectCores() - 1), setup_strategy = "sequential")
      doParallel::registerDoParallel(cl, cores = (parallel::detectCores() - 1))
    }
  }
  # pass each scenario to run
  for (i in seq_along(scen_list)) {
    tmp_scen <- scen_list[[i]]
    # run for each scenario
    return_df <- run_SSMSE_scen(
      scen_name = names(scen_list)[i],
      nscen = i,
      out_dir_scen = tmp_scen[["out_dir_scen"]],
      iter = tmp_scen[["iter"]],
      OM_name = tmp_scen[["OM_name"]],
      OM_in_dir = tmp_scen[["OM_in_dir"]],
      EM_name = tmp_scen[["EM_name"]],
      EM_in_dir = tmp_scen[["EM_in_dir"]],
      MS = tmp_scen[["MS"]],
      custom_MS_source = custom_MS_source,
      use_SS_boot = tmp_scen[["use_SS_boot"]],
      nyrs = tmp_scen[["nyrs"]],
      nyrs_assess = tmp_scen[["nyrs_assess"]],
      # rec_devs_scen = tmp_scen[["rec_devs"]],
      # impl_error = tmp_scen[["impl_error"]],
      scen_seed = tmp_scen[["scen_seed"]],
      sample_struct = tmp_scen[["sample_struct"]],
      future_om_list = future_om_list,
      sample_struct_hist = tmp_scen[["sample_struct_hist"]],
      sample_catch = tmp_scen[["sample_catch"]],
      interim_struct = tmp_scen[["interim_struct"]],
      run_EM_last_yr = run_EM_last_yr,
      n_F_search_loops = n_F_search_loops,
      tolerance_F_search = tolerance_F_search,
      verbose = verbose,
      run_parallel = run_parallel,
      n_cores = n_cores
    )
    scen_list[[i]][["errored_iterations"]] <- return_df
  }
  if (run_parallel) {
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  }
  message("Completed all SSMSE scenarios")
  invisible(scen_list)
}

#' Run an MSE scenario using SS3 OM
#'
#' High level function to run 1 scenario (but potentially many iterations) for
#' a management strategy evaluation using Stock Synthesis as the Operating Model
#' @template run_SSMSE_params_all
#' @template run_SSMSE_scen_iter_params
#' @template parallel
#' @param out_dir_scen The directory to which to write output. IF NULL, will
#' @template OM_EM_in_dir
#' @template OM_name
#' @param iter The number of iterations for the scenario. A single integer
#'  value.
#' @template future_om_list
#' @param scen_seed List containing fixed seeds for this scenario and its iterations.
#' @template MS
#' @template sample_struct
#' @template sample_struct_hist
#' @template verbose
#' @template F_search_loops
#' @template sample_catch
#' @author Kathryn Doering & Nathan Vaughan
#' @examples
#' \dontrun{
#' # Create a temporary folder for the output and set the working directory:
#' temp_path <- file.path(tempdir(), "run_SSMSE_scen-example")
#' dir.create(temp_path, showWarnings = FALSE)
#'
#' # run 2 iteration and 1 scenario of SSMSE
#' run_SSMSE_scen(
#'   scen_name = "no_catch",
#'   iter = 1:2,
#'   OM_name = "cod",
#'   MS = "no_catch",
#'   out_dir_scen = temp_path,
#'   nyrs = 6,
#'   nyrs_assess = 3
#' )
#' unlink(temp_path, recursive = TRUE)
#' }
#'
run_SSMSE_scen <- function(scen_name = "scen_1",
                           nscen = 1,
                           out_dir_scen = NULL,
                           iter = 2,
                           OM_name = "cod",
                           OM_in_dir = NULL,
                           EM_name = NULL,
                           EM_in_dir = NULL,
                           run_EM_last_yr = FALSE,
                           MS = "no_catch",
                           custom_MS_source = NULL,
                           use_SS_boot = TRUE,
                           nyrs = 100,
                           nyrs_assess = 3,
                           scen_seed = NULL,
                           sample_struct = NULL,
                           future_om_list = NULL,
                           sample_struct_hist = NULL,
                           sample_catch = FALSE,
                           interim_struct = NULL,
                           verbose = FALSE,
                           run_parallel = FALSE,
                           n_cores = NULL,
                           n_F_search_loops = 20,
                           tolerance_F_search = 0.001) {
  # input checks
  assertive.types::assert_is_a_string(scen_name)
  assertive.properties::assert_is_atomic(iter)
  assertive.types::assert_is_any_of(iter, classes = c("numeric", "integer"))
  if (!is.null(OM_name)) assertive.types::assert_is_a_string(OM_name)
  assertive.types::assert_is_a_bool(use_SS_boot)
  if (!is.null(EM_name)) assertive.types::assert_is_a_string(EM_name)
  if (!is.null(EM_in_dir)) assertive.types::assert_is_a_string(EM_in_dir)
  assertive.types::assert_is_a_string(MS)
  if (!is.null(out_dir_scen)) assertive.types::assert_is_a_string(out_dir_scen)
  assertive.types::assert_is_any_of(nyrs, classes = c("numeric", "integer"))
  assertive.properties::assert_is_of_length(nyrs, 1)
  assertive.types::assert_is_any_of(nyrs_assess, classes = c("numeric", "integer"))
  assertive.properties::assert_is_of_length(nyrs_assess, 1)
  if (!is.null(sample_struct)) assertive.types::assert_is_list(sample_struct)
  if (!is.null(sample_struct_hist)) assertive.types::assert_is_list(sample_struct_hist)
  assertive.types::assert_is_a_bool(verbose)

  # create the out_dir to store all files for all iter in the scenario.
  if (is.null(out_dir_scen)) {
    out_dir_iter <- scen_name
  } else {
    out_dir_iter <- file.path(out_dir_scen, scen_name)
  }
  dir.create(out_dir_iter, showWarnings = FALSE)
  # find the first iteration, in case other iterations have been run previously
  # for this scenario.
  all_folds <- list.dirs(path = out_dir_iter, full.names = FALSE)
  existing_iters <- grep("^[0-9]+$", all_folds, value = TRUE)
  if (length(existing_iters) > 0) {
    max_prev_iter <- max(as.integer(existing_iters))
    message(
      "Previous iterations found in folder for scenario ", scen_name, ".",
      " First iteration folder will be ", max_prev_iter + 1, "."
    )
  } else {
    max_prev_iter <- 0
  }
  # make a dataframe to store dataframes that error
  return_val <- vector(mode = "list", length = iter)
  if (run_parallel) {
    return_val <- foreach::`%dopar%`(
      foreach::foreach(i = seq_len(iter), .errorhandling = "pass"),
      {
        iter_seed <- vector(mode = "list", length = 3)
        names(iter_seed) <- c("global", "scenario", "iter")
        iter_seed[["global"]] <- scen_seed[["global"]]
        iter_seed[["scenario"]] <- scen_seed[["scenario"]]
        iter_seed[["iter"]] <- scen_seed[["iter"]][i]
        run_SSMSE_iter(
          out_dir = out_dir_iter,
          OM_name = OM_name,
          OM_in_dir = OM_in_dir,
          EM_name = EM_name,
          EM_in_dir = EM_in_dir,
          MS = MS,
          custom_MS_source = custom_MS_source,
          use_SS_boot = use_SS_boot,
          nyrs = nyrs,
          nyrs_assess = nyrs_assess,
          nscen = nscen,
          scen_name = scen_name,
          niter = max_prev_iter + i,
          run_EM_last_yr = run_EM_last_yr,
          iter_seed = iter_seed,
          sample_struct = sample_struct,
          future_om_list = future_om_list,
          sample_struct_hist = sample_struct_hist,
          sample_catch = sample_catch,
          interim_struct = interim_struct,
          n_F_search_loops = n_F_search_loops,
          tolerance_F_search = tolerance_F_search,
          verbose = verbose
        )
      }
    )
  } else {
    for (i in seq_len(iter)) {
      # for(i in seq_len(iter)){
      iter_seed <- vector(mode = "list", length = 3)
      names(iter_seed) <- c("global", "scenario", "iter")
      iter_seed[["global"]] <- scen_seed[["global"]]
      iter_seed[["scenario"]] <- scen_seed[["scenario"]]
      iter_seed[["iter"]] <- scen_seed[["iter"]][i]
      return_val[[i]] <- tryCatch(run_SSMSE_iter(
        out_dir = out_dir_iter,
        OM_name = OM_name,
        OM_in_dir = OM_in_dir,
        EM_name = EM_name,
        EM_in_dir = EM_in_dir,
        MS = MS,
        custom_MS_source = custom_MS_source,
        use_SS_boot = use_SS_boot,
        nyrs = nyrs,
        nyrs_assess = nyrs_assess,
        nscen = nscen,
        scen_name = scen_name,
        niter = max_prev_iter + i,
        run_EM_last_yr = run_EM_last_yr,
        iter_seed = iter_seed,
        sample_struct = sample_struct,
        future_om_list = future_om_list,
        sample_struct_hist = sample_struct_hist,
        sample_catch = sample_catch,
        interim_struct = interim_struct,
        n_F_search_loops = n_F_search_loops,
        tolerance_F_search = tolerance_F_search,
        verbose = verbose
      ), error = function(e) e)
    }
  }
  # summarize the errors from runs.
  run_failed_df <- NULL
  for (i in seq_len(iter)) {
    if ("error" %in% class(return_val[[i]])) {
      message(
        "Iteration ", max_prev_iter + i, " failed in directory ",
        out_dir_iter,
        ". Please delete folders with failed runs before running summary functions."
      )
      tmp_df <- data.frame(
        iteration = max_prev_iter + i,
        scenario = basename(out_dir_iter),
        out_dir = out_dir_iter,
        error = paste(return_val[[i]][["message"]])
      )
      # todo: add more info on why the run failed.
      run_failed_df <- rbind(run_failed_df, tmp_df)
    }
  }
  if (is.null(run_failed_df)) {
    run_failed_df <- "No errored iterations"
  }
  message("Completed all iterations for scenario ", scen_name)
  invisible(run_failed_df)
}

#' Run one iteration of an MSE using SS3 OM
#'
#' High level function to run 1 iteration of a scenario for a management
#' strategy evaluation using Stock Synthesis as the Operating model.
#' @template run_SSMSE_params_all
#' @template run_SSMSE_scen_iter_params
#' @template OM_EM_in_dir
#' @template OM_name
#' @param out_dir The directory to which to write output. IF NULL, will default
#'   to the working directory.
#' @param nyrs_lag Number of years of lag in obtaining data (i.e., the number of years
#'  post EM assessment end yr before advice can be implemented). Defaults to 0.
#' @param niter The iteration number, which is also the name of the folder the
#'  results will be written to.
#' @template future_om_list
#' @template sample_struct
#' @template sample_struct_hist
#' @template sample_catch
#' @param iter_seed List containing fixed seeds for this iteration.
#' @template MS
#' @template verbose
#' @template F_search_loops
#' @author Kathryn Doering & Nathan Vaughan
#' @examples
#' \dontrun{
#' # Create a temporary folder for the output
#' temp_path <- file.path(tempdir(), "run_SSMSE_iter-ex")
#' dir.create(temp_path)
#'
#' # run 1 iteration and 1 scenario of SSMSE
#' run_SSMSE_iter(
#'   OM_name = "cod",
#'   MS = "no_catch",
#'   out_dir = temp_path,
#'   nyrs = 6,
#'   nyrs_assess = 3
#' )
#' unlink(file.path(temp_path, "1"), recursive = TRUE)
#' # run 1 iteration and 1 scenario of SSMSE using an EM.
#' run_SSMSE_iter(
#'   OM_name = "cod",
#'   MS = "EM",
#'   out_dir = temp_path,
#'   EM_name = "cod",
#'   nyrs = 6,
#'   nyrs_assess = 3,
#'   sample_struct = list(
#'     catch = data.frame(Yr = 101:106, Seas = 1, FltSvy = 1, SE = 0.05),
#'     CPUE = data.frame(Yr = c(102, 105), Seas = 7, FltSvy = 2, SE = 0.01),
#'     lencomp = data.frame(
#'       Yr = c(102, 105), Seas = 1, FltSvy = 1,
#'       Sex = 0, Part = 0, Nsamp = 100
#'     ),
#'     agecomp = data.frame(
#'       Yr = c(102, 105), Seas = 1, FltSvy = 2,
#'       Sex = 0, Part = 0, Ageerr = 1,
#'       Lbin_lo = -1, Lbin_hi = -1, Nsamp = 50
#'     )
#'   )
#' )
#' unlink(temp_path, recursive = TRUE)
#' }
run_SSMSE_iter <- function(out_dir = NULL,
                           OM_name = "cod",
                           OM_in_dir = NULL,
                           EM_name = NULL,
                           EM_in_dir = NULL,
                           run_EM_last_yr = FALSE,
                           MS = "last_yr_catch",
                           custom_MS_source = NULL,
                           use_SS_boot = TRUE,
                           nyrs = 100,
                           nyrs_assess = 3,
                           nyrs_lag = 0,
                           nscen = 1,
                           scen_name = NULL,
                           niter = 1,
                           iter_seed = NULL,
                           sample_struct = NULL,
                           future_om_list = NULL,
                           sample_struct_hist = NULL,
                           sample_catch = FALSE,
                           interim_struct = NULL,
                           n_F_search_loops = 20,
                           tolerance_F_search = 0.001,
                           verbose = FALSE) {
  # input checks ----
  # checks for out_dir, OM_name, OM_in_dir, EM_name, EM_in_dir done in create_out_dirs
  assertive.types::assert_is_a_bool(use_SS_boot)
  assertive.types::assert_is_a_string(MS)
  assertive.types::assert_is_any_of(nyrs, c("integer", "numeric"))
  assertive.types::assert_is_any_of(nyrs_assess, c("integer", "numeric"))
  assertive.types::assert_is_any_of(niter, c("integer", "numeric"))
  assertive.types::assert_is_any_of(nscen, c("integer", "numeric"))
  if (!is.null(sample_struct)) {
    assertive.types::assert_is_list(sample_struct)
    sample_struct <- check_sample_struct(sample_struct)
  } else {
    if (MS == "EM") {
      stop(
        "sample_struct cannot be NULL when using an EM. Please specify. The ",
        "helper function create_sample_struct can be used to specify."
      )
    }
  }
  assertive.types::assert_is_a_bool(verbose)

  if (!is.null(custom_MS_source)) {
    source(custom_MS_source)
  }

  message("Starting iteration ", niter, ".")
  set.seed((iter_seed[["iter"]][1] + 123))
  # get and create directories, copy model files ----
  # assign or reassign OM_dir and OM_in_dir in case they weren't specified
  # as inputs
  out_loc <- create_out_dirs(
    out_dir = out_dir, niter = niter, OM_name = OM_name,
    OM_in_dir = OM_in_dir, EM_name = EM_name, EM_in_dir = EM_in_dir
  )
  if (is.null(out_loc)) { # out loc returns false if the iteration was already run.
    return(FALSE)
  }
  OM_out_dir <- out_loc[["OM_out_dir"]]
  OM_in_dir <- out_loc[["OM_in_dir"]]
  EM_in_dir <- out_loc[["EM_in_dir"]]
  EM_out_dir <- out_loc[["EM_out_dir"]]
  if (!is.null(EM_out_dir)) {
    EM_out_dir_basename <- strsplit(EM_out_dir, "_init$")[[1]][1]
  } else {
    EM_out_dir_basename <- NULL
  }
  copy_model_files(
    OM_in_dir = OM_in_dir, OM_out_dir = OM_out_dir,
    EM_in_dir = EM_in_dir, EM_out_dir = EM_out_dir
  )
  # clean model files ----
  # want to do this as soon as possible to "fail fast"
  # for now, this just gets rid of -year value observations, but could do
  # other things.
  clean_init_mod_files(
    OM_out_dir = OM_out_dir, EM_out_dir = EM_out_dir,
    overwrite = TRUE
  )

  # convert sample_struct names ----
  # get the full sampling structure for components that the user didn't specify.
  # if meaning is ambiguous, then this will exit on error.
  if (!is.null(sample_struct)) {
    sample_struct <- get_full_sample_struct(
      sample_struct = sample_struct,
      OM_out_dir = OM_out_dir
    )
    # convert to r4ss names
    sample_struct <- convert_to_r4ss_names(sample_struct)
    sample_struct_hist <- convert_to_r4ss_names(sample_struct_hist)
  }
  # Convert the user input parameter modifications into vectors of annual additive deviations
  future_om_dat <- convert_future_om_list_to_devs_df(future_om_list = future_om_list, scen_name = scen_name, niter = niter, om_mod_path = OM_out_dir, nyrs = nyrs, global_seed = (iter_seed[["iter"]][1] + 1234))

  # MSE first iteration ----
  # turn the stock assessment model into an OM
  init_mod <- create_OM(
    OM_out_dir = OM_out_dir, overwrite = TRUE,
    sample_struct_hist = sample_struct_hist,
    sample_struct = sample_struct,
    verbose = verbose, writedat = TRUE, nyrs = nyrs,
    nyrs_assess = nyrs_assess, nscen = nscen,
    scen_name = scen_name, niter = niter,
    future_om_dat = future_om_dat,
    seed = (iter_seed[["iter"]][1] + 1234)
  )
  impl_error <- init_mod[["impl_error"]]
  # Complete the OM run so it can be use for expect values or bootstrap
  if (use_SS_boot == TRUE) {
    OM_dat <- run_OM(
      OM_dir = OM_out_dir, boot = use_SS_boot, nboot = 1,
      sample_catch = sample_catch,
      verbose = verbose, init_run = TRUE, seed = (iter_seed[["iter"]][1] + 12345)
    )
  }
  if (use_SS_boot == FALSE) {
    stop(
      "Currently, only sampling can be done using the bootstrapping ",
      "capabilities within SS3"
    )
    # TODO: add sampling functions then run a future sampling function that would
    # make it into a dataset.
    # The plan is to implement custom random data sampling similarly to how random
    # parameter changes in the OM are achieved. So that users can implement not only
    # sample uncertainty but also systematic biases possible in real sampling.
  }
  message(
    "Finished running and sampling OM for the historical period for ",
    "iteration ", niter, "."
  )
  # get catch/discard using the chosen management strategy ----
  # This can use an estimation model or EM proxy, or just be a simple management
  # strategy

  # nyrs_lag <- 0
  # first_catch_yr <- (OM_dat[["endyr"]] - nyrs + 1)
  # n_yrs_catch <- nyrs_assess + nyrs_lag
  # EM_avail_dat <- OM_dat
  # EM_avail_dat[[endyr]] <- first_catch_yr - (1 + nyrs_lag)

  # TODO: If we want to add in data lag then we will need to add a remove data
  # function I think as well as the ability to put in fixed catches for the
  # interim years using the Forecatch section.

  new_catch_list <- parse_MS(
    MS = MS, EM_out_dir = EM_out_dir, init_loop = TRUE,
    OM_dat = OM_dat, OM_out_dir = OM_out_dir,
    verbose = verbose, nyrs_assess = nyrs_assess,
    interim_struct = interim_struct,
    dat_yrs = (init_mod[["dat"]][["endyr"]] - nyrs + 1):(init_mod[["dat"]][["endyr"]] - nyrs + nyrs_assess),
    seed = (iter_seed[["iter"]][1] + 123456)
  )

  message(
    "Finished getting catch (years ",
    min(new_catch_list[["catch"]][, "year"]), " to ", max(new_catch_list[["catch"]][, "year"]),
    ") to feed into OM for iteration ", niter, "."
  )

  # Next iterations of MSE procedure ----
  # set up all the years when the assessment will be done.
  # remove first value, because done in the initialization stage.
  styr_MSE <- OM_dat[["endyr"]] - nyrs
  assess_yrs <- seq(styr_MSE, styr_MSE + nyrs, nyrs_assess)
  assess_yrs <- assess_yrs[-1]
  # calculate years after the last assessment. The OM will need to run
  # 1 additional time if this is the case.
  extra_yrs <- nyrs - (assess_yrs[length(assess_yrs)] - styr_MSE)
  if (extra_yrs == 0) {
    test_run_EM_yr <- assess_yrs[length(assess_yrs)]
  } else {
    test_run_EM_yr <- NULL
  }
  # Loop over the assessment years.
  for (yr in assess_yrs) {
    if (verbose) {
      message(
        "Extending, running, and sampling from the OM though year ", yr,
        "."
      )
    }
    update_OM(
      OM_dir = OM_out_dir,
      catch = new_catch_list[["catch"]],
      harvest_rate = new_catch_list[["catch_F"]],
      catch_basis = NULL,
      F_limit = NULL,
      EM_pars = new_catch_list[["EM_pars"]],
      write_dat = TRUE,
      impl_error = impl_error,
      verbose = verbose,
      n_F_search_loops = n_F_search_loops,
      tolerance_F_search = tolerance_F_search,
      seed = (iter_seed[["iter"]][1] + 234567 + yr)
    )
    # rerun OM (without estimation), get samples (or expected values)
    if (use_SS_boot == TRUE) {
      if (!is.null(interim_struct)) {
        if (MS == "Interim" & !is.null(interim_struct[["auto_corr"]])) {
          new_OM_dat <- run_OM(
            OM_dir = OM_out_dir, boot = FALSE, nboot = 1,
            sample_catch = sample_catch,
            verbose = verbose, seed = (iter_seed[["iter"]][1] + 345678 + yr)
          )
        } else {
          new_OM_dat <- run_OM(
            OM_dir = OM_out_dir, boot = use_SS_boot, nboot = 1,
            sample_catch = sample_catch,
            verbose = verbose, seed = (iter_seed[["iter"]][1] + 345678 + yr)
          )
        }
      } else {
        new_OM_dat <- run_OM(
          OM_dir = OM_out_dir, boot = use_SS_boot, nboot = 1,
          sample_catch = sample_catch,
          verbose = verbose, seed = (iter_seed[["iter"]][1] + 345678 + yr)
        )
      }
    } else {
      stop(
        "Currently, only sampling can be done using the bootstrapping ",
        "capabilities within SS"
      )
    }
    message(
      "Finished running and sampling OM through year ", max(new_catch_list[["catch"]][, "year"]),
      "."
    )
    if (run_EM_last_yr == FALSE && isTRUE(yr == test_run_EM_yr)) {
      skip_EM_run <- TRUE
    } else {
      skip_EM_run <- FALSE
    }
    if (skip_EM_run == FALSE) {
      # if using an EM, want to save results to a new folder
      if (!is.null(EM_out_dir) & MS != "Interim") {
        new_EM_out_dir <- paste0(EM_out_dir_basename, "_", yr)
        dir.create(new_EM_out_dir)
        success <- copy_model_files(
          EM_in_dir = EM_out_dir,
          EM_out_dir = new_EM_out_dir
        )
        EM_out_dir <- new_EM_out_dir
      }
      # Only want data for the new years: (yr+nyrs_assess):yr
      # create the new dataset to input into the EM
      # loop EM and get management quantities.
      if (!is.null(EM_out_dir_basename)) {
        tmp_EM_init_dir <- paste0(EM_out_dir_basename, "_init")
      } else {
        tmp_EM_init_dir <- NULL
      }

      new_catch_list <- parse_MS(
        MS = MS,
        EM_out_dir = EM_out_dir,
        EM_init_dir = tmp_EM_init_dir,
        OM_dat = new_OM_dat,
        init_loop = FALSE, verbose = verbose,
        nyrs_assess = nyrs_assess,
        OM_out_dir = OM_out_dir,
        dat_yrs = (yr + 1):(yr + nyrs_assess),
        sample_struct = sample_struct,
        interim_struct = interim_struct,
        seed = (iter_seed[["iter"]][1] + 5678901 + yr)
      )
      message(
        "Finished getting catch (years ", (yr + 1), " to ",
        (yr + nyrs_assess), ") to feed into OM for iteration ", niter, "."
      )
    }
  }
  if (extra_yrs > 0) {
    message("Running the OM 1 final time, because last year extends past the last
    assessment.")
    yr <- assess_yrs[length(assess_yrs)] + extra_yrs
    subset_catch_list <- lapply(new_catch_list,
      function(x, yr) new_catch <- x[x[["year"]] <= yr, ],
      yr = yr
    )
    update_OM(
      OM_dir = OM_out_dir,
      catch = subset_catch_list[["catch"]],
      harvest_rate = subset_catch_list[["catch_F"]],
      catch_basis = NULL,
      F_limit = NULL,
      EM_pars = subset_catch_list[["EM_pars"]],
      impl_error = impl_error,
      verbose = verbose,
      n_F_search_loops = n_F_search_loops,
      tolerance_F_search = tolerance_F_search,
      seed = (iter_seed[["iter"]][1] + 6789012)
    )

    # Don't need bootstrapping, b/c not samplling
    run_OM(
      OM_dir = OM_out_dir, boot = FALSE, verbose = verbose,
      sample_catch = sample_catch,
      seed = (iter_seed[["iter"]][1] + 7890123)
    )
  }
  message("Finished iteration ", niter, ".")
  invisible(TRUE)
}
