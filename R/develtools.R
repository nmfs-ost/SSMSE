# tools for SSMSE developer troubleshooting

#' Change a model from running with par to running without par
#'
#' The intention of this function is to help troubleshooting issues with the
#' par file. It is intended mostly to help troubleshooting while developing
#' the SSMSE package, but may also be helpful with runtime testing.
#' @param orig_mod_dir The original model directory
#' @param new_mod_dir The new model directory (folder need not exist)
test_no_par <- function(orig_mod_dir, new_mod_dir) {
  message(
    "Trying to run model in ", orig_mod_dir, " in a new directory ",
    new_mod_dir, " and starting from values in the ctl file with no",
    " estimation to diagnose problem with the model."
  )
  if (!dir.exists(new_mod_dir)) {
    dir.create(new_mod_dir)
  }
  r4ss::copy_SS_inputs(orig_mod_dir, new_mod_dir, verbose = FALSE)
  start <- r4ss::SS_readstarter(file.path(new_mod_dir, "starter.ss"),
    verbose = FALSE
  )
  start[["init_values_src"]] <- 0 # read inits from ctl instead of par.
  r4ss::SS_writestarter(start,
    dir = file.path(new_mod_dir), overwrite = TRUE,
    verbose = FALSE
  )
  try(run_ss_model(new_mod_dir, "-maxfn 0 -phase 50 -nohess", verbose = FALSE))
  # read in the 2 par files.
  orig_par <- readLines(get_ss_par_file(orig_mod_dir))

  dat_file <- list.files(
    new_mod_dir,
    pattern = "^data(\\.ss_new|_echo\\.ss_new|_expval\\.ss|_boot_[0-9]{3}\\.ss)$"
  )

  if (length(dat_file) > 0) {
    if (file.exists(file.path(new_mod_dir, dat_file))) {
      new_par <- readLines(get_ss_par_file(new_mod_dir))
      if (length(orig_par) != length(new_par)) {
        new_par_names <- grep("^# [^N]", new_par, value = TRUE)
        orig_par_names <- grep("^# [^N]", orig_par, value = TRUE)
        missing_vars <- setdiff(new_par_names, orig_par_names)
        if (length(orig_par) < length(new_par)) {
          msg <- " is missing the "
        }
        if (length(orig_par) > length(new_par)) {
          msg <- " has added "
        }
        stop(
          "Problem with the ss.par file - different number of lines. ",
          "The original par file in ", orig_mod_dir, msg, " parameters: ",
          paste0(missing_vars, collapse = ", ")
        )
      } else {
        stop(
          "Problem with the ss.par file - same number of lines. ",
          "The original par file in ", orig_mod_dir,
          "has the same number of values as the new ",
          "par in ", new_mod_dir, ", so not sure what the issue is."
        )
      }
      # TODO: develop some more sophisticated ways to look at par file diffs.
    }
  } else {
    # problem is not with the par file, but with some other model misspecification.
    stop(
      "Problem with model - not ss.par related. ",
      "Model originally in ", orig_mod_dir, " could not be run from values ",
      "in the ctl file, suggesting the issue is not with bad par file ",
      " specification. Please look at ",
      file.path(orig_mod_dir, "echoinput.sso"), " to see what went wrong."
    )
  }
  invisible(orig_mod_dir)
}

#' function for package developers to update SS3 input files
#' 
#' @param dir_models directory containing the model input files

update_ss3_version <- function(dir_models = "inst/extdata/models") {
  # get list of models
  models <- dir(dir_models, full.names = TRUE)

  # get current executable
  dir_exe <- r4ss::get_ss3_exe(dir = tempdir())
  
  # run with current SS3 version without estimation
  for (i in 1:length(models)) {
    r4ss::run(
      dir = models[i],
      exe = dir_exe,
      skipfinished = FALSE,
      extras = "-nohess -stopph 0"
    )
  }

  # replace original files with renamed ss_new files
  for (i in 1:length(models)) {
    r4ss::copy_SS_inputs(
      dir.old = models[i],
      dir.new = models[i],
      use_ss_new = TRUE,
      overwrite = TRUE
    )
  }
}
