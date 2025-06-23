## EnvirEM is based off of BiasEM and is aimed at improving
## the inclusion of environment driven mortality in the SSMSE 
## process.  

add_new_dat_envir <- function (OM_dat,
                              EM_datfile,
                              sample_struct,
                              EM_dir,
                              nyrs_assess,
                              do_checks = TRUE,
                              new_datfile_name = NULL,
                              verbose = FALSE)
  
{
  if (do_checks) {
    if (OM_dat[["type"]] != "Stock_Synthesis_data_file") {
      r4ss_obj_err("OM_dat", "data list")
    }
  }
  EM_dat <- SS_readdat(file.path(EM_dir, EM_datfile), verbose = FALSE)
  new_EM_dat <- EM_dat
  new_EM_dat[["endyr"]] <- new_EM_dat[["endyr"]] + nyrs_assess
  extracted_dat <- mapply(function(df, df_name, OM_dat) {
    OM_df <- OM_dat[[df_name]]
    if (is.integer(OM_df[1, 3]) | is.numeric(OM_df[1, 3])) {
      OM_df[, 3] <- base::abs(OM_df[, 3])
    } else if (is.character(OM_df[1, 3])) {
      OM_df[, 3] <- as.character(base::abs(as.integer(OM_df[, 
                                                            3])))
    } # end else if
    by_val <- switch(df_name, catch = c("year", "seas", "fleet"),
                     CPUE = c("year", "seas", "index"),
                     discard_data = c("Yr", "Seas","Flt"), 
                     lencomp = c("Yr", "Seas", "FltSvy", "Gender", "Part"), 
                     agecomp = c("Yr","Seas", "FltSvy", "Gender", "Part", 
                                 "Ageerr", "Lbin_lo", "Lbin_hi"), 
                     meanbodywt = c("Year", "Seas", "Fleet", "Partition", "Type"),
                     MeanSize_at_Age_obs = c("Yr", "Seas", "FltSvy", "Gender", "Part", "AgeErr"))
    new_dat <- merge(df, OM_df, by = by_val, all.x = TRUE, all.y = FALSE)
    if ("catch_se.y" %in% colnames(new_dat)) {
      new_dat[["catch_se.x"]] <- NULL
      colnames(new_dat)[which(colnames(new_dat) == "catch_se.y")] <- "catch_se"
    }
    if ("se_log.y" %in% colnames(new_dat)) {
      new_dat[["se_log.x"]] <- NULL
      colnames(new_dat)[which(colnames(new_dat) == "se_log.y")] <- "se_log"
    }
    if ("Nsamp.y" %in% colnames(new_dat)) {
      new_dat[["Nsamp.x"]] <- NULL
      colnames(new_dat)[which(colnames(new_dat) == "Nsamp.y")] <- "Nsamp"
    }
    if ("Std_in.y" %in% colnames(new_dat)) {
      new_dat[["Std_in.x"]] <- NULL
      colnames(new_dat)[which(colnames(new_dat) == "Std_in.y")] <- "Std_in"
    }
    if ("Ignore.y" %in% colnames(new_dat)) {
      new_dat[["Ignore.y"]] <- NULL
      colnames(new_dat)[which(colnames(new_dat) == "Ignore.x")] <- "Ignore"
    }
    if ("N_" %in% colnames(new_dat)) {
      n_col <- which(colnames(new_dat) == "N_")
      new_dat <- new_dat[, -n_col]
    }
    if (any(is.na(new_dat))) {
      warning("Some values specified in sample_struct (list component ", 
              df_name, ") were not found in OM_dat, so they will not be added to ", 
              "the EM_dat.")
      new_dat <- na.omit(new_dat)
    }
    new_dat
  }, df = sample_struct, df_name = names(sample_struct), MoreArgs = list(OM_dat = OM_dat), 
  SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  #extracted_dat takes observations from the OM; then below we use the EM2OM multiplier to put catch back into EM units and remove the EM2OMcatch_basis element so that it doesn't get added to the datafile. 
  
  if (!is.null(extracted_dat$catch)) {
    extracted_dat$catch <- extracted_dat$catch[order(
      base::abs(extracted_dat$catch$fleet),
      base::abs(extracted_dat$catch$year),
      base::abs(extracted_dat$catch$seas)
    ), ]
    tmp_catch <- merge(extracted_dat$catch, sample_struct$EM2OMcatch_bias)
    tmp_catch <- tmp_catch[order(
      base::abs(tmp_catch$fleet),
      base::abs(tmp_catch$year),
      base::abs(tmp_catch$seas)
    ), ]
    extracted_dat$catch$catch <- tmp_catch$catch / tmp_catch$bias ### EM2OM edits to get EM catch back to OM
  }
  
  if (!is.null(extracted_dat$discard_data)) {
    extracted_dat$discard_data <- extracted_dat$discard_data[order(
      base::abs(extracted_dat$discard_data$Flt),
      base::abs(extracted_dat$discard_data$Yr),
      base::abs(extracted_dat$discard_data$Seas)
    ), ]
    tmp_discard <- merge(extracted_dat$discard_data,
                         sample_struct$EM2OMdiscard_bias)
    tmp_discard <- tmp_discard[order(
      base::abs(tmp_discard$Flt),
      base::abs(tmp_discard$Yr),
      base::abs(tmp_discard$Seas)
    ), ]
    extracted_dat$discard_data$Discard <- tmp_discard$Discard / tmp_discard$bias
  }
  
  if(!is.null(sample_struct$FixedCatchEM)){
    extracted_dat[["catch"]] <- rbind(extracted_dat[["catch"]], sample_struct[["FixedCatchEM"]])
  }
  
  extracted_dat[["EM2OMcatch_bias"]] <- NULL
  extracted_dat[["FixedCatch"]] <- NULL

  for (n in names(extracted_dat)) {
    new_EM_dat[[n]] <- rbind(new_EM_dat[[n]], extracted_dat[[n]])
  }
  if  (!is.null(new_datfile_name)) {
    SS_writedat(new_EM_dat, file.path(EM_dir, new_datfile_name), 
                overwrite = TRUE, verbose = FALSE)
  }
  new_EM_dat
}
# assignInNamespace("add_new_dat_BIAS", add_new_dat_BIAS, "SSMSE") ## ERROR


sample_struct_hist_func<-function(name, dat) {
  df <- dat[[name]]
  if (is.null(df)) {
    return(NA)
  }
  # get year, seas, fleet combo, ignoring -999 values.
  yr_col <- grep("year|yr|Yr", colnames(df), ignore.case = TRUE, value = TRUE)
  seas_col <- grep("seas|season|Seas", colnames(df), ignore.case = TRUE, value = TRUE)
  flt_col <- grep("FltSvy|fleet|index|Flt", colnames(df),
                  ignore.case = TRUE,
                  value = TRUE
  )
  input_SE_col <- grep("_se|se_|Std_in", colnames(df),
                       ignore.case = TRUE,
                       value = TRUE
  ) # catch sample size
  Nsamp_col <- grep("Nsamp", colnames(df),
                    ignore.case = TRUE,
                    value = TRUE
  ) # input sample size
  # sanity checks. should match with 1 (or 0 in some cases) cols. Failing these
  # checks indicate a bug in the code (invalid assuptions of how to match the
  # cols.)
  assertive.properties::assert_is_of_length(yr_col, 1)
  assertive.properties::assert_is_of_length(seas_col, 1)
  assertive.properties::assert_is_of_length(flt_col, 1)
  # b/c only Nsamp or SE should exist for a df
  assertive.base::assert_is_identical_to_true(
    (length(input_SE_col) == 0 & length(Nsamp_col) == 1) |
      (length(input_SE_col) == 1 & length(Nsamp_col) == 0) |
      (length(input_SE_col) == 0 & length(Nsamp_col) == 0)
  )
  # remove equilibrium catch -- do not remove equilibrium catch years in hist
  # df <- df[df[[yr_col]] != -999, ]
  
  sex_col <- grep("Sex|Gender", colnames(df), # sex
                  ignore.case = TRUE,
                  value = TRUE)
  part_col <- grep("part", colnames(df), # partition
                   ignore.case = TRUE,
                   value = TRUE)
  Nsamp_col <- grep("Nsamp", colnames(df), #Nsamp
                    ignore.case = TRUE,
                    value = TRUE)
  # age_err
  ageerr_col <- grep("ageerr", colnames(df), #age err
                     ignore.case = TRUE,
                     value = TRUE)
  Lbin_lo_col <- grep("Lbin_lo", colnames(df), # age err
                      ignore.case = TRUE,
                      value = TRUE)
  Lbin_hi_col <- grep("Lbin_hi", colnames(df), #age err
                      ignore.case = TRUE,
                      value = TRUE)
  #meanbodywt
  type_col <- grep("Type", colnames(df), # type
                   ignore.case = TRUE,
                   value = TRUE) 
  #meanbodywt
  type_col <- grep("Type", colnames(df), # type
                   ignore.case = TRUE,
                   value = TRUE) 
  #mean size at age
  ageerr_col <- grep("Ageerr", colnames(df), # Ageerr
                     ignore.case = TRUE,
                     value = TRUE)
  if(name=="MeanSize_at_Age_obs"){
    n_col <- grep("N_", colnames(df), # N_
                  ignore.case = TRUE,
                  value = TRUE)
    tmp_n <- as.data.frame(rowSums(type.convert(df[,n_col], as.is=TRUE)))
    colnames(tmp_n)<-N_col<-"N_"
    df<-cbind(df,tmp_n)}else{
      N_col<-character(0)
    }
  
  df<-df[,c(yr_col, seas_col, flt_col, type_col, input_SE_col, sex_col,part_col, ageerr_col, Lbin_lo_col, Lbin_hi_col,  Nsamp_col, N_col)] 
  
  #rename to match sample_struct
  if(name=="catch") colnames(df)<- c("Yr", "Seas", "FltSvy", "SE")
  # if(name=="EM2OMcatch_bias") colnames(df) <- c("Yr", "Seas", "FltSvy", "bias")
  if(name=="CPUE") colnames(df) <- c("Yr", "Seas", "FltSvy", "SE")
  if(name=="discard_data") colnames(df)<- c("Yr", "Seas", "FltSvy", "SE")
  # if(name=="EM2OMdiscard_bias") colnames(df)<- c("Yr", "Seas", "FltSvy", "bias")
  if(name=="lencomp") colnames(df)<- c("Yr", "Seas", "FltSvy", "Sex", "Part", "Nsamp")
  if(name=="agecomp") colnames(df)<- c("Yr", "Seas", "FltSvy", "Sex", "Part", "Ageerr", "Lbin_lo", "Lbin_hi", "Nsamp")
  if(name=="meanbodywt") colnames(df)<- c("Yr", "Seas", "FltSvy", "Part", "Type", "Std_in")
  if(name=="MeanSize_at_Age_obs") colnames(df)<- c("Yr", "Seas", "FltSvy", "Sex", "Part", "Ageerr", "N_")
  
  
  return(df) 
}


## COPIED FROM SSMSE GITHUB ###
# to use as reference

# contains the specific MS functions available in the SSMSE package

#' Use EM as the management strategy option.
#' @template EM_out_dir
#' @param init_loop Logical. If this is the first initialization loop of the
#'   MSE, \code{init_loop} should be TRUE. If it is in further loops, it should
#'   be FALSE.
#' @param OM_dat An valid SS data file read in using r4ss. In particular,
#'   this should be sampled data.
#' @template verbose
#' @param nyrs_assess The number of years between assessments. E.g., if an
#'   assessment is conducted every 3 years, put 3 here. A single integer value.
#' @param dat_yrs Which years should be added to the new model? Ignored if
#'  init_loop is TRUE.
#' @template OM_out_dir
#' @template sample_struct
#' @param sample_struct_hist historical sample structure object
#' @template seed
#' @param ... Any additional parameters

#For example, if EM has a positive bias (e.g., EM catch = 1 when true OM catch = 0.5), the EM2OM multiplier should be less than 1 (EM2OM = 0.5)


EnvirEM <- function(EM_out_dir = NULL,
                    init_loop = TRUE,
                    OM_dat,
                    verbose = FALSE,
                    nyrs_assess,
                    dat_yrs,
                    sample_struct = NULL,
                    sample_struct_hist = NULL,
                    seed = NULL,
                    OM_out_dir,
                    ...) {
  
  SSMSE:::check_dir(EM_out_dir)
  # TODO: change this name to make it less ambiguous
  new_datfile_name <- "init_dat.ss"
  # change the name of data file.
  start <- SS_readstarter(file.path(EM_out_dir, "starter.ss"),
                          verbose = FALSE
  )
  
  if (init_loop) {
    
    # copy over raw data file from the OM to EM folder
    SS_writedat(OM_dat,
                file.path(EM_out_dir, new_datfile_name),
                overwrite = TRUE,
                verbose = FALSE
    )
    orig_datfile_name <- start[["datfile"]] # save the original data file name.
    start[["datfile"]] <- new_datfile_name
    start[["seed"]] <- seed
    SS_writestarter(start, file.path(EM_out_dir),
                    verbose = FALSE,
                    overwrite = TRUE
    )
    # make sure the data file has the correct formatting (use existing data
    # file in the EM directory to make sure)??
    # TODO: is this necessary, given we have sample structures? 
    #Yes, now it is necessary because of biased sample struct hist happening within this. 
    new_EM_dat <- biasEM_change_dat(  ###--------------------------------- NEW FUNCTION FOR BIASEM #
      OM_datfile = new_datfile_name,
      EM_datfile = orig_datfile_name,
      EM_dir = EM_out_dir,
      do_checks = TRUE,
      verbose = verbose,
      sample_struct_hist = sample_struct_hist
    )
    ctl <- SS_readctl(file.path(EM_out_dir, start[["ctlfile"]]),
                      datlist = new_EM_dat
    )
    if (ctl[["EmpiricalWAA"]] == 1) {
      message(
        "EM uses weight at age, so copying over wtatage file from OM.",
        "\nNote wtatage data is not sampled."
      )
      file.copy(
        from = file.path(OM_out_dir, "wtatage.ss"),
        to = file.path(EM_out_dir, "wtatage.ss"),
        overwrite = TRUE
      )
    }
    if (!all(ctl[["time_vary_auto_generation"]] == 1)) {
      warning("Turning off autogeneration of time varying lines in the control file of the EM")
      ctl[["time_vary_auto_generation"]] <- rep(1, times = 5)
      r4ss::SS_writectl(ctl, file.path(EM_out_dir, start[["ctlfile"]]),
                        overwrite = TRUE
      )
    }
    
    
  } else {
    
    if (!is.null(sample_struct)) {
      sample_struct_sub <- lapply(sample_struct,
                                  function(df, y) df[df[, 1] %in% y, ],
                                  y = dat_yrs - nyrs_assess
      )
    } else {
      sample_struct_sub <- NULL
    }
    
    new_EM_dat <- add_new_dat_envir( ######## NEW FUNCTION TO BUILD IN CONVERSION FACTOR
      OM_dat = OM_dat,  # why is this OM?  
      EM_datfile = new_datfile_name,
      sample_struct = sample_struct_sub,
      EM_dir = EM_out_dir,
      nyrs_assess = nyrs_assess,
      do_checks = TRUE,
      new_datfile_name = new_datfile_name,
      verbose = verbose
    )
    
    #Increment main recruitment phase end year by number of assessment years
    #So the model can continue to estimate rec devs 
    ctl <- SS_readctl(file.path(EM_out_dir, start[["ctlfile"]]),
                      datlist = new_EM_dat
    )
    ctl$MainRdevYrLast <- ctl$MainRdevYrLast + nyrs_assess
    r4ss::SS_writectl(ctl, file.path(EM_out_dir, start[["ctlfile"]]),
                      overwrite = TRUE
    )
  } # end else not first iteration
  
  # Update SS random seed
  start <- SS_readstarter(file.path(EM_out_dir, "starter.ss"),
                          verbose = FALSE
  )
  start[["seed"]] <- seed
  SS_writestarter(start, file.path(EM_out_dir),
                  verbose = FALSE,
                  overwrite = TRUE
  )
  # manipulate the forecasting file.
  # make sure enough yrs can be forecasted.
  
  fcast <- SS_readforecast(file.path(EM_out_dir, "forecast.ss"),
                           readAll = TRUE,
                           verbose = FALSE
  )
  # check that it can be used in the EM. fleets shoul
  SSMSE:::check_EM_forecast(fcast,
                            n_flts_catch = length(which(new_EM_dat[["fleetinfo"]][, "type"] %in%
                                                          c(1, 2)))
  )
  fcast <- SSMSE:::change_yrs_fcast(fcast,
                                    make_yrs_rel = (init_loop == TRUE),
                                    nyrs_fore = nyrs_assess,
                                    nyrs_increment = nyrs_assess,
                                    mod_styr = new_EM_dat[["styr"]],
                                    mod_endyr = new_EM_dat[["endyr"]]
  )
  # # Need to update the year of forecast assignments where allocations exist
  # if (fcast[["N_allocation_groups"]] > 0) {
  #   fcast$allocation_among_groups$Year <-(new_EM_dat[["endyr"]]+1)
  # }
  
  SS_writeforecast(fcast,
                   dir = EM_out_dir, writeAll = TRUE, overwrite = TRUE,
                   verbose = FALSE
  )
  # given all checks are good, run the EM
  # check convergence (figure out way to error if need convergence)
  # get the future catch using the management strategy used in the SS model.
  run_EM(EM_dir = EM_out_dir, verbose = verbose, check_converged = TRUE)
  
  # get the forecasted catch.
  new_EM_catch_list <- get_EM_catch_df(EM_dir = EM_out_dir, dat = new_EM_dat) 

  ## COnSIDER MAKING A get_OM_catch_df if structure of OM =/= EM. May need to do this now
  # For the simple approach, we can just apply a series of scalars from the EM catch list to create an OM catch list
  new_OM_catch_list = new_EM_catch_list
  
  ## IF FixedCatches==TRUE
  # Will need to be mindful about units -- so far, assuming is presented in same units as historical OM
  # come back to this section if want to use different units
  if(!is.null(sample_struct$FixedCatch)){
    
    # tmp_ss<- sample_struct$FixedCatch[sample_struct$FixedCatch$year==dat_yrs,]
    tmp_ss<- sample_struct$FixedCatch[sample_struct$FixedCatch$year %in% dat_yrs,] # create obj of fixed catches
    colnames(tmp_ss)[4]<- "Fcatch" # rename fixed catches to allow for merge
    
    if(nrow(tmp_ss)>0){
      if(!is.null(new_OM_catch_list$catch)){
        tmp_ss_catch <- tmp_ss[tmp_ss$units!=99,]
        if(nrow(tmp_ss_catch)>0){
          tmp_merge <- base::merge(base::abs(new_OM_catch_list$catch), base::abs(tmp_ss_catch), all.x=TRUE, all.y=FALSE) # merge 
          tmp_merge$catch[which(!is.na(tmp_merge$Fcatch))] <- tmp_merge$Fcatch[which(!is.na(tmp_merge$Fcatch))] # replace fixed catches with 
          tmp_merge <- tmp_merge[base::order(base::abs(tmp_merge$fleet),base::abs(tmp_merge$year),base::abs(tmp_merge$seas)),]
          new_OM_catch_list$catch<-tmp_merge[,c(1:5)] #reorder columns of merged
        }
      }#end if catch exists
      
      if(!is.null(new_OM_catch_list$catch_bio)){
        tmp_ss_catch <- tmp_ss[tmp_ss$units==1,]
        if(nrow(tmp_ss_catch)>0){
          tmp_merge <- base::merge(base::abs(new_OM_catch_list$catch), base::abs(tmp_ss_catch), all.x=TRUE, all.y=FALSE) # merge 
          tmp_merge$catch[which(!is.na(tmp_merge$Fcatch))] <- tmp_merge$Fcatch[which(!is.na(tmp_merge$Fcatch))] # replace fixed catches with
          tmp_merge <- tmp_merge[base::order(base::abs(tmp_merge$fleet),base::abs(tmp_merge$year),base::abs(tmp_merge$seas)),]
          new_OM_catch_list$catch_bio<-tmp_merge[,c(1:5)] #reorder columns of merged
        }
      }else{
        new_OM_catch_list$catch_bio <- NULL }#end if catch_bio exists
      
      if(!is.null(new_OM_catch_list$catch_F)){
        tmp_ss_catch <- tmp_ss[tmp_ss$units==99,]
        if(nrow(tmp_ss_catch)>0){
          tmp_merge <- base::merge(base::abs(new_OM_catch_list$catch_F), base::abs(tmp_ss_catch), all.x=TRUE, all.y=FALSE) # merge 
          tmp_merge$catch[which(!is.na(tmp_merge$Fcatch))] <- tmp_merge$Fcatch[which(!is.na(tmp_merge$Fcatch))] # replace fixed catches with 
          tmp_merge <- tmp_merge[base::order(base::abs(tmp_merge$fleet),base::abs(tmp_merge$year),base::abs(tmp_merge$seas)),]
          # replace the catch_se with the values from FixedCatchEM
          # Create a key to match on
          main_key <- paste(tmp_merge$year, tmp_merge$seas, tmp_merge$fleet, sep = "_")
          update_key <- paste(sample_struct$FixedCatchEM$year, sample_struct$FixedCatchEM$seas, sample_struct$FixedCatchEM$fleet, sep = "_")
          
          # Get the row in df_update for each row in df_main
          match_idx <- match(main_key, update_key)
          
          # Replace catch_se using ifelse inside bracket assignment
          tmp_merge$catch_se <- ifelse(
            !is.na(match_idx) & !is.na(sample_struct$FixedCatchEM$catch_se[match_idx]),
            sample_struct$FixedCatchEM$catch_se[match_idx],
            tmp_merge$catch_se
          )
          new_OM_catch_list$catch_F<-tmp_merge[,c(1:5)] #reorder columns of merged
        }
      }#end if catch exists
    }# end if fixed catches in this mgmt cycle. 
  }# end fixed catches

  # Address EM2OM Catch Bias
  sample_struct$EM2OMcatch_bias <- sample_struct$EM2OMcatch_bias[base::order(base::abs(sample_struct$EM2OMcatch_bias$fleet),base::abs(sample_struct$EM2OMcatch_bias$year),base::abs(sample_struct$EM2OMcatch_bias$seas)),]
  sample_struct$EM2OMcatch_bias <- sample_struct$EM2OMcatch_bias[!duplicated(sample_struct$EM2OMcatch_bias),]
  if(!is.null(new_OM_catch_list$catch)){
    tmp_catch <- base::merge(base::abs(new_OM_catch_list$catch), base::abs(sample_struct$EM2OMcatch_bias), all.x=TRUE, all.y=FALSE)
    tmp_catch <- tmp_catch[base::order(base::abs(tmp_catch$fleet),base::abs(tmp_catch$year),base::abs(tmp_catch$seas)),]
    new_OM_catch_list$catch$catch <- new_OM_catch_list$catch$catch * tmp_catch$bias 
  }
  if(!is.null(new_OM_catch_list$catch_bio)){
    tmp_catch_bio <- base::merge(base::abs(new_OM_catch_list$catch_bio), base::abs(sample_struct$EM2OMcatch_bias), all.x=TRUE)
    tmp_catch_bio <- tmp_catch_bio[base::order(base::abs(tmp_catch_bio$fleet),base::abs(tmp_catch_bio$year),base::abs(tmp_catch_bio$seas)),]
    new_OM_catch_list$catch_bio$catch <- new_OM_catch_list$catch_bio$catch * tmp_catch_bio$bias
  }
  # new_OM_catch_list$catch_bio <- NULL
  if(!is.null(new_OM_catch_list$discards)){
    sample_struct$EM2OMdiscard_bias <- sample_struct$EM2OMdiscard_bias[base::order(base::abs(sample_struct$EM2OMdiscard_bias$Flt),base::abs(sample_struct$EM2OMdiscard_bias$Yr),base::abs(sample_struct$EM2OMdiscard_bias$Seas)),]
    sample_struct$EM2OMdiscard_bias <- sample_struct$EM2OMdiscard_bias[!duplicated(sample_struct$EM2OMdiscard_bias),]
    if(!is.null(new_OM_catch_list$discards)){
      tmp_discards <- base::merge(base::abs(new_OM_catch_list$discards), base::abs(sample_struct$EM2OMdiscard_bias), all.x=TRUE) #need to sort to figure this out
      tmp_discards <- tmp_discards[base::order(base::abs(tmp_discards$Flt),base::abs(tmp_discards$Yr),base::abs(tmp_discards$Seas)),]
      new_OM_catch_list$discards$Discard <- new_OM_catch_list$discards$Discard * tmp_discards$bias 
    }
  }
  
  
  if(2 %in% OM_dat$fleetinfo$type){ # if bycatch fleet -- remove from catch list and keep only bycatch fleets in catch_F list
    # byc_f <- as.numeric(row.names(OM_dat$fleetinfo[which(OM_dat$fleetinfo$type==2),]))
    # new_OM_catch_list$catch<- new_OM_catch_list$catch[new_OM_catch_list$catch$fleet!=byc_f,]
    # new_OM_catch_list$catch_bio<- new_OM_catch_list$catch_bio[new_OM_catch_list$catch$fleet!=byc_f,]
    # new_OM_catch_list$catch_F<- new_OM_catch_list$catch_F[new_OM_catch_list$catch_F$fleet==byc_f,]
    byc_f <- which(OM_dat$fleetinfo$type==2)#as.numeric(row.names(OM_dat$fleetinfo[which(OM_dat$fleetinfo$type==2),]))
    new_OM_catch_list$catch<- new_OM_catch_list$catch[which(!is.element(new_OM_catch_list$catch$fleet,byc_f)),]
    new_OM_catch_list$catch_bio<- new_OM_catch_list$catch_bio[which(!is.element(new_OM_catch_list$catch$fleet,byc_f)),]
    new_OM_catch_list$catch_F<- new_OM_catch_list$catch_F[which(is.element(new_OM_catch_list$catch_F$fleet,byc_f)),]
  } else{
    new_OM_catch_list$catch_F <- NULL
  }
  
  new_catch_list<-new_OM_catch_list
  new_catch_list
  
}



#' Create the sample_struct list for biased catch and discard data
#'
#' Create a sampling structure list using the pattern in a data file and a year
#' range. NAs are added if no pattern is found (and rm_NAs = FALSE). The types
#' of structure that are added to this list (given their presence in the dat file)
#' with their names as called in the list object in parentheses are:
#'  catch (catch), relative indices (CPUE), length composition (lencomp),
#' age composition (agecomp), mean body weight (meanbodywt), and mean size at
#' age (MeanSize_at_Age_obs). Details for creating the sample structure list are
#' available in the [sampling options section of the SSMSE user manual](https://nmfs-fish-tools.github.io/SSMSE/manual/SSMSE.html#sampling-options).
#'
#' @param dat An r4ss list object read in using r4ss::SS_readdat() or the path
#'  (relative or absolute) to an SS data file to read in.
#' @param nyrs Number of years beyond the years included in the dat file to run
#'  the MSE. A single integer value.
#' @param rm_NAs Should all NAs be removed from dataframes? Defaults to FALSE.
#' @export
#' @return A sample_struct list object, where each list element is a dataframe
#'   containing sampling values. If there were no data for the type, NA is
#'   returned for the element.
#' @author Kathryn Doering modified by Cassidy Peterson
#' @examples
#' OM_path <- system.file("extdata", "models", "cod", "ss3.dat", package = "SSMSE")
#' # note there is a warning for lencomp because it does not have a consistent pattern
#' sample_struct <- create_sample_struct(OM_path, nyrs = 20)
#' print(sample_struct)
#' @param FixedCatches T/F defines whether you want to manually specify catches in the future (e.g., for an environmental or bycatch fleet). When switched on, default values to catch in terminal year of OM and historical units. 
#' 
create_sample_struct_envir <- function(dat, nyrs, rm_NAs = FALSE, FixedCatches = FALSE, FixedCatchesEM = FALSE) { ### edited to include EM2OMcatch_bias
  assertive.types::assert_is_a_number(nyrs)
  if (length(dat) == 1 & is.character(dat)) {
    dat <- SS_readdat(dat, verbose = FALSE)
  }
  
  list_name <- c(
    "catch", "EM2OMcatch_bias", "FixedCatch", "CPUE", "discard_data", "EM2OMdiscard_bias", 
    "lencomp", "agecomp", "meanbodywt", "MeanSize_at_Age_obs"
  )
  sample_struct <- lapply(list_name,
                          function(name, dat) {
                            df <- dat[[name]]
                            if (is.null(df)) {
                              return(NA)
                            }
                            # get year, seas, fleet combo, ignoring -999 values.
                            yr_col <- grep("year|yr|Yr", colnames(df), ignore.case = TRUE, value = TRUE)
                            seas_col <- grep("seas|season|Seas", colnames(df), ignore.case = TRUE, value = TRUE)
                            flt_col <- grep("FltSvy|fleet|index|Flt", colnames(df),
                                            ignore.case = TRUE,
                                            value = TRUE
                            )
                            input_SE_col <- grep("_se|se_|Std_in", colnames(df),
                                                 ignore.case = TRUE,
                                                 value = TRUE
                            ) # catch sample size
                            Nsamp_col <- grep("Nsamp", colnames(df),
                                              ignore.case = TRUE,
                                              value = TRUE
                            ) # input sample size
                            # sanity checks. should match with 1 (or 0 in some cases) cols. Failing these
                            # checks indicate a bug in the code (invalid assuptions of how to match the
                            # cols.)
                            assertive.properties::assert_is_of_length(yr_col, 1)
                            assertive.properties::assert_is_of_length(seas_col, 1)
                            assertive.properties::assert_is_of_length(flt_col, 1)
                            # b/c only Nsamp or SE should exist for a df
                            assertive.base::assert_is_identical_to_true(
                              (length(input_SE_col) == 0 & length(Nsamp_col) == 1) |
                                (length(input_SE_col) == 1 & length(Nsamp_col) == 0) |
                                (length(input_SE_col) == 0 & length(Nsamp_col) == 0)
                            )
                            # remove equilibrium catch
                            df <- df[df[[yr_col]] != -999, ]
                            # find combinations of season and fleet in the df.
                            df_combo <- unique(df[, c(seas_col, flt_col), drop = FALSE])
                            fill_vec <- vector(mode = "list", length = nrow(df_combo))
                            for (i in seq_len(nrow(df_combo))) {
                              tmp_seas <- df_combo[i, seas_col]
                              tmp_flt <- df_combo[i, flt_col]
                              tmp_yrs <- df[df[[seas_col]] == tmp_seas &
                                              df[[flt_col]] == tmp_flt, yr_col]
                              tmp_yrs <- as.numeric(unique(tmp_yrs))
                              tmp_yrs <- tmp_yrs[order(tmp_yrs)]
                              
                              
                              # figure out diff between first and second yr.
                              tmp_diff <- tmp_yrs[2] - tmp_yrs[1]
                              # reconstruct the pattern
                              pat <- seq(tmp_yrs[1], by = tmp_diff, length.out = length(tmp_yrs))
                              if (all(!is.na(pat)) && all(pat == tmp_yrs)) { # a pattern was found
                                future_pat <- seq(pat[length(pat)], dat[["endyr"]] + nyrs, by = tmp_diff)
                                future_pat <- future_pat[future_pat > dat[["endyr"]]]
                                if (length(future_pat) > 0) {
                                  future_pat <- data.frame(
                                    Yr = future_pat,
                                    Seas = tmp_seas,
                                    FltSvy = tmp_flt,
                                    stringsAsFactors = FALSE
                                  )
                                } else {
                                  message(
                                    "Pattern found for ", name, ": FltSvy ", tmp_flt,
                                    ", Seas ", tmp_seas, ", but no data to add for the ",
                                    "timeframe specified. Returning NA for Yr in this ",
                                    "dataframe."
                                  )
                                  future_pat <- data.frame(
                                    Yr = NA,
                                    Seas = tmp_seas,
                                    FltSvy = tmp_flt,
                                    stringsAsFactors = FALSE
                                  )
                                }
                              } else {
                                # the pattern was not found
                                warning(
                                  "Pattern not found for ", name, ": FltSvy ", tmp_flt,
                                  ", Seas ", tmp_seas, ". Returning NA for Yr in this dataframe."
                                )
                                future_pat <- data.frame(
                                  Yr = NA,
                                  Seas = tmp_seas,
                                  FltSvy = tmp_flt,
                                  stringsAsFactors = FALSE
                                )
                              }
                              if (name %in% c("lencomp", "agecomp", "MeanSize_at_Age_obs")) {
                                # Sex
                                sex_col <- grep("Sex|Gender", colnames(df),
                                                ignore.case = TRUE,
                                                value = TRUE
                                )
                                tmp_sex <- unique(df[df[[seas_col]] == tmp_seas &
                                                       df[[flt_col]] == tmp_flt, sex_col])
                                if (length(tmp_sex) == 1) {
                                  future_pat[["Sex"]] <- tmp_sex
                                } else {
                                  future_pat[["Sex"]] <- NA
                                }
                              }
                              if (name %in% c("lencomp", "agecomp", "meanbodywt", "MeanSize_at_Age_obs")) {
                                # partition
                                part_col <- grep("part", colnames(df),
                                                 ignore.case = TRUE,
                                                 value = TRUE
                                )
                                tmp_part <- unique(df[df[[seas_col]] == tmp_seas &
                                                        df[[flt_col]] == tmp_flt, part_col])
                                if (length(tmp_part) == 1) {
                                  future_pat[["Part"]] <- tmp_part
                                } else {
                                  future_pat[["Part"]] <- NA
                                }
                              }
                              if (name %in% c("agecomp", "MeanSize_at_Age_obs")) {
                                # Ageerr
                                ageerr_col <- grep("ageerr", colnames(df),
                                                   ignore.case = TRUE,
                                                   value = TRUE
                                )
                                tmp_err <- unique(df[df[[seas_col]] == tmp_seas &
                                                       df[[flt_col]] == tmp_flt, ageerr_col])
                                
                                if (length(tmp_err) == 1) {
                                  future_pat[["Ageerr"]] <- tmp_err
                                } else {
                                  future_pat[["Ageerr"]] <- NA
                                }
                              }
                              if (name == "agecomp") {
                                # Lbin_lo (expect should be -1)
                                tmp_lo <- unique(df[df[[seas_col]] == tmp_seas &
                                                      df[[flt_col]] == tmp_flt, "Lbin_lo"])
                                if (length(tmp_lo) == 1) {
                                  future_pat[["Lbin_lo"]] <- tmp_lo
                                } else {
                                  future_pat[["Lbin_lo"]] <- NA
                                }
                                # Lbin_hi (expect should be -1)
                                tmp_hi <- unique(df[df[[seas_col]] == tmp_seas &
                                                      df[[flt_col]] == tmp_flt, "Lbin_hi"])
                                if (length(tmp_hi) == 1) {
                                  future_pat[["Lbin_hi"]] <- tmp_hi
                                } else {
                                  future_pat[["Lbin_hi"]] <- NA
                                }
                              }
                              if (name == "meanbodywt") {
                                tmp_type <- unique(df[df[[seas_col]] == tmp_seas &
                                                        df[[flt_col]] == tmp_flt, "Type"])
                                if (length(tmp_type) == 1) {
                                  future_pat[["Type"]] <- tmp_type
                                } else {
                                  future_pat[["Type"]] <- NA
                                }
                              }
                              # add sample size, if possible
                              # see if se or Nsamp is the same across years for the seas/flt. If so,
                              # add to the df. If not, add NA's.
                              if (length(input_SE_col) == 1) {
                                tmp_SE <- unique(df[df[[seas_col]] == tmp_seas &
                                                      df[[flt_col]] == tmp_flt &
                                                      df[[yr_col]] != -999, input_SE_col])
                                if (length(tmp_SE) == 1) {
                                  future_pat[["SE"]] <- tmp_SE
                                } else {
                                  future_pat[["SE"]] <- NA
                                  warning("NA included in column SE for ", name, ".")
                                }
                              }
                              if (length(Nsamp_col) == 1) {
                                tmp_Nsamp <- unique(df[df[[seas_col]] == tmp_seas &
                                                         df[[flt_col]] == tmp_flt &
                                                         df[[yr_col]] != -999, Nsamp_col])
                                if (length(tmp_Nsamp) == 1) {
                                  future_pat[["Nsamp"]] <- tmp_Nsamp
                                } else {
                                  future_pat[["Nsamp"]] <- NA
                                  warning("NA included in column Nsamp for ", name, ".")
                                }
                              }
                              if (name == "MeanSize_at_Age_obs") {
                                # Ageerr
                                n_col <- grep("N_", colnames(df),
                                              ignore.case = FALSE,
                                              value = TRUE
                                )
                                tmp_n <- unique(df[df[[seas_col]] == tmp_seas &
                                                     df[[flt_col]] == tmp_flt, n_col])
                                tmp_n <- unlist(tmp_n, use.names = FALSE)
                                tmp_n <- unique(as.numeric(tmp_n))
                                if (length(tmp_n) == 1) {
                                  future_pat[["N_"]] <- tmp_n
                                } else {
                                  future_pat[["N_"]] <- NA
                                }
                              }
                              fill_vec[[i]] <- future_pat
                            }
                            future_pat_all <- do.call("rbind", fill_vec)
                          },
                          dat = dat
  )
  sample_struct <- lapply(
    sample_struct,
    function(x) utils::type.convert(x, as.is = TRUE)
  )
  if (rm_NAs == TRUE) {
    sample_struct <- lapply(
      sample_struct,
      function(x) {
        x <- na.omit(x)
        if (!is.data.frame(x)) {
          x <- NA
        }
        x
      }
    )
  }
  names(sample_struct) <- list_name
  
  ## ADD EM2OMcatch_bias
  sample_struct$EM2OMcatch_bias<- sample_struct$catch
  names(sample_struct$EM2OMcatch_bias)[4] = "bias"
  sample_struct$EM2OMcatch_bias$bias= rep(1, length=nrow(sample_struct$catch))
  
  ## Add FixedCatches
  if(FixedCatches==TRUE){
    sample_struct$FixedCatch <- sample_struct$catch
    sample_struct$FixedCatch$Units <- NA
    names(sample_struct$FixedCatch)[4] = "Catch"
    
    for(f in unique(sample_struct$FixedCatch$FltSvy)){
      
      sample_struct$FixedCatch[sample_struct$FixedCatch$FltSvy==f,]$Catch = rep(dat$catch[dat$catch$year==dat$endyr & dat$catch$fleet==f,]$catch, nyrs)
      if(dat$fleetinfo$type[f] == 1){
        sample_struct$FixedCatch[sample_struct$FixedCatch$FltSvy==f,]$Units = rep(dat$fleetinfo$units[f], nyrs)
      } else{
        sample_struct$FixedCatch[sample_struct$FixedCatch$FltSvy==f,]$Units = rep(99, nyrs)
      }
    }
    sample_struct$FixedCatch
  } else{# end if FixedCatches==TRUE
    FixedCatches <- NULL
  }
  
  if(FixedCatchesEM == TRUE){
    sample_struct$FixedCatchEM <- sample_struct$FixedCatch
    names(sample_struct$FixedCatchEM)[5] = "catch_se"
    
    for(f in unique(sample_struct$FixedCatchEM$FltSvy)){
      sample_struct$FixedCatchEM[sample_struct$FixedCatchEM$FltSvy==f,]$catch_se = rep(dat$catch[dat$catch$year==dat$endyr & dat$catch$fleet==f,]$catch_se, nyrs)
    }
    sample_struct$FixedCatchEM
  } else{# end if FixedCatches==TRUE
    FixedCatchesEM <- NULL
  }
  
  
  ## ADD EM2OMdiscard_bias
  if(!is.null(ncol(sample_struct$discard_data))){
    sample_struct$EM2OMdiscard_bias<- sample_struct$discard_data
    names(sample_struct$EM2OMdiscard_bias)[4] = "bias"
    sample_struct$EM2OMdiscard_bias$bias= rep(1, length=nrow(sample_struct$discard_data))
  } else{
    sample_struct$EM2OMdiscard_bias<-NA
  }
  
  
  
  
  sample_struct
}


#' Create the sample_struct list for historical data

create_sample_struct_hist_envir <- function(dat, rm_NAs = FALSE) { ### edited to include EM2OMcatch_bias
  # assertive.types::assert_is_a_number(nyrs)
  if (length(dat) == 1 & is.character(dat)) {
    dat <- SS_readdat(dat, verbose = FALSE)
  }
  
  list_name <- c(
    "catch",  "CPUE", "discard_data", 
    "lencomp", "agecomp", "meanbodywt", "MeanSize_at_Age_obs"
  )
  
  sample_struct <- lapply(list_name,
                          sample_struct_hist_func,
                          dat = dat
  )
  
  sample_struct <- lapply(
    sample_struct,
    function(x) utils::type.convert(x, as.is = TRUE)
  )
  
  if (rm_NAs == TRUE) {
    sample_struct <- lapply(
      sample_struct,
      function(x) {
        x <- na.omit(x)
        if (!is.data.frame(x)) {
          x <- NA
        }
        x
      }
    )
  }
  names(sample_struct) <- list_name ## END HERE
  
  # ## ADD EM2OMcatch_bias
  # sample_struct$EM2OMcatch_bias<- sample_struct$catch
  # names(sample_struct$EM2OMcatch_bias)[4] = "bias"
  # sample_struct$EM2OMcatch_bias$bias= rep(1, length=nrow(sample_struct$catch))
  # 
  # ## Add FixedCatches
  # if(FixedCatches==TRUE){
  #   sample_struct$FixedCatch <- sample_struct$catch
  #   sample_struct$FixedCatch$Units <- NA
  #   names(sample_struct$FixedCatch)[4] = "Catch"
  #   
  #   for(f in unique(sample_struct$FixedCatch$FltSvy)){
  #     sample_struct$FixedCatch[sample_struct$FixedCatch$FltSvy==f,]$Catch = rep(dat$catch[dat$catch$year==dat$endyr & dat$catch$fleet==f,]$catch, nyrs)
  #     sample_struct$FixedCatch[sample_struct$FixedCatch$FltSvy==f,]$Units = rep(dat$fleetinfo$units[f], nyrs)
  #   }
  #   sample_struct$FixedCatch
  # } else{# end if FixedCatches==TRUE
  #   FixedCatches <- NULL
  # }
  # 
  # ## ADD EM2OMdiscard_bias
  # if(!is.null(ncol(sample_struct$discard_data))){
  #   sample_struct$EM2OMdiscard_bias<- sample_struct$discard_data
  #   names(sample_struct$EM2OMdiscard_bias)[4] = "bias"
  #   sample_struct$EM2OMdiscard_bias$bias= rep(1, length=nrow(sample_struct$discard_data))
  # } else{
  #   sample_struct$EM2OMdiscard_bias<-NA
  # }
  # 
  
  
  
  sample_struct
} 

# catch: c("Yr", "Seas", "FltSvy", "SE")
# EM2OMcatch_bias: c("Yr", "Seas", "FltSvy", "bias")
# CPUE: c("Yr", "Seas", "FltSvy", "SE")
# discard_data: c("Yr", "Seas", "FltSvy", "SE")
# EM2OMdiscard_bias: c("Yr", "Seas", "FltSvy", "bias")
# FixedCatch: c("Yr", "Seas", "FltSvy", "Catch", "Units")
# lencomp: c("Yr", "Seas", "FltSvy", "Sex", "Part", "Nsamp")
# agecomp: c("Yr", "Seas", "FltSvy", "Sex", "Part", "Ageerr", "Lbin_lo", "Lbin_hi", "Nsamp")
# meanbodywt: c("Yr", "Seas", "FltSvy", "Part", "Type", "Std_in")
# MeanSize_at_Age_obs: c("Yr", "Seas", "FltSvy", "Sex", "Part", "Ageerr", "N_")


#' Create the sample_struct list for historical data

create_sample_struct_hist_biased <- function(dat, rm_NAs = FALSE, FixedCatches = FALSE) { ### edited to include EM2OMcatch_bias
  # assertive.types::assert_is_a_number(nyrs)
  if (length(dat) == 1 & is.character(dat)) {
    dat <- SS_readdat(dat, verbose = FALSE)
  }
  
  list_name <- c(
    "catch", "EM2OMcatch_bias", "CPUE", "discard_data", "EM2OMdiscard_bias", 
    "lencomp", "agecomp", "meanbodywt", "MeanSize_at_Age_obs"
  )
  

  sample_struct <- lapply(list_name,
                          sample_struct_hist_func,
                          dat = dat
  )
  
  sample_struct <- lapply(
    sample_struct,
    function(x) utils::type.convert(x, as.is = TRUE)
  )
  
  if (rm_NAs == TRUE) {
    sample_struct <- lapply(
      sample_struct,
      function(x) {
        x <- na.omit(x)
        if (!is.data.frame(x)) {
          x <- NA
        }
        x
      }
    )
  }
  names(sample_struct) <- list_name ## END HERE
  
  ## ADD EM2OMcatch_bias
  sample_struct$EM2OMcatch_bias<- sample_struct$catch
  names(sample_struct$EM2OMcatch_bias)[4] = "bias"
  sample_struct$EM2OMcatch_bias$bias= rep(1, length=nrow(sample_struct$catch))
  
  ## Add FixedCatches
  if(FixedCatches==TRUE){
    sample_struct$FixedCatch <- sample_struct$catch
    sample_struct$FixedCatch$Units <- NA
    names(sample_struct$FixedCatch)[4] = "Catch"

    for(f in unique(sample_struct$FixedCatch$FltSvy)){
      sample_struct$FixedCatch[sample_struct$FixedCatch$FltSvy==f,]$Catch = rep(dat$catch[dat$catch$year==dat$endyr & dat$catch$fleet==f,]$catch, nyrs)
      sample_struct$FixedCatch[sample_struct$FixedCatch$FltSvy==f,]$Units = rep(dat$fleetinfo$units[f], nyrs)
    }
    sample_struct$FixedCatch
  } else{# end if FixedCatches==TRUE
    FixedCatches <- NULL
  }
  
  ## ADD EM2OMdiscard_bias
  if(!is.null(ncol(sample_struct$discard_data))){
    sample_struct$EM2OMdiscard_bias<- sample_struct$discard_data
    names(sample_struct$EM2OMdiscard_bias)[4] = "bias"
    sample_struct$EM2OMdiscard_bias$bias= rep(1, length=nrow(sample_struct$discard_data))
  } else{
    sample_struct$EM2OMdiscard_bias<-NA
  }
  
  
  sample_struct
} 




#' Change dataset from OM into format for EM
#' @param OM_datfile Filename of the datfile produced by the OM within the
#'  EM_dir.
#' @param EM_datfile Filename of the datfile from the original EM within the
#'  EM_dir.
#' @param EM_dir Absolute or relative path to the Estimation model directory.
#' @param do_checks Should checks on the data be performed? Defaults to TRUE.
#' @template verbose
#' @param sample_struct_hist historical sample structure object
#' @author Kathryn Doering modified by Cassidy Peterson
#' @importFrom r4ss SS_readstarter SS_readdat SS_writedat SS_writestarter
#' @return the new EM data file. Side effect is saving over the OM_dat file in
#'   EM_dir.
#' @examples
#' \dontrun{
#' # TODO: Add example
#' }
biasEM_change_dat <- function(OM_datfile,
                              EM_datfile,
                              EM_dir,
                              do_checks = TRUE,
                              verbose = FALSE,
                              sample_struct_hist = NULL) {
  
  EM_dir <- normalizePath(EM_dir)
  
  # checks
  assertive.types::assert_is_a_string(OM_datfile)
  assertive.types::assert_is_a_string(EM_dir)
  check_dir(EM_dir)
  assertive.types::assert_is_a_bool(do_checks)
  assertive.types::assert_is_a_bool(verbose)
  
  # read in the dat files
  EM_dat <- SS_readdat(file.path(EM_dir, EM_datfile), verbose = FALSE)
  OM_dat <- SS_readdat(file.path(EM_dir, OM_datfile), verbose = FALSE)
  
  # remove extra years of data in the OM data file to create EM.
  new_EM_dat <- get_EM_dat(
    OM_dat = OM_dat, EM_dat = EM_dat,
    do_checks = do_checks
  )
  
  # Add bias-correction in historical period: 
  if(!is.null(sample_struct_hist$EM2OMcatch_bias)){
    new_EM_dat$catch$catch<-new_EM_dat$catch$catch/sample_struct_hist$EM2OMcatch_bias$bias
    if(new_EM_dat$N_discard_fleets>0){
      new_EM_dat$discard_data$Discard<-new_EM_dat$discard_data$Discard / sample_struct_hist$EM2OMdiscard_bias$bias
    }
  }# end if sample_struct_hist exists
  # TODO: Add sample_struct_hist code here to alter the dat file.  
  
  # write out the modified files that can be used in future EM run
  SS_writedat(new_EM_dat, file.path(EM_dir, OM_datfile),
              verbose = FALSE,
              overwrite = TRUE
  )
  
  return(new_EM_dat)
}



