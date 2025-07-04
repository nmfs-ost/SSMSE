% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runSSMSE.R
\name{run_SSMSE_iter}
\alias{run_SSMSE_iter}
\title{Run one iteration of an MSE using SS3 OM}
\usage{
run_SSMSE_iter(
  out_dir = NULL,
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
  verbose = FALSE
)
}
\arguments{
\item{out_dir}{The directory to which to write output. IF NULL, will default
to the working directory.}

\item{OM_name}{Name of the operating model (OM).}

\item{OM_in_dir}{Relative or absolute path to the operating model, if using a
model outside of the SSMSE package. Should be a string.}

\item{EM_name}{Name of a valid Stock Synthesis stock assessment model to use
as an EM. If the value of EM_name is NULL and  \code{MS = "EM"}, then SSMSE
will look for the estimation model in the path specified in EM_in_dir.
or \code{NULL}.}

\item{EM_in_dir}{Relative or absolute path to the estimation model,}

\item{run_EM_last_yr}{Should the MS be implemented to get future catch if the
last year is an assessment year? TRUE/FALSE option, so the same for all
scenarios and iterations. Defaults to FALSE.}

\item{MS}{The management strategy to use. Current options
are:  \code{"last_yr_catch"} which uses the previous year's catch;
\code{"no_catch"} which uses 0 catch; \code{"EM"} which uses an stock
synthesis model as the estimation method and the management strategy as
defined in the forecast file of the stock synthesis estimation method;
\code{"Interim"} to modify catch based on survey predictions between
assessments. Users can also specify their own management strategies as a
function. For example, if the function is
called "my_ms" then the user should specify MS = "my_ms" and specify the
path to the file containing the function in \code{custom_MS_source}.}

\item{custom_MS_source}{A file location with the source code for any custom MS
models to be used in the simulation. SSMSE will source this file which should
contain a function/s whose name/s match each custom MS models included in MS_vec.
To learn more about using custom management strategies, see the
\href{https://nmfs-ost.github.io/SSMSE/manual/custom-MS.html}{using a Custom Management Stratey/Procedure section}
in the \href{https://nmfs-ost.github.io/SSMSE/manual}{SSMSE User Manual}.}

\item{use_SS_boot}{Should a bootstrapped data set generated by SS3 be used?
Defaults to TRUE.}

\item{nyrs}{Number of years beyond the years included in the OM to run the
MSE. A single integer value.}

\item{nyrs_assess}{The number of years between assessments. E.g., if an
assessment is conducted every 3 years, put 3 here. A single integer value.
(NOTE: we could make this more flexible by instead reading in a vector of
assessment years, so users could specify irregular numbers of yrs between
assessments.)}

\item{nyrs_lag}{Number of years of lag in obtaining data (i.e., the number of years
post EM assessment end yr before advice can be implemented). Defaults to 0.}

\item{nscen}{Which scenario is this. Integer value >=1}

\item{scen_name}{Name of the scenario. The directory containing all the model
runs the scenario will be stored within a folder of this name.}

\item{niter}{The iteration number, which is also the name of the folder the
results will be written to.}

\item{iter_seed}{List containing fixed seeds for this iteration.}

\item{sample_struct}{A optional list including which years, seasons, and fleets
should be  added from the OM into the EM for different types of data.
If NULL, the data structure will try to be inferred from the pattern found
for each of the datatypes within the EM datafiles. Include this structure
for the number of years to extend the model out. Note that the data should
be specified using the list component names and column names as in would be
used in \code{r4ss::SS_readdat()}. The \code{run_SSMSE_iter} function
examples give an example of what this structure should be. Running the
function \code{create_sample_struct()} will also produce a sample_struct
object in the correct form. Can be NULL only when MS is not EM.}

\item{future_om_list}{An optional list of lists including changes that should
be made after the end year of the input model. Each first level list element
outlines 1 change to be made to the operating model. To see an example, try
running \code{\link{create_future_om_list}}. Defaults to NULL, which implies
that the model will be extended forward in time assuming the original model
structure.}

\item{sample_struct_hist}{An optional list including which years should be
sampled for the historical period for the data generated from the OM. If
this is left as NULL, then the same sampling scheme will be used as in the
OM's data file. If it is not NULL, then each year.}

\item{sample_catch}{Should catch be sampled or fixed at the OM values?
Defaults to FALSE.}

\item{interim_struct}{A optional list of parameters to control an interim
assessment with an example structure below, where Beta=a positive value
that is inversely proportional to risk, MA_years= the number of years to
average index observations of when calculating deviations, assess_freq=the
number of years between full assessments during with an interim assessment
will happen every year, and Index_weights is a vector of length n indexes
that weights all indexes for multi index inference.

\if{html}{\out{<div class="sourceCode">}}\preformatted{interim_struct<-list(Beta = 1,
                     MA_years = 3,
                     assess_freq = 5,
                     Index_weights = rep(1, max(ref_index[, 3])))
}\if{html}{\out{</div>}}}

\item{n_F_search_loops}{Number of times to try to find an F that achieves the
catches input in the OM. Defaults to 20.}

\item{tolerance_F_search}{How far apart the input catch and achieved
catch can be in tried to find an F that achieves the catch input in the
OM. Defaults to 0.001.}

\item{verbose}{Want verbose output? Defaults to FALSE.}
}
\description{
High level function to run 1 iteration of a scenario for a management
strategy evaluation using Stock Synthesis as the Operating model.
}
\examples{
\dontrun{
# Create a temporary folder for the output
temp_path <- file.path(tempdir(), "run_SSMSE_iter-ex")
dir.create(temp_path)

# run 1 iteration and 1 scenario of SSMSE
run_SSMSE_iter(
  OM_name = "cod",
  MS = "no_catch",
  out_dir = temp_path,
  nyrs = 6,
  nyrs_assess = 3
)
unlink(file.path(temp_path, "1"), recursive = TRUE)
# run 1 iteration and 1 scenario of SSMSE using an EM.
run_SSMSE_iter(
  OM_name = "cod",
  MS = "EM",
  out_dir = temp_path,
  EM_name = "cod",
  nyrs = 6,
  nyrs_assess = 3,
  sample_struct = list(
    catch = data.frame(Yr = 101:106, Seas = 1, FltSvy = 1, SE = 0.05),
    CPUE = data.frame(Yr = c(102, 105), Seas = 7, FltSvy = 2, SE = 0.01),
    lencomp = data.frame(
      Yr = c(102, 105), Seas = 1, FltSvy = 1,
      Sex = 0, Part = 0, Nsamp = 100
    ),
    agecomp = data.frame(
      Yr = c(102, 105), Seas = 1, FltSvy = 2,
      Sex = 0, Part = 0, Ageerr = 1,
      Lbin_lo = -1, Lbin_hi = -1, Nsamp = 50
    )
  )
)
unlink(temp_path, recursive = TRUE)
}
}
\author{
Kathryn Doering & Nathan Vaughan
}
