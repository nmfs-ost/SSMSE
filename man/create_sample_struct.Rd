% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_struct.R
\name{create_sample_struct}
\alias{create_sample_struct}
\title{Create the sample_struct list}
\usage{
create_sample_struct(dat, nyrs, rm_NAs = FALSE)
}
\arguments{
\item{dat}{An r4ss list object read in using r4ss::SS_readdat() or the path
(relative or absolute) to an SS3 data file to read in.}

\item{nyrs}{Number of years beyond the years included in the dat file to run
the MSE. A single integer value.}

\item{rm_NAs}{Should all NAs be removed from dataframes? Defaults to FALSE.}
}
\value{
A sample_struct list object, where each list element is a dataframe
containing sampling values. If there were no data for the type, NA is
returned for the element.
}
\description{
Create a sampling structure list using the pattern in a data file and a year
range. NAs are added if no pattern is found (and rm_NAs = FALSE). The types
of structure that are added to this list (given their presence in the dat file)
with their names as called in the list object in parentheses are:
catch (catch), relative indices (CPUE), length composition (lencomp),
age composition (agecomp), mean body weight (meanbodywt), and mean size at
age (MeanSize_at_Age_obs). Details for creating the sample structure list are
available in the \href{https://nmfs-ost.github.io/SSMSE/manual/SSMSE.html#sampling-options}{sampling options section of the SSMSE user manual}.
}
\examples{
OM_path <- system.file("extdata", "models", "cod", "ss3.dat", package = "SSMSE")
# note there is a warning for lencomp because it does not have a consistent pattern
sample_struct <- create_sample_struct(OM_path, nyrs = 20)
print(sample_struct)
}
\author{
Kathryn Doering
}
