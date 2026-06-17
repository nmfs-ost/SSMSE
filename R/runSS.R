# Utility functions for running SS, taken from ss3sim.
# https://github.com/ss3sim/ss3sim

#' Get SS3 binary/executable location in package
#'
#' Get the binary/executable location in the package SSMSE. This function
#' is from \href{https://github.com/ss3sim/ss3sim}{ss3sim}.
#'
#' @param bin_name Name of SS3 binary, defaults to "ss3"
#' @return The path to an SS3 binary. If using the GitHub version of the
#'   package, this will be an internal binary. Otherwise, this function
#'   will search for a version of the binary in your path. See the
#'   ss3sim vignette.
#' @export
#' @examples
#' \dontrun{
#' get_bin()
#' }
get_bin <- function(bin_name = "ss3") {
  # code inspiration from glmmADMB package:
  if (.Platform[["OS.type"]] == "windows") {
    platform <- "Windows64"
    bin_name <- paste0(bin_name, ".exe")
    bit <- R.version[["arch"]]
    if (grepl("3", bit)) {
      if (!grepl("86", bit)) {
        platform <- "Windows32"
        warning(
          "SS3 binary is not available for 32-bit ",
          .Platform[["OS.type"]], " within the package. ",
          "You must have an appropriate SS3 binary in your path. ",
          "See the ss3sim vignette."
        )
      }
    }
  } else {
    if (substr(R.version[["os"]], 1, 6) == "darwin" && R.version[["arch"]] == "x86_64") {
      platform <- "MacOS"
    } else {
      if (substr(R.version[["os"]], 1, 6) == "darwin" && R.version[["arch"]] == "aarch64") {
        platform <- "MacOS_arm64"
      } else {
        if (R.version[["os"]] == "linux-gnu") {
          platform <- "Linux64"
        } else {
          warning(
            "SS3 binary is not available for OS ", R.version[["os"]], " ", R.version[["arch"]],
            " within the package. You must have an appropriate SS3 binary in your ",
            "path. See the ss3sim vignette."
          )
        }
      }
    }
  }
  loc <- system.file("bin", package = "SSMSE")
  if (loc != "") { # we found binaries in the package
    bin <- file.path(loc, platform, bin_name)
    if (!file.exists(bin)) bin <- ""
  } else {
    bin <- ""
  }
  if (bin == "") { # resort to binaries in path
    bin <- Sys.which(bin_name)[[1]]
    if (bin == "") {
      stop(
        "The expected SS3 executable, ", bin_name, ", was not found in your",
        " path. See the ss3sim vignette and ?ss3sim::run_ss3model for ",
        "instructions."
      )
    }
  }
  if (grepl("[[:space:]]", bin)) {
    bin <- shQuote(bin)
  }
  bin
}
