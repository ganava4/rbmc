#' Summarize one and two sample randomization tests
#'
#' Summarize the results produced by RT1SAMP and RT2SAMP
#'
#' @param x an object of class RT
#' @seealso \code{\link[rbmc]{RT1SAMP}}, \code{\link[rbmc]{RT2SAMP}}, and the
#' main plotting function \code{\link[rbmc]{plot.RT}}
#' @references Manly, B.F.J. and Navarro-Alberto, J.A. (2021) Randomization,
#' Bootstrap and Monte Carlo Methods in Biology. 4th Edition. Chapman and
#' Hall/ CRC Press.
#' @examples
#' # Example in Manly and Navarro Alberto (2021), Section 1.5
#' Dif <- cornheight[,"Cross"] - cornheight[,"Self"]
#' corn1s.c <- RT1SAMP(Dif, complete=TRUE) # Complete enumeration
#' summary.RT(corn1s.c)
#'
#' @export summary.RT
summary.RT <- function(x) {
  stopifnot(inherits(x, "RT"))
  if (x$name == "Two sample RT for mean difference") {
  cat(" Results of Randomization Testing of the Mean Difference\n")
  cat("   Observed Mean Difference =", formatC(x$D.obs, digits = 3, width = 6,
                                               format = "f"), "\n")
  cat("   With", x$NRAND, "permutations: \n")
  cat("   Less than or equal to observed = ",
          formatC(x$p.perm[1], digits = 3, width = 7, format = "f"), "% \n")
  cat("   Greater than or equal to observed = ",
          formatC(x$p.perm[2], digits = 3, width = 7, format = "f"), "% \n")
  cat("   As far as or further from zero than observed = ",
          formatC(x$p.perm[3], digits = 3, width = 7, format = "f"), "% \n\n")
  if (!is.null(x$CI9599)) {
   cat(" Approximate confidence limits of mean difference (by interpolation)\n")
   cat("   95%:", formatC(x$CI9599[1], digits = 3, width = 7, format = "f"),
        " to  ", formatC(x$CI9599[2], digits = 3, width = 7, format = "f"),
        ";  99%:", formatC(x$CI9599[3], digits = 3, width = 7, format = "f"),
        " to  ", formatC(x$CI9599[4], digits = 3, width = 7, format = "f"))
    }
  }
  else
  if (x$name == "Fisher (or paired sample) randomization test") {
    cat("Results of Fisher (or paired sample) randomization test\n")
    cat(x$RDtype,"\n")
    cat("Observed sum =", formatC(x$Sobs, digits = 2,
                                width=7, format = "f"), "\n")
    cat("   With", x$NRAND, "permutations: \n")
    cat("   Less than or equal to observed = ",
        formatC(x$p.perm[1], digits = 3, width = 7, format = "f"), "% \n")
    cat("   Greater than or equal to observed = ",
        formatC(x$p.perm[2], digits = 3, width = 7, format = "f"), "% \n")
    cat("   As far as or further from zero than observed = ",
        formatC(x$p.perm[3], digits = 3, width = 7, format = "f"), "% \n\n")
    if (!is.null(x$CI9599)) {
      cat(" Approximate confidence limits (by interpolation)\n")
      cat("   95%:", formatC(x$CI9599[1], digits = 3, width = 7, format = "f"),
          " to  ", formatC(x$CI9599[2], digits = 3, width = 7, format = "f"),
          ";  99%:", formatC(x$CI9599[3], digits = 3, width = 7, format = "f"),
          " to  ", formatC(x$CI9599[4], digits = 3, width = 7, format = "f"))
    }
  }
}
