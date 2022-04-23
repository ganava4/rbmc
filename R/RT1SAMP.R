##' One-sample or Paired comparison randomization test
##'
##' RT1SAMP carries out the Fisher randomization test or the paired comparison
##' randomization test. The permutational one-sample randomization test does not
##' require normality of the distribution. A randomization confidence interval
##' for the mean can also be determined.
##'
##' The procedure (1) calculates the mean scores for the difference between the
##' observed data \code{x} and the value of mu under the null hypothesis
##' \code{mu}. The user has the option to enumerate all possible randomizations
##' or to sample the randomization distribution. In this latter case, the
##' function \code{sampleRD} is invoked. In addition, it is possible to
##' produce an approximate 95% / 99% confidence interval by interpolation
##' (useful for simulations).  If confidence limits are needed, the number of
##' randomizations should be a large enough number (probably 4999 or more). The
##' upper percentage points (percentages of randomization differences greater
##' than or equal to the observed \code{x - mu} means) can be determined for a
##' range of trial values for L and U, which, when subtracted from the \code{x},
##' just avoid giving a significant difference between the two sample means. The
##' upper percentage points (percentages of randomization differences greater
##' than or equal to the observed difference between \code{x} and \code{mu}) can
##' be determined for a range of trial values for L, and linear interpolation is
##' used to determine the value of L to be substracted from the \code{x}
##' measurements producing a difference between \code{x} and \code{mu} that is
##' on the borderline of being significantly large at about the 2.5\% or 0.5\%
##' level. Analogously, the lower percentage points (the percentages of
##' randomization differences less than or equal to the observed difference
##' \code{x - mu} can be determined for some trial values of U. Again, linear
##' interpolation is used to determine the value of U to be substracted from the
##' \code{x - mu} measurements producing a difference that is on the borderline
##' of being significantly small at about the 2.5\% or 0.5\%level.
##'
##' \code{seed} is a way to call the \code{\link[base]{set.seed}} function,
##' "the recommended way to specify seeds" in random number generation.
##'
##' The function \code{summary.RT} is used to obtain and print a summary of the
##' results, and a \code{plot.RT} method is available for displaying the
##' randomization distribution of mean differences.
##'
##' @param x a (non-empty) numeric vector of data values.
##' @param mu an scalar; it is the value of \code{mu} under the null hypothesis
##' @param NRAND numeric; the number of randomizations (permutations). It is
##' ignored when \code{complete=TRUE}.
##' @param alt a character string specifying the alternative hypothesis, must
##' be one of "\code{two}" (two-sided, the default), "\code{greater}" or
##' "\code{less}".
##' @param complete a logical variable; when \code{TRUE}, all possible
##' randomizations are enumerated. It is ignored when \code{CI=TRUE}.
##' @param CI a logical variable indicating whether approximate 95\% and 99\%
##' randomization confidence intervals will be calculated (\code{TRUE}) or not
##' (\code{FALSE}).
##' @param silent a logical variable indicating whether calculation results are
##' printed to the R console (\code{silent = FALSE}). If \code{TRUE} then
##' calculation results are not printed to the R console (useful for
##' simulations)
##' @param seed a single value, interpreted as an integer, or \code{NULL} (see
##' "Details").
##' @return The function returns a \code{RT} result object (list)
##' @author Jorge Navarro-Alberto
##' @seealso \code{\link[rbmc]{sampleRD}}, \code{\link[rbmc]{count}},
##' \code{\link[rbmc]{summary.RT}} and the main plotting function
##' \code{\link[rbmc]{plot.RT}}
##' @references Manly, B.F.J. and Navarro-Alberto, J.A. (2021) Randomization,
##' Bootstrap and Monte Carlo Methods in Biology. 4th Edition. Chapman and
##' Hall/ CRC Press.
##' @examples
##'
##' # Example in Manly and Navarro Alberto (2021), Section 1.5
##' Dif <- cornheight[,"Cross"] - cornheight[,"Self"]
##' corn1s.c <- RT1SAMP(Dif, complete=TRUE) # Complete enumerations
##' #
##' corn1s.s <- RT1SAMP(Dif, CI=TRUE) # Sampled distribution, CI for x - mu.
##'
##'
##' @export

RT1SAMP <- function(x, mu=0, NRAND = 4999, alt="two",
                    complete=FALSE, CI = FALSE, silent = FALSE, seed = NULL)
# This function computes the Fisher randomization test
# or the paired comparison randomization test. The
# permutational one-sample randomization test does not require
# normality of the distributions.
# Parameters of the function.
# x: vector of observations
# mu: Value of mu under the null hypothesis
# NRAND = number of permutations (default value: 4999)
# complete = TRUE: enumerates all possible randomizations
# complete = FALSE: samples the randomization distribution
# CI = FALSE: Does not calculate an approximate 95% / 99% confidence interval
# CI = TRUE: Calculates an approximate 95% / 99% confidence interval by
#      interpolation
# silent = FALSE: calculation results are printed to the R console.
# silent = TRUE: calculation results are not printed to
# the R console (useful for simulations).  Values returned:
{
  rt1 <- list()
  CI9599 <- NULL
  # Compute the sumD-statistic for the observed data
  sumD0 <- sum(x - mu)
  abs.sumD0 <- abs(sumD0)
  n <- length(x)
  RDtype <- ifelse(complete == TRUE, "Complete enumerations",
                   "Randomization distribution was sampled")
  if (complete == TRUE) {
  if (!silent) {
    if (alt == "less")
      cat("Lower tailed test\n")
    if (alt == "greater")
      cat("Upper tailed test\n")
    if (alt == "two")
      cat("Two tailed test\n")
  }
# Function count() performs the test using complete enumeration
  complete.enum <- count(x, mu)
  P <- complete.enum$P
  NRAND <- complete.enum$nperm
  type <- ifelse(alt == "less", 1, ifelse(alt == "greater", 2, 3))
  if (!silent) {
    cat("Complete enumeration of distribution\n")
    cat("Total number of possibilities =", complete.enum$nperm, "\n")
    cat("Observed sum =", formatC(complete.enum$Sobs, digits = 2,
                                width=7, format = "f"), "\n")
    cat("Sig. level (%) =", formatC(complete.enum$P[type], digits = 2,
                                    width = 5, format = "f"), "\n")
  }
}
if (complete == FALSE) {
  if (!silent) {
    if (alt == "less")
      cat("Lower tailed test\n")
    if (alt == "greater")
      cat("Upper tailed test\n")
    if (alt == "two")
      cat("Two tailed test\n")
  }
# Function sampleRD() performs the test, by sampling the randomization
# distribution
  sample.dist <- sampleRD(x, mu, NRAND)
  P <- sample.dist$P
  type <- ifelse(alt == "less", 1, ifelse(alt == "greater", 2, 3))
  if (!silent) {
    cat("The randomization distribution has been sampled\n")
    cat("Randomizations (observed included) =", sample.dist$nperm, "\n")
    cat("Observed sum", formatC(sample.dist$Sobs, digits = 2,
                                width=7, format = "f"), "\n")
    cat("Sig. level (%) =", formatC(sample.dist$P[type], digits = 2,
                                    width = 5, format = "f"), "\n")
  }
}
  # If CI=TRUE, approximate CIs are calculated
  if (CI == TRUE) {
    # Find limits to try the Delta
    (VMEAN <- mean(x))
    (SE <- sd(x)/sqrt(n))
    (DMIN <- VMEAN - 3.5 * SE)
    (DINC <- 0.2 * SE)
    # Try Delta values
    DELTA <- numeric(36)
    Q <- matrix(0, nrow = 36, ncol = 2)
    for (i in 1:36) {
      DELTA[i] <- DMIN + (i - 1) * DINC
      x.CI <- x - DELTA[i]
      P.CI <- rep(NA, 2)
      if (complete == TRUE) {
        comp.Delta <- count(x.CI, mu)
        P.CI <- comp.Delta$P[1:2]
      }
      else {
        samp.Delta <- sampleRD(x.CI, mu, NRAND)
        P.CI <- samp.Delta$P[1:2]
      }
      Q[i, 1:2] <- P.CI
      print(c(DELTA[i], Q[i, 1], Q[i, 2]))
      DELTA0 <- DELTA[i] - DINC
      if (i == 1) {next}
      else {
        if ((Q[i, 2] >= 2.5) && (Q[i - 1, 2] < 2.5))
          DL25 <- approx(c(Q[i - 1, 2], Q[i, 2]), c(DELTA0, DELTA[i]),
                         xout = 2.5)$y
          else if ((Q[i, 1] < 2.5) && (Q[i - 1, 1] >= 2.5))
                    DU25 <- approx(c(Q[i - 1, 1], Q[i, 1]), c(DELTA0, DELTA[i]),
                                          xout = 2.5)$y
          else if ((Q[i, 2] >= 0.5) && (Q[i - 1, 2] < 0.5))
                    DL05 <- approx(c(Q[i, 2], Q[i - 1, 2]), c(DELTA0, DELTA[i]),
                                          xout = 0.5)$y
          else if ((Q[i, 1] < 0.5) && (Q[i - 1, 1] >= 0.5))
                    DU05 <- approx(c(Q[i - 1, 1], Q[i, 1]), c(DELTA0, DELTA[i]),
                                          xout = 0.5)$y
      }
    }
    cat("  Approximate confidence limits by interpolation\n")
    CI9599 <- c(DL25, DU25, DL05, DU05)
    print(CI9599)
 }
 rt1 <- list(Sobs = sumD0, name="Fisher (or paired sample) randomization test",
              RDtype = RDtype, p.perm = P, NRAND = NRAND, alt = alt,
              CI9599 = CI9599)
 class(rt1) <- "RT"
 return(rt1)
}
