#' Randomization Test of the Mean Difference between Two Samples
#'
#' RT2SAMP carries out a two sample randomization test on the mean difference
#' between two samples. A randomization confidence interval for the mean
#' difference between two source populations can also be determined.
#'
#' The procedure (1) calculates the mean scores for vectors \code{x} and
#' \code{y}, and the difference \code{D[0]} between these; then (2)
#' \code{length(x)} and \code{length(y)} observations are randomly reallocated
#' to the first and second group, respectively, using the \code{\link{sample}}
#' function. Step (2) is repeated \code{NRAND} times to find the randomization
#' distribution of \code{D} differences. \code{alt = "greater"} is the
#' alternative that x has a larger mean than y.
#'
#' The randomization test allows confidence intervals to be placed on
#' treatment effects (when \code{CI = TRUE}), as described in Manly and
#' Navarro (2021, Section 1.4). If confidence limits are needed, \code{NRAND}
#' should be a large enough number (probably 4999 or more). The upper
#' percentage points (percentages of randomization differences greater than or
#' equal to the observed between \code{x} and \code{y} means) can be
#' determined for a range of trial values for L and U, which, when subtracted
#' from the \code{x}, just avoid giving a significant difference between the
#' two sample means. The upper percentage points (percentages of randomization
#' differences greater than or equal to the observed difference between
#' \code{x} and \code{y}) can be determined for a range of trial values for L,
#' and linear interpolation is used to determine the value of L to be
#' substracted from the \code{x} measurements producing a difference between
#' \code{x} and \code{y} that is on the borderline of being significantly
#' large at about the 2.5% or 0.5%. Analogously, the lower percentage
#' points (the percentages of randomization differences less than or equal to
#' the observed difference between \code{x} and \code{y} means) can be
#' determined for some trial values of U. Again, linear interpolation is used
#' to determine the value of U to be substracted from the \code{x}
#' measurements producing a difference between \code{x} and \code{y} that is
#' on the borderline of being significantly small at about the 2.5% or
#' 0.5% level.
#'
#' \code{seed} is a way to call the \code{\link[base]{set.seed}} function,
#' "the recommended way to specify seeds" in random number generation.
#'
#' The function \code{\link[rbmc]{summary.RT}} is used to obtain and print a summary of the
#' results, and a \code{plot.RT} method is available for displaying the
#' randomization distribution of mean differences.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y a (non-empty) numeric vector of data values.
#' @param NRAND numeric; the number of randomizations (permutations).
#' @param alt a character string specifying the alternative hypothesis, must
#' be one of "\code{two}" (two-sided, the default), "\code{greater}" or
#' "\code{less}".
#' @param CI a logical variable indicating whether approximate 95% and
#' 99% randomization confidence intervals will be calculated (\code{TRUE})
#' or not (\code{FALSE}).
#' @param silent a logical variable indicating whether calculation results are
#' printed to the R console (\code{silent = FALSE}). If \code{TRUE} then
#' calculation results are not printed to the R console (useful for
#' simulations)
#' @param seed a single value, interpreted as an integer, or \code{NULL} (see
#' "Details").
#' @return The function returns a \code{RT} result object (list)
#' @author Jorge Navarro-Alberto
#' @seealso \code{\link[rbmc]{summary.RT}} and the main plotting function
#' \code{\link[rbmc]{plot.RT}}
#' @references Manly, B.F.J. and Navarro-Alberto, J.A. (2021) Randomization,
#' Bootstrap and Monte Carlo Methods in Biology. 4th Edition. Chapman and
#' Hall/ CRC Press.
#' @examples
#'
#' # Example in Manly and Navarro Alberto (2021), Section 1.1
#' male <- jackals$Mand.length[jackals$Sex=="M"]
#' female <- jackals$Mand.length[jackals$Sex=="F"]
#' jackals.RT2 <- RT2SAMP(male, female, alt="greater")
#' @export RT2SAMP
RT2SAMP <- function(x, y, NRAND = 4999, alt = "two", CI = FALSE,
                    silent = FALSE, seed = NULL)
# This function computes a permutation test of comparison of the means
# of two vectors corresponding to independent samples.  The
# permutational two-sample randomization test requires equality or
# near-equality of the within-group variances. It does not require
# normality of the distributions.
# Parameters of the function.
# x, y: the two vectors to be compared
# NRAND = number of permutations (default value: 4999)
# alt = "greater": one-tailed test in the upper tail
# alt = "lower": one-tailed test in the lower tail
# alt = "two": two-tailed test (default value)
# CI = FALSE: Does not calculate an approximate 95% and 99% confidence interval
# CI = TRUE: Calculates an approximate 95% and 99% CI's by interpolation
# silent = FALSE: calculation results are printed to the R console.
# silent = TRUE: calculation results are not printed to
# the R console (useful for simulations).  Values returned:
# Adapted from a source file function by Pierre Legendre, November 2005
{
  rt2 <- list()
  D.stat <- function(n1, n2, x, y) # Compute the D-statistic
  {
    mean1 <- mean(x)
    mean2 <- mean(y)
    D <- mean1 - mean2
    return(list(mean1 = mean1, mean2 = mean2, stat = D))
  }
  # Find the lengths of the two vectors
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  # If CI=TRUE, approximate CIs are calculated
  D <- NULL
  # Compute the D-statistic for the observed data
  D.obs <- D.stat(n1, n2, x, y)
  D[1] <- D.obs$stat  # Start the vector containing the list of D statistics

  # Perform the permutation test
  nTLG <- rep(1, 3)
  P <- rep(NA, 3)
  for (index in 2:(NRAND + 1)) {
    samp <- sample(n, n1)
    vec <- c(x, y)
    vec1.perm <- vec[samp]
    vec2.perm <- vec[-samp]
    D.perm <- D.stat(n1, n2, vec1.perm, vec2.perm)
    D[index] <- D.perm$stat
    EPS <- 1e-05 * D[1]
    if (D[index] <= D[1] + EPS)
      nTLG[1] <- nTLG[1] + 1
    if (D[index] >= D[1] - EPS)
      nTLG[2] <- nTLG[2] + 1
    if (abs(D[index]) >= abs(D[1]) - EPS)
      nTLG[3] <- nTLG[3] + 1
  }

  # Compute and print the permutational p-value
  P <- nTLG * 100/(NRAND + 1)
  names(P) <- c("Less", "Greater", "Two")
  type <- ifelse(alt == "less", 1, ifelse(alt == "greater", 2, 3))
  if (!silent)
    cat("Sig. level (%) (", NRAND, "permutations) =", formatC(P[type], digits = 2,
                                                    width = 5, format = "f"), "\n")
  if (!silent) {
    if (alt == "less")
      cat("Lower tailed test\n")
    if (alt == "greater")
      cat("Upper tailed test\n")
    if (alt == "two")
      cat("Two tailed test\n")
  }
  # If CI=TRUE, approximate CIs are calculated
  CI9599 <- NULL
  if (CI) {
    # Compute the D-statistic for the observed data
    D.data <- D.stat(n1, n2, x, y)
    (DMEAN <- D.data$stat)
    (SED <- sqrt(var(x)/n1 + var(y)/n2))
    (DMIN <- DMEAN - 3.5 * SED)
    (DINC <- 0.2 * SED)
    # Try Delta values
    DELTA <- numeric(36)
    Q <- matrix(0, nrow = 36, ncol = 2)
    for (i in 1:36) {
      DELTA[i] <- DMIN + (i - 1) * DINC
      vec1.CI <- x - DELTA[i]
      vec2.CI <- y
      D.obsCI <- D.stat(n1, n2, vec1.CI, vec2.CI)
      D.CI <- numeric(NRAND + 1)
      D.CI[1] <- D.obsCI$stat
      nTLG.CI <- rep(1, 2)
      P.CI <- rep(NA, 2)
      set.seed(seed)
      for (index in 2:(NRAND + 1)) {
        EPS <- 1e-05 * D.CI[1]
        sampCI <- sample(n, n1)
        vecCI <- c(vec1.CI, vec2.CI)
        vec1.CIperm <- vecCI[sampCI]
        vec2.CIperm <- vecCI[-sampCI]
        D.perm.CI <- D.stat(n1, n2, vec1.CIperm, vec2.CIperm)
        D.CI[index] <- D.perm.CI$stat
        if (D.CI[index] <= (D.CI[1] + EPS)) {
          nTLG.CI[1] <- nTLG.CI[1] + 1
        }
        if (D.CI[index] >= (D.CI[1] - EPS)) {
          nTLG.CI[2] <- nTLG.CI[2] + 1
        }
      }
      P.CI <- 100 * nTLG.CI/(NRAND + 1)
      Q[i, 1:2] <- P.CI[1:2]
      print(c(DELTA[i], Q[i, 1], Q[i, 2]))
      DELTA0 <- DELTA[i] - DINC
      if (i == 1) {
        DL25 <- DELTA0
        DL05 <- DL25
        DU25 <- DMEAN + 3.5 * SED
        DU05 <- DU25
      } else {
        if ((Q[i, 2] >= 2.5) && (Q[i - 1, 2] < 2.5))
          DL25 <- approx(c(Q[i - 1, 2], Q[i, 2]), c(DELTA0, DELTA[i]),
                         xout = 2.5)$y else if ((Q[i, 1] < 2.5) && (Q[i - 1, 1] >= 2.5))
                           DU25 <- approx(c(Q[i - 1, 1], Q[i, 1]), c(DELTA0, DELTA[i]),
                                          xout = 2.5)$y else if ((Q[i, 2] >= 0.5) && (Q[i - 1, 2] < 0.5))
                                            DL05 <- approx(c(Q[i, 2], Q[i - 1, 2]), c(DELTA0, DELTA[i]),
                                                           xout = 0.5)$y else if ((Q[i, 1] < 0.5) && (Q[i - 1, 1] >= 0.5))
                                                             DU05 <- approx(c(Q[i - 1, 1], Q[i, 1]), c(DELTA0, DELTA[i]),
                                                                            xout = 0.5)$y
      }
    }
    cat("  Approximate confidence limits by interpolation\n")
    CI9599 <- c(DL25, DU25, DL05, DU05)
    print(CI9599)
  }
  rt2 <- list(Sobs = D[1], name="Two sample RT for mean difference",
              p.perm = P, NRAND = NRAND, alt = alt, CI9599 = CI9599)
  class(rt2) <- "RT"
  return(rt2)
}
