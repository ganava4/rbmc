#' Sampling a randomization distribution
#'
#' sampleRD samples the randomization distribution for Fisher or paired sample
#' test
#'
#' The procedure calculates first the sumD-statistic for the observed data,
#' where \code{sumD} stands for the differences \code{x - mu}, and \code{mu} is
#' a hypothesized value for mu under a null hypothesis (usually set to 0, the
#' default. Then, the randomization procedure decides to change sign of each
#' \code{x - mu} or not, after a coin flip. The number of coin flips are
#' controlled by the parameter \code{NRAND} thus producing the randomization
#' distribution of \code{sumD} differences.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param mu an scalar; it is the value of \code{mu} under the null hypothesis
#' @param NRAND numeric; the number of randomizations (permutations).
#'
#' @return The function returns a list containing the value of the \code{sumD}
#' statistic and the p-values of the sampled randomization distribution for the
#' three types of alternatives, and the number of randomizations + 1
#' @author Jorge Navarro-Alberto
#' @seealso \code{RT1SAMP}
#' @references Manly, B.F.J. and Navarro-Alberto, J.A. (2021) Randomization,
#' Bootstrap and Monte Carlo Methods in Biology. 4th Edition. Chapman and
#' Hall/ CRC Press.
#' @examples
#' data(cornheight)
#' Dif <- cornheight[,"Cross"] - cornheight[,"Self"]
#' Results.RD <- sampleRD(Dif)
#'
#' @export
sampleRD <- function(x, mu=0, NRAND=4999) {
  # Fisher's (paired-sample) RT, by sampling the randomization distribution
  # Find the length of the observed vector
  n <- length(x)
  # Compute the sumD-statistic for the observed data
  sumD0 <- sum(x - mu)
  abs.sumD0 <- abs(sumD0)
  nTLG <- rep(1, 3)
  P <- vector(mode="numeric", length=3)
  EPS <- 1e-05
  for (j in 2:NRAND) {
  # Decide to change sign or not, after a coin flip
  for (i in 1:n) {
    if (runif(1) < 0.5) x[i] <- -x[i]
  }
  # Compute the sumD-statistic
  sumD <- sum(x - mu)
  if (sumD < sumD0 + EPS)
    nTLG[1] <- nTLG[1] + 1
  if (sumD > sumD0 - EPS)
    nTLG[2] <- nTLG[2] + 1
  if (abs(sumD) > abs.sumD0 - EPS)
    nTLG[3] <- nTLG[3] + 1
}
  # Compute the p-value of the sampled randomization distribution
  P <- nTLG * 100/NRAND
  names(P) <- c("Less", "Greater", "Two")
  sample.stats <- list(Sobs=sumD0,
                      name="Fisher (or paired sample) randomization test",
                      P = P, nperm = NRAND+1)
  return(sample.stats)
}
