plot.RT <- function(x, col = "grey", xlab = expression(paste(italic(D))),
                  ylab = "Frequency",
                  main = "Randomization distribution\nTwo sample RT", ...) {
  type.tst <- ifelse(x$alt == "less", "Lower tailed test",
                     ifelse(x$alt == "greater", "Upper tailed test",
                            "Two tailed test"))
  hist(x$D, col = col, xlab = xlab, ylab = ylab, main = main)
  abline(v = x$D[1], col = "red")
  if (x$alt == "two") {abline(v=-x$D[1], col="red")}
  Dform <- format(round(x$D, digits = 2), nsmall = 1)
  mtext(bquote({italic(D)}[0] == .(Dform)), side = 1, at = x$D[1], line = 0.5,
        cex = 0.8)
  if (x$alt == "two") {
    Dform <- format(round(-x$D, digits = 2), nsmall = 1)
    mtext(bquote({italic(-D)}[0] == .(Dform)), side = 1, at = -x$D[1], line = 0.5,
        cex = 0.8)
}
  p.value <- ifelse(x$alt == "less", x$p.perm[1], ifelse(x$alt == "greater",
                                                         x$p.perm[2],
                                                         x$p.perm[3]))
  pform <- format(round(p.value/100, digits = 4), nsmall = 4)
  mtext(bquote({italic(P)} == .(pform)), side = 1, at = x$D[1], line = 1.5,
        cex = 0.8)
  npform <- as.character(x$NRAND)
  mtext(paste(npform, "randomizations"), side = 1, at = x$D[1], line = 2.5,
        cex = 0.7)
  mtext(type.tst, side = 1, at = x$D[1], line = 3.5, cex = 0.7)
}
