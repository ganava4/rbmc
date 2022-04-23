print.RT <- function(x) {
  cl <- oldClass(x)
  oldClass(x) <- cl[cl != "RT"]
  NextMethod("print")
  invisible(x)
 }
