#' Create a vector of 'temperature' (from blue over white to red) colors.
#'
#' @param mn when \code{mx} is not specified: number of colors (>1) in palette
#' @param mx maximum
#' @param intensity saturation of the most extreme color(s) 
#'
#' @examples
#' image(as.matrix(1:7), z=as.matrix(1:7), col=temperature.colors(7))
#' image(as.matrix(1:7), z=as.matrix(1:7), col=temperature.colors(-4, 2))
#' @export
temperature.colors <- function(mn, mx=NULL, intensity=1) {
  if (is.null(mx)) {
    mx <- floor(mn/2)
    mn <- ceiling(-mn/2)
  }
  grDevices::hsv(c(rep(0.65, abs(mn)), FALSE, rep(0, abs(mx))),
    intensity*abs(mn:mx)/max(abs(c(mn, mx))))
}

#' Convert p-values into strings with inequalities.
#'
#' Format p value into a string representation that can be used for publishing.
#'
#' @param ps a vector (or single scalar) of significance levels
#' @param digits round maximum digits displayed
#' @param thresholds vector of increasing significance levels at which exact
#'   p values should be replaced with "<threshold" labels
#'
#' @examples
#' pvalue.str(1:8/101)
#' @export
pvalue.str <- function(ps, digits=3, thresholds=c(.001, .01, .05)) {
  lvls <- apply(outer(as.vector(ps), thresholds, '<'), 1, function(r) match(TRUE, r))
  out <- paste('<', thresholds, sep='')[lvls]
  # could add = too?
  out[is.na(lvls)] <- round(ps[is.na(lvls)], digits)
  out
}
