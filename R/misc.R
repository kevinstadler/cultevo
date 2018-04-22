#' Create a vector of 'temperature' colors (from blue over white to red).
#'
#' @param mn integer: when \code{mx} is not specified, total number of colors
#'   (>1) in the palette. when \code{mx} is specified: 'coldest' temperature (see examples)
#' @param mx integer: 'warmest' temperature (see examples)
#' @param intensity saturation of the most extreme color(s), in the range `[0,1]`.
#'
#' @examples
#' # full intensity
#' image(as.matrix(1:7), z=as.matrix(1:7), col=temperature.colors(7))
#' # half intensity
#' image(as.matrix(1:7), z=as.matrix(1:7), col=temperature.colors(7, intensity=0.5))
#' # skewed palette with more negative than positive temperature colors
#' image(as.matrix(1:7), z=as.matrix(1:7), col=temperature.colors(-4, 2))
#' @seealso \code{\link[grDevices]{gray}}, \code{\link[grDevices]{hsv}}, \code{\link[grDevices]{rainbow}}
#' @export
temperature.colors <- function(mn, mx=NULL, intensity=1.0) {
  if (is.null(mx)) {
    mx <- floor(mn/2)
    mn <- ceiling(-mn/2)
  }
  grDevices::hsv(c(rep(0.65, abs(mn)), FALSE, rep(0, abs(mx))),
    intensity*abs(mn:mx)/max(abs(c(mn, mx))))
}

# Convert p-values into strings with inequalities.
#
# Format p value into a string representation that can be used for publishing.
#
# @param ps a vector (or single scalar) of significance levels
# @param digits round maximum digits displayed
# @param thresholds vector of increasing significance levels at which exact
#   p values should be replaced with "<threshold" labels
#
# @examples
# pvalue.str(1:8/101)
pvalue.str <- function(ps, digits=3, thresholds=c(.001, .01)) {
  lvls <- apply(outer(as.vector(ps), thresholds, '<'), 1, function(r) match(TRUE, r))
  out <- paste('<', thresholds, sep='')[lvls]
  # TODO add option for prepending "=" to values greater than all thresholds
  out[is.na(lvls)] <- round(ps[is.na(lvls)], digits)
  out
}