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
