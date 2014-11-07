# functions for visualising the outcome of consecutive Mantel tests

# for plotting error bars
if (!suppressWarnings(require(Hmisc, quietly=TRUE))) {
  message("For prettier error bars, invoke: install.packages('Hmisc')")
  # ingenious solution thanks to http://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
  errbar <- function(x, y, yplus, yminus, ...) {
    plot(x, y, ...)
    arrows(x, yplus, x, yminus, length=0.05, angle=90, code=3)
  }
}

# lookup table for plotting functions for different mantel test measures
mantel.plotfuns <- list()

mantel.plotfuns[["z"]] <- function(mantels, ylim=range(0, mantels[,"z"])+c(0,1), ...) {
  plot(1:nrow(mantels), mantels[,"z"], ylab="z score", ylim=ylim, ...)
  abline(h=0, lty=2, col="grey")
}

mantel.plotfuns[["r"]] <- function(mantels, ylim=range(0.1, mantels[,"veridical"], unlist(mantels[,"mean"])+unlist(mantels[,"sd"]), unlist(mantels[,"mean"])-unlist(mantels[,"sd"]))+c(-0.1,0.1), ...) {
  errbar(1:nrow(mantels), mantels[,"mean"], unlist(mantels[,"mean"])+unlist(mantels[,"sd"]), unlist(mantels[,"mean"])-unlist(mantels[,"sd"]), ylab="r", ylim=ylim, ...)
  # plot correlation coefficient reference points
  abline(h=-1:1, lty=c(3,2,3), col="grey")
  points(1:nrow(mantels), mantels[,"veridical"], col="red")
  # label the veridical rs with their z scores
  text(1:nrow(mantels), mantels[,"veridical"], labels=paste("z", round(unlist(mantels[,"z"]), digits=2), sep="="), pos=2+sign(unlist(mantels[,"z"])))
}

#' Plot outcome measures of consecutive Mantel tests.
#' 
#' @param mantels a matrix or data frame of Mantel test results
#' @param plot the measure to be visualised, currently supported: 'r' and 'z'
#' @param ... additional parameters to be passed on to the plotting function
#' @seealso \link{\code{mantel.development}}
plot.mantels <- function(mantels, plot, xlab="generation", ...) {
  # suppress x axis drawing in the subroutine
  mantel.plotfuns[[plot]](mantels, xaxt="n", xlab=xlab, ...)
  axis(1, at=1:nrow(mantels), labels=rownames(mantels))
}
