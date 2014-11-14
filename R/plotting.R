# functions for visualising the outcome of consecutive Mantel tests

# for plotting error bars
if (suppressWarnings(require(Hmisc, quietly=TRUE))) {
  errbar <- Hmisc::errbar
} else {
  message("For prettier error bars, invoke: install.packages('Hmisc')")
  # ingenious solution thanks to http://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
  errbar <- function(x, y, yplus, yminus, errbar.col="black", ...) {
    plot(x, y, ...)
    arrows(x, yplus, x, yminus, length=0.05, angle=90, code=3, col=errbar.col)
  }
}

red.if.na <- function(x) {
  c("black", "red")[1+sapply(x, is.na)]
}

# lookup table for plotting functions for different mantel test measures
mantel.plotfuns <- list()

mantel.plotfuns[["z"]] <- function(mantels, ylim=range(0, mantels[,"z"])+c(0,1), ...) {
  # TODO add msample histogram when nrow(mantels) == 1
  plot(1:nrow(mantels), mantels[,"z"], ylab="z score", ylim=ylim, col=red.if.na(mantels[,"p.smoothed"]), ...)
  abline(h=0, lty=2, col="grey")
}

mantel.plotfuns[["r"]] <- function(mantels, ylim=range(0.1, mantels[,"veridical"], unlist(mantels[,"mean"])+unlist(mantels[,"sd"]), unlist(mantels[,"mean"])-unlist(mantels[,"sd"]))+c(-0.1,0.1), ...) {
  # TODO add msample histogram when nrow(mantels) == 1
  errbar(1:nrow(mantels), mantels[,"mean"], unlist(mantels[,"mean"])+unlist(mantels[,"sd"]), unlist(mantels[,"mean"])-unlist(mantels[,"sd"]), ylab="r", ylim=ylim, errbar.col=red.if.na(mantels[,"p.smoothed"]), ...)
  # plot correlation coefficient reference points
  abline(h=-1:1, lty=c(3,2,3), col="grey")
  # blue points signify that no single larger r value has been sampled
  points(1:nrow(mantels), mantels[,"veridical"], col=c("black", "blue")[1+mantels[,"is.unique.max"]])
  # label the veridical rs with their z scores
  text(1:nrow(mantels), mantels[,"veridical"], labels=paste("z", round(unlist(mantels[,"z"]), digits=2), sep="="), pos=2+sign(unlist(mantels[,"z"])))
}

#' Plot a Mantel test result.
#' 
#' Plots the result of one or more Mantel permutation tests.
#' 
#' @param mantels a Mantel test result (a data frame of class 'mantel')
#' @param plot the measure to be visualised, currently supported: 'r' and 'z'
#' @param ... additional parameters to be passed on to the plotting function
#' @examples
#' plot(mantel.test(hammingdists(allmeaningcombinations(c(2,2,2,2))), suppressWarnings(as.dist(1:16))))
#' @seealso \link{\code{mantel.development}}, \link{\code{plot.default}}
#' @export
plot.mantel <- function(mantels, plot="r", xlab="generation", ...) {
  # suppress x axis drawing in the subroutine
  mantel.plotfuns[[plot]](mantels, xaxt="n", xlab=xlab, ...)
  axis(1, at=1:nrow(mantels), labels=rownames(mantels))
}
