# General functions for generating and manipulating matrices

#' Repeat the entries in the matrix along both dimensions.
#' 
#' The resulting matrix will have \code{times*each} times more rows and columns
#' than the original matrix.
#' 
#' @param m a matrix
#' @param times how often the matrix should be replicated next to itself
#' @param each how often individual cells should be replicated next to themselves
#' @param times.row number of vertical repetitions of the matrix, overrides \code{times}
#' @param times.col number of horizontal repetitions of the matrix, overrides \code{times}
#' @param each.row number of vertical repetitions of individual elements, overrides \code{each}
#' @param each.col number of horizontal repetitions of individual elements, overrides \code{each}
#' @return a matrix
#' @examples
#' rep.matrix(diag(4))
#' rep.matrix(diag(4), times=2)
#' rep.matrix(diag(4), each=2)
#' rep.matrix(diag(3), times=2, each=2)
#' rep.matrix(diag(4), each.row=2)
#' rep.matrix(diag(4), times.row=2)
#' @seealso \code{\link{rep}}
#' @export
rep.matrix <- function(m, times=1, each=1, times.row=times, times.col=times, each.row=each, each.col=each) {
  # replicate individual elements
  m <- matrix(rep(m, each=each.row), ncol=each.row*dim(m)[1], byrow=TRUE)
  m <- matrix(rep(m, each=each.col), ncol=each.col*dim(m)[1], byrow=TRUE)
  # replicate entire matrix
  do.call(cbind, replicate(times.col, do.call(rbind, replicate(times.row, m, simplify=FALSE)), simplify=FALSE))
}

#' Enumerate meaning combinations.
#'
#' Enumerates all possible combinations of meanings for a meaning space of the given dimensionality.
#'
#' The resulting matrix can be passed straight on to \code{\link{hammingdists}}.
#'
#' @param dimensionality either a) a vector of integers specifying the number
#'   of different possible values for every meaning dimension, or b) a list or
#'   other (potentially ragged) 2-dimensional data structure listing the
#'   possible meaning values for every dimension
#' @param uniquelabels logical, determines whether the same integers can be
#'   reused across meaning dimensions or not. When \code{uniquelabels==FALSE},
#'   the resulting matrix will be very reminiscent of tables listing all binary
#'   combinations of factors. Ignored when \code{dimensionality} specifies the
#'   meaning values
#' @param offset a constant that is added to all meaning specifiers. Ignored
#'   when \code{dimensionality} specifies the meaning values
#' @return a matrix that has as many columns as there are dimensions, with every
#'   row specifying one of the possible meaning combinations, and the first
#'   dimension cycling slowest.
#' @examples
#' allmeaningcombinations(c(2,2))
#' allmeaningcombinations(c(3,4))
#' allmeaningcombinations(c(2,2,2,2))
#' allmeaningcombinations(8) # trivial
#' allmeaningcombinations(list(shape=c("square", "circle"), color=c("red", "blue")))
#' @seealso \code{\link{hammingdists}}
#' @export
allmeaningcombinations <- function(dimensionality, uniquelabels=TRUE, offset=0) {
  if (is.list(dimensionality)) {
    # argument is a list or other 2-dimensional data structure
    # specifying the possible meanings per dimension
    combs <- expand.grid(rev(dimensionality), stringsAsFactors=FALSE)
  } else {
    # generate cell values
    offsets <- if (uniquelabels) offset+c(rev(cumsum(dimensionality))[-1], 0) else rep(offset, length(dimensionality))
    combs <- do.call(expand.grid, mapply(function(nvalues, off)if (is.na(nvalues)) NA else off+1:nvalues, rev(dimensionality), offsets, SIMPLIFY=FALSE))
  }
  as.matrix(combs[,ncol(combs):1])
}
