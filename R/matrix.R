#' Extend a matrix by repetition of elements.
#' 
#' Returns a new matrix, where the entries of the original matrix are
#' repeated along both dimensions.
#' 
#' @param x a matrix
#' @param times how often the matrix should be replicated next to itself
#' @param each how often individual cells should be replicated next to themselves
#' @param times.row number of vertical repetitions of the matrix, overrides \code{times}
#' @param times.col number of horizontal repetitions of the matrix, overrides \code{times}
#' @param each.row number of vertical repetitions of individual elements, overrides \code{each}
#' @param each.col number of horizontal repetitions of individual elements, overrides \code{each}
#' @param ... not used
#' @return A matrix, which will have \code{times*each} times more rows and
#' columns than the original matrix.
#' @examples
#' rep.matrix(diag(4))
#' rep.matrix(diag(4), times=2)
#' rep.matrix(diag(4), each=2)
#' rep.matrix(diag(3), times=2, each=2)
#' rep.matrix(diag(4), each.row=2)
#' rep.matrix(diag(4), times.row=2)
#' @seealso \code{\link[base]{rep}}
#' @export
rep.matrix <- function(x, times=1, each=1, times.row=times, times.col=times, each.row=each, each.col=each, ...) {
  # replicate individual elements
  x <- matrix(rep(x, each=each.row), ncol=each.row*nrow(x), byrow=TRUE)
  x <- matrix(rep(x, each=each.col), ncol=each.col*nrow(x), byrow=TRUE)
  # replicate entire matrix
  do.call(cbind, replicate(times.col, do.call(rbind, replicate(times.row, x, simplify=FALSE)), simplify=FALSE))
}

#' Enumerate meaning combinations.
#'
#' Enumerates all possible combinations of meanings for a meaning space of the
#' given dimensionality.
#'
#' The resulting matrix can be passed straight on to
#' \code{\link{hammingdists}} and other meaning distance functions created by
#' \code{\link{wrap.meaningdistfunction}}.
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
#' @return A matrix that has as many columns as there are dimensions, with
#'   every row specifying one of the possible meaning combinations. The entries
#'   of the first dimension cycle slowest (see examples).
#' @examples
#' enumerate.meaningcombinations(c(2, 2))
#' enumerate.meaningcombinations(c(3, 4))
#' enumerate.meaningcombinations(c(2, 2, 2, 2))
#' enumerate.meaningcombinations(8) # trivial
#' enumerate.meaningcombinations(list(shape=c("square", "circle"), color=c("red", "blue")))
#' @seealso \code{\link{hammingdists}}
#' @export
enumerate.meaningcombinations <- function(dimensionality, uniquelabels=TRUE, offset=0) {
  if (is.list(dimensionality)) {
    # argument is a list or other 2-dimensional data structure
    # specifying the possible meanings per dimension
    combs <- expand.grid(rev(dimensionality), stringsAsFactors=FALSE)
  } else {
    # generate cell values
    offsets <- if (uniquelabels) offset+c(rev(cumsum(dimensionality))[-1], 0) else rep(offset, length(dimensionality))
    combs <- do.call(expand.grid, mapply(function(nvalues, off)if (is.na(nvalues)) NA else off+1:nvalues, rev(dimensionality), offsets, SIMPLIFY=FALSE))
  }
  structure(as.matrix(combs[,ncol(combs):1]),
    dimnames=list(NULL, colnames(combs)))
}

#' Convert a categorical meaning matrix to a binary feature matrix.
#' 
#' @param meanings a matrix with meaning dimensions along columns and different
#'   meaning combinations along rows (such as created by
#'   \code{\link{enumerate.meaningcombinations}})
#' @return A matrix of \code{TRUE}/\code{FALSE} values.
#' @examples
#' enumerate.meaningcombinations(c(2, 2))
#' binaryfeaturematrix(enumerate.meaningcombinations(c(2, 2)))
#' @export
binaryfeaturematrix <- function(meanings) {
  if (is.logical(meanings)) {
    meanings
  } else {
    if (is.null(colnames(meanings)))
      mdimnames <- paste("V", seq(ncol(meanings)), sep="")
    else
      mdimnames <- colnames(meanings)

    do.call(cbind, lapply(1:ncol(meanings), function(i) {
      features <- unique(meanings[, i])
      cols <- matrix(FALSE, nrow=nrow(meanings), ncol=length(features))
      cols[cbind(1:nrow(cols), match(meanings[, i], features))] <- TRUE
      structure(cols,
        dimnames=list(NULL, paste(mdimnames[i], features, sep="=")))
    }))
  }
}
