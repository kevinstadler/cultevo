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
#' repmatrix(diag(4))
#' repmatrix(diag(4), times=2)
#' repmatrix(diag(4), each=2)
#' repmatrix(diag(3), times=2, each=2)
#' repmatrix(diag(4), each.row=2)
#' repmatrix(diag(4), times.row=2)
#' @seealso \code{\link[base]{rep}}
#' @export repmatrix
repmatrix <- function(x, times=1, each=1, times.row=times, times.col=times, each.row=each, each.col=each, ...) {
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
#'   reused across meaning dimensions or not. When \code{uniquelabels = FALSE},
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
    dimnames=list(NULL, rev(colnames(combs))))
}

#' Convert a meaning matrix to a binary 'meaning-feature-present' matrix.
#' 
#' Transforms a meaning matrix to 'wide' format where, instead of having
#' a column for every meaning dimension store all possible meaning values,
#' every possible value for any dimension is treated as its own categorical
#' 'meaning feature' whose presence or absence is represented by a logical
#' \code{TRUE}/\code{FALSE} value in its own meaning feature column.
#'
#' Given a matrix or data frame with meaning dimensions along columns and
#' different combinations of meaning feature values along rows, creates a
#' a matrix with the same number of rows but with one column for every
#' possible value for every meaning dimension.
#'
#' All meaning dimensions and values are treated \emph{categorically}, i.e. as
#' factors with no gradual notion of meaning feature similarity, neither
#' within nor across the original meaning dimensions. Information about which
#' feature values correspond to which meaning dimensions is essentially
#' discarded in this representation, but could in principle be recovered
#' through the patterns of (non)-co-occurrence of different meaning features.
#'
#' In order for the resulting meaning columns to be interpretable, the column
#' names of the result are of the structure \code{columnname=value}, based on
#' the column names of the input meaning matrix (see Examples).
#'
#' @param meanings a matrix or data frame with meaning dimensions along columns
#'   and different meaning combinations along rows (such as created by
#'   \code{\link{enumerate.meaningcombinations}}).
#' @param rownames optional character vector of the same length as the number
#'   of rows of \code{meanings}.
#' @return A matrix of \code{TRUE}/\code{FALSE} values with as many rows as
#'   \code{meanings} and one column for every column-value combination in
#'   \code{meanings}.
#' @examples
#' enumerate.meaningcombinations(c(2, 2))
#' binaryfeaturematrix(enumerate.meaningcombinations(c(2, 2)))
#' @export
binaryfeaturematrix <- function(meanings, rownames=NULL) {
  if (is.logical(meanings)) {
    if (!is.null(row.names))
      row.names(meanings) <- rownames
    meanings
  } else {
    if (is.null(colnames(meanings)))
      mdimnames <- paste("V", seq(ncol(meanings)), sep="")
    else
      mdimnames <- colnames(meanings)

    do.call(cbind, lapply(1:ncol(meanings), function(i) {
      # TODO special case for logical columns?
      features <- unique(meanings[, i])
      cols <- matrix(FALSE, nrow=nrow(meanings), ncol=length(features))
      cols[cbind(1:nrow(cols), match(meanings[, i], features))] <- TRUE
      structure(cols,
        dimnames=list(rownames, paste(mdimnames[i], features, sep="=")))
    }))
  }
}
