# Functions for creating and manipulating distance matrices

#' Compute all pairwise hamming distances between matrix rows.
#'
#' Returns a distance matrix of all meanings. The meanings argument should be
#' a matrix with the different dimensions encoded along columns, and all
#' combinations of meanings specified along rows. The data type of the cells
#' does not matter since distance is simply based on equality - in fact
#' specifying a meaning component as NA allows you to ignore that dimension for
#' the given row/meaning combinations (see examples).
#'
#' @examples
#' # example #1: a 2x2 design using strings (the character redundantly
#' # specifies the dimension for clarity)
#' hammingdists(matrix(c("a1", "b1", "a1", "b2", "a2", "b1", "a2", "b2"), ncol=2, byrow=TRUE))
#' # example #2: a 2x3 design using integers for encoding
#' hammingdists(matrix(c(0,0,0,1,0,2,1,0,1,1,1,2), ncol=2, byrow=TRUE))
#' # example #3: a 2x2x2 design using factors (ncol is always the number of dimensions)
#' hammingdists(matrix(c(0,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,0,1,1,1), ncol=3, byrow=TRUE))
#' # example #4: if some meaning dimension is not relevant for some
#' # combinations of meanings (e.g. optional arguments), specifying them as NA
#' # in the matrix will make them not be counted towards the hamming distance!
#' # in this example the value of the second dimension does not matter (and
#' # does not count towards the distance) when the the first dimension has
#' # value '1'
#' hammingdists(matrix(c(0,0,0,1,1,NA), ncol=2, byrow=TRUE))
#' @seealso \code{\link{stats::dist}}
#' @export
hammingdists <- function(meanings) {
  nmeanings <- dim(meanings)[1]
  matrix <- matrix(nrow=nmeanings, ncol=nmeanings)
  combs <- combn(1:nmeanings, 2)
  # fill only the lower half to pass to as.dist
  matrix[cbind(combs[2,], combs[1,])] <- apply(combs, 2, function(x) sum(meanings[x[1],] != meanings[x[2],], na.rm=TRUE))
  # TODO if the meaning dimensions are named then unlist(meanings[x[1],])
  # might have to be used internally
  as.dist(matrix)
}

#' Compute the normalised Levenshtein distances between strings.
#' 
#' @param strings a vector or list of strings
#' @return a distance matrix of normalised Levenshtein distances between the strings
#' @examples normalisedlevenshteindists(c("abd", "absolute", "asdasd", "casd"))
#' @seealso \code{\link{stats::dist}}
#' @export
normalisedlevenshteindists <- function(strings) {
  # if you want more than just simple Levenshtein distance, have a look at the
  # stringdist package: http://cran.r-project.org/web/packages/stringdist/index.html
  levs <- adist(strings)
  lens <- sapply(strings, nchar)
  as.dist(levs / outer(lens, lens, pmax))
}

#' Create a symmetric distance matrix.
#' 
#' If \code{m} is a matrix, check whether it is a valid specification of a
#' distance matrix and return it, making it symmetric if it isn't already.
#' For all other obect types, try to coerce \code{m} to a \code{dist} object
#' and return the corresponding distance matrix.
#' 
#' @return a symmetric matrix with 0s in the diagonal
#' @seealso \code{\link{stats::dist}}
check.dist <- function(m) {
  if (is.matrix(m)) {
    if (any(diag(m) != 0)) {
      stop("Not a valid distance matrix: nonzero diagonal entries")
    } else if (any(m[upper.tri(m)] != 0, na.rm = TRUE) && !isSymmetric(m)) {
      if (any(m[lower.tri(m)] != 0, na.rm = TRUE)) {
        stop("Not a valid distance matrix: distances are not symmetric. If only specifying one side then enter distances in the lower triangle of the matrix.")
      } else {
        m <- t(m)
      }
    }
  }
  as.matrix(tryCatch(as.dist(m)), warning=function(w)stop(w))
}

#' Permute the rows/columns of a matrix.
#'
#' Returns the given matrix with rows and columns permuted in the same order.
#'
#' @param dist a matrix with an equal number of rows and columns
#' @param perm vector of indices specifying the new order of rows/columns
#' @return a matrix of the same size as \code{m}
shuffle.locations <- function(m, perm = sample.int(dim(m)[1])) {
  m[perm, perm]
}

#' Aggregate rows and columns of a symmetric (distance) matrix.
#' 
#' Conflates the given matrix by either averaging across or discarding
#' rows which are marked by the same group index/id.
#' 
#' @param m a symmetric matrix, can be a matrix or \code{dist} object
#' @param groups a list of group indices/ids of the same length as the
#'   dimension of the matrix, indicating which rows/columns should be
#'   aggregated over
#' @return a matrix with as many rows and columns as there are unique ids in
#'   the \code{groups} parameter
conflate.rows <- function(m, groups, discard=FALSE) {
  newsamples <- unique(groups)

  if (discard) {
    # drop identical (zero-distance) entries - add a safety check!
    as.matrix(d)[newsamples,newsamples]
  } else { # average
    # discard rows which will be conflated
    d <- as.matrix(d)[newsamples,]
    # calculate average distance to other datapoints
    d <- tapply(d, rep(groups, each=nrow(d)), function(collated)rowMeans(matrix(collated,nrow=nrow(d))))
    # turn back into a matrix
    t(sapply(d,I))
  }
}
