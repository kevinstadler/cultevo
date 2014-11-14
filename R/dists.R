# Functions for creating and manipulating distance matrices

#' Compute all pairwise hamming distances between matrix rows.
#'
#' Returns a distance matrix of all meanings. The meanings argument should be
#'
#' @param meanings a matrix with the different dimensions encoded along columns,
#'   and all combinations of meanings specified along rows. The data type of the
#'   cells does not matter since distance is simply based on equality - in fact
#'   specifying a meaning component as NA allows you to ignore that dimension
#'   for the given row/meaning combinations (see examples).
#' @return a distance matrix of type \code{\link{dist}} with \code{n*(n-1)/2}
#'   rows/columns, where n is the number of rows in \code{meanings}.
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
#' @seealso \code{\link{dist}}
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
#' @seealso \code{\link{dist}}
#' @export
normalisedlevenshteindists <- function(strings) {
  # if you want more than just simple Levenshtein distance, have a look at the
  # stringdist package: http://cran.r-project.org/web/packages/stringdist/index.html
  levs <- adist(strings)
  lens <- sapply(strings, nchar)
  as.dist(levs / outer(lens, lens, pmax))
}

#' Construct a distance matrix by reading distances from a file or dataframe
#'
#' @param data a filename, dataframe or matrix
#' @param el1.column the column name or id specifying the first element
#' @param el2.column the column name or id specifying the second element
#' @param dist.columns the column name(s) or id(s) specifying the distance(s)
#'   between the two corresponding elements
#' @return a distance matrix (or list of distance matrixes when there is more
#'   than one \code{dist.columns}) of type \code{matrix}
#' @examples
#' read.dist(cbind(c(1,1,1,2,2,3), c(2,3,4,3,4,4), 1:6, 6:1), dist.columns=c(3,4))
#' @export
read.dist <- function(data, el1.column=1, el2.column=2, dist.columns=3) {
  if (is.character(data)) {
    data <- read.csv(filename, header=TRUE)
  }
  els <- sort(unique(c(data[,el1.column], data[,el2.column])))
  # create as many matrices as there are dist.columns
  d <- replicate(length(dist.columns), matrix(0, nrow=length(els), ncol=length(els)))
  el1s <- match(data[,el1.column], els)
  el2s <- match(data[,el2.column], els)
  indices <- cbind(rep.matrix(cbind(el1s, el2s), each.row=length(dist.columns)), 1:length(dist.columns))
  d[indices] <- as.vector(t(data[,dist.columns]))
  if (length(dist.columns) == 1) {
    return(d[,,1])
  } else {
    # turn 3d array into named list
    d <- lapply(1:length(dist.columns), function(i)d[,,i])
    names(d) <- dist.columns
    return(d)
  }
}

#' Check or fix a distance matrix.
#' 
#' Checks or fixes the given distance matrix specification and, if possible,
#' returns an equivalent symmetric \code{matrix} object with 0s in the diagonal.
#' 
#' If the argument is a matrix, check whether it is a valid specification of a
#' distance matrix and return it, making it symmetric if it isn't already.
#' 
#' If the argument is a list, calls \code{check.dist} on every of its elements
#' and returns a list of the results.
#' 
#' For all other object types, attempts to coerce \code{m} to a \code{dist}
#' object and return the corresponding distance matrix.
#' 
#' @param d an object (or list of objects) specifying a distance matrix
#' @return a symmetric \code{matrix} object (or list of such objects) of the
#'   same dimension as \code{d}
#' @seealso \code{\link{dist}}
#' @export
check.dist <- function(d) {
  UseMethod("check.dist", d)
}

#' @S3method check.dist dist
check.dist.dist <- as.matrix

#' @export
check.dist.default <- function(d) {
  as.matrix(tryCatch(as.dist(d), warning=function(w)stop(w)))
}

#' @export
check.dist.matrix <- function(m) {
  if (any(diag(m) != 0)) {
    stop("Not a valid distance matrix: nonzero diagonal entries")
  } else if (any(m[upper.tri(m)] != 0, na.rm = TRUE) && !isSymmetric(m)) {
    if (any(m[lower.tri(m)] != 0, na.rm = TRUE)) {
      stop("Not a valid distance matrix: distances are not symmetric. The matrix either has to be symmetric or have one of its triangles unspecified (0 or NA).")
    } else {
      m <- t(m)
    }
  }
  check.dist.default(m)
}

check.dist.list <- function(ms) {
  ms <- lapply(ms, check.dist)
  # TODO check if all elements have the same dimension?
  attr(ms, "dim") <- dim(ms[1])
  return(ms)
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
    as.matrix(m)[newsamples,newsamples]
  } else { # average
    # discard rows which will be conflated
    m <- as.matrix(m)[newsamples,]
    # calculate average distance to other datapoints
    m <- tapply(m, rep(groups, each=nrow(m)), function(collated)rowMeans(matrix(collated,nrow=nrow(m))))
    # turn back into a matrix
    t(sapply(m,I))
  }
}
