# Functions for creating and manipulating distance matrices

#' Compute all pairwise hamming distances between matrix rows.
#'
#' Returns a distance matrix of all meanings. The meanings argument should be
#'
#' @param meanings a matrix with the different dimensions encoded along
#'   columns, and all combinations of meanings specified along rows. The data
#'   type of the cells does not matter since distance is simply based on
#'   equality - in fact specifying a meaning component as NA allows you to
#'   ignore that dimension for the given row/meaning combinations (see
#'   examples). Vectors are treated as matrices with a single row, so the
#'   resulting distances can only be 0 or 1.
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
#' @seealso \code{\link[stats]{dist}}
#' @export
hammingdists <- function(meanings) {
  if (is.vector(meanings))
    meanings <- as.matrix(meanings)
  nmeanings <- dim(meanings)[1]
  matrix <- matrix(nrow=nmeanings, ncol=nmeanings)
  combs <- utils::combn(1:nmeanings, 2)
  # fill only the lower half to pass to as.dist
  matrix[cbind(combs[2,], combs[1,])] <- apply(combs, 2, function(x) sum(meanings[x[1],] != meanings[x[2],], na.rm=TRUE))
  # TODO if the meaning dimensions are named then unlist(meanings[x[1],])
  # might have to be used internally
  stats::as.dist(matrix)
}

#' Compute the normalised Levenshtein distances between strings.
#' 
#' @param strings a vector or list of strings
#' @return a distance matrix of normalised Levenshtein distances between the strings
#' @examples normalisedlevenshteindists(c("abd", "absolute", "asdasd", "casd"))
#' @seealso \code{\link[stats]{dist}}
#' @export
normalisedlevenshteindists <- function(strings) {
  levs <- utils::adist(strings)
  lens <- sapply(strings, nchar)
  stats::as.dist(levs / outer(lens, lens, pmax))
}

#' Split one or more strings into their constituent segments.
#'
#' Vectorisable.
#' @param x one or more strings to be split
#' @param split the boundary character or sequence at which to segment the strings
#' @return a list (of the same length as \code{x}) containing character vectors
#' @export
segment.string <- function(x, split=NULL)
  strsplit(x, split, fixed=TRUE)

#' Calculate the frequency of individual characters in one or more strings.
#'
#' @param x one or more strings for which segments should be counted
#' @param split the boundary character or sequence at which to segment the strings
#' @return a matrix with one row for every string in \code{x}
#' @export
segment.counts <- function(x, split=NULL) {
  chars <- segment.string(x, split)
  set <- unique(unlist(chars))
  structure(t(sapply(chars, function(s) {
    tabulate(match(s, set), nbins=length(set))
  })), dimnames=list(x, set))
}

#' Calculate the bag-of-characters similarity between the given strings.
#'
#' @param strings a vector or list of strings
#' @param split boundary sequency at which to segment the strings (default
#'   splits the string into all its constituent characters)
#' @param segmentcounts if custom segmentation is required, the pre-segmented
#'   strings can be passed as this argument (which is a list of lists)
#' @return a distance matrix
#' @examples orderinsensitivedists(c("xxxx", "asdf", "asd", "dsa"))
#' @seealso \code{\link[stats]{dist}}
#' @export
orderinsensitivedists <- function(strings=NULL, split=NULL, segmentcounts=segment.counts(strings, split))
  stats::as.dist(cbind(sapply(1:(length(strings)-1), function(i) {
    c(rep(0, i),
      colSums(abs(t(segmentcounts[(i+1):length(strings),,drop=FALSE]) - segmentcounts[i,])))
  }), 0))

#' Read a distance matrix from a file or data frame.
#'
#' @param data a filename, data frame or matrix
#' @param el1.column the column name or id specifying the first element
#' @param el2.column the column name or id specifying the second element
#' @param dist.columns the column name(s) or id(s) specifying the distance(s)
#'   between the two corresponding elements
#' @return a distance matrix (or list of distance matrixes when there is more
#'   than one \code{dist.columns}) of type \code{matrix}
#' @examples
#' read.dist(cbind(c(1,1,1,2,2,3), c(2,3,4,3,4,4), 1:6, 6:1), dist.columns=c(3,4))
#' @export
#' @importFrom utils read.csv
read.dist <- function(data, el1.column=1, el2.column=2, dist.columns=3) {
  if (is.character(data))
    data <- read.csv(data, header=TRUE)

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
#' @param x an object (or list of objects) specifying a distance matrix
#' @return a symmetric \code{matrix} object (or list of such objects) of the
#'   same dimension as \code{d}
#' @seealso \code{\link[stats]{dist}}
#' @export
check.dist <- function(x)
  UseMethod("check.dist")

#' @rdname check.dist
#' @export
check.dist.dist <- function(x)
  as.matrix(x)

#' @rdname check.dist
#' @export
check.dist.default <- function(x)
  as.matrix(tryCatch(stats::as.dist(x), warning=function(w) stop(w)))

#' @rdname check.dist
#' @export
check.dist.matrix <- function(x) {
  if (any(diag(x) != 0)) {
    stop("Not a valid distance matrix: nonzero diagonal entries")
  } else if (any(x[upper.tri(x)] != 0, na.rm = TRUE) && !isSymmetric(x)) {
    if (any(x[lower.tri(x)] != 0, na.rm = TRUE)) {
      stop("Not a valid distance matrix: distances are not symmetric. The matrix either has to be symmetric or have one of its triangles unspecified (0 or NA).")
    } else {
      x <- t(x)
    }
  }
  check.dist.default(x)
}

#' @rdname check.dist
#' @export
check.dist.list <- function(x) {
  x <- lapply(x, check.dist)
  # TODO check if all elements have the same dimension?
  attr(x, "dim") <- dim(x[1])
  return(x)
}

#' Permute the rows and columns of a square matrix.
#'
#' Returns the given matrix with rows and columns permuted in the same order.
#'
#' @param m a matrix with an equal number of rows and columns
#' @param perm vector of indices specifying the new order of rows/columns
#' @return a matrix of the same size as \code{m}
#' @export
shuffle.locations <- function(m, perm = sample.int(dim(m)[1]))
  m[perm, perm]
