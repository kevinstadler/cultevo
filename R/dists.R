#' Make a meaning distance function vectorisable.
#'
#' This function takes as its only argument a function \code{f(m1, m2)} which
#' returns a single numeric indicating the distance between two 'meanings'
#' \code{m1, m2} (which are themselves most likely vectors or lists). Based
#' on \code{f}, this function returns a function \code{g(mm)} which takes as
#' its only argument a matrix or data frame \code{mm} with the meaning
#' elements (equivalent to the ones in \code{m1, m2}) along columns and
#' different meaning combinations (like \code{m1, m2, ...}) along rows. This
#' function returns a distance matrix of class \code{\link[stats]{dist}}
#' containing all pairwise distances between the rows of \code{mm}. The
#' resulting function \code{g} can be passed to other functions in this
#' package, in particular \code{\link{mantel.test}}.
#'
#' The meaning distance function should be commutative, i.e.
#' \code{f(a,b) = f(b,a)}, and meanings should have a distance of zero to
#' themselves, i.e. \code{f(a,a) = 0}.
#' @param pairwisemeaningdistfun a function of two arguments returning a
#'   single numeric indicating the semantic distance between its arguments
#' @return A function that takes a meaning matrix and returns a corresponding
#'   distance matrix of class \code{\link[stats]{dist}}.
#' @examples
#' trivialdistance <- function(a, b) return(a - b)
#'
#' trivialmeanings <- as.matrix(3:1)
#' trivialdistance(trivialmeanings[1], trivialmeanings[2])
#' trivialdistance(trivialmeanings[1], trivialmeanings[3])
#' trivialdistance(trivialmeanings[2], trivialmeanings[3])
#'
#' distmatrixfunction <- wrap.meaningdistfunction(trivialdistance)
#' distmatrixfunction(trivialmeanings)
#' @export
wrap.meaningdistfunction <- function(pairwisemeaningdistfun) {
  f <- function(meanings) {
    if (is.vector(meanings))
      meanings <- as.matrix(meanings)
    nmeanings <- dim(meanings)[1]
    stats::as.dist(cbind(sapply(1:(nmeanings-1), function(i) {
      c(rep(0, i),
        sapply((i+1):nmeanings, function(j)
          pairwisemeaningdistfun(meanings[i,], meanings[j,])))
    }), 0))
#    matrix <- matrix(nrow=nmeanings, ncol=nmeanings)
#    combs <- utils::combn(1:nmeanings, 2)
    # fill only the lower half to pass to as.dist
#    matrix[cbind(combs[2,], combs[1,])] <- apply(combs, 2,
#      function(i) pairwisemeaningdistfun(meanings[i[1],], meanings[i[2],]))
  }
  # cache distance matrix computations
  if (requireNamespace("memoise", quietly = TRUE))
    f <- memoise::memoise(f)
  invisible(f)
}

#' Pairwise Hamming distances between matrix rows.
#'
#' Returns a distance matrix giving all pairwise Hamming distances between the
#' rows of its argument \code{meanings}, which can be a matrix, data frame or
#' vector. Vectors are treated as matrices with a single column, so the
#' distances in its return value can only be 0 or 1.
#'
#' This function behaves differently from calling
#' \code{\link[stats]{dist}(meanings, method="manhattan")} in how \code{NA}
#' values are treated: specifying a meaning component as \code{NA} allows you
#' to \emph{ignore} that dimension for the given row/meaning combinations,
#' (instead of counting a difference between \code{NA} and another value as a
#' distance of 1).
#'
#' @param meanings a matrix with the different dimensions encoded along
#'   columns, and all combinations of meanings specified along rows. The data
#'   type of the cells does not matter since distance is simply based on
#'   equality (with the exception of \code{NA} values, see below.
#' @return A distance matrix of type \code{\link{dist}} with \code{n*(n-1)/2}
#'   rows/columns, where n is the number of rows in \code{meanings}.
#'
#' @examples
#' # a 2x2 design using strings
#' print(strings <- matrix(c("a1", "b1", "a1", "b2", "a2", "b1", "a2", "b2"),
#'   ncol=2, byrow=TRUE))
#' hammingdists(strings)
#'
#' # a 2x3 design using integers
#' print(integers <- matrix(c(0, 0, 0, 1, 0, 2, 1, 0, 1, 1, 1, 2), ncol=2, byrow=TRUE))
#' hammingdists(integers)
#'
#' # a 3x2 design using factors (ncol is always the number of dimensions)
#' print(factors <- data.frame(colour=c("red", "red", "green", "blue"),
#'                             animal=c("dog", "cat", "dog", "cat")))
#' hammingdists(factors)
#'
#' # if some meaning dimension is not relevant for some combinations of
#' # meanings (e.g. optional arguments), specifying them as NA in the matrix
#' # will make them not be counted towards the hamming distance! in this
#' # example the value of the second dimension does not matter (and does not
#' # count towards the distance) when the the first dimension has value '1'
#' print(ignoredimension <- matrix(c(0, 0, 0, 1, 1, NA), ncol=2, byrow=TRUE))
#' hammingdists(ignoredimension)
#'
#' # trivial case of a vector: first and last two elements are identical,
#' # otherwise a difference of one
#' hammingdists(c(0, 0, 1, 1))
#' @seealso \code{\link[stats]{dist}}
#' @export
hammingdists <-
  wrap.meaningdistfunction(function(m1, m2) sum(m1 != m2, na.rm=TRUE))

#' Compute the normalised Levenshtein distances between strings.
#' 
#' @param strings a vector or list of strings
#' @return A distance matrix specifying all pairwise normalised Levenshtein
#'   distances between the strings.
#' @examples normalisedlevenshteindists(c("abd", "absolute", "asdasd", "casd"))
#' @seealso \code{\link[stats]{dist}}
#' @export
normalisedlevenshteindists <- function(strings) {
  levs <- utils::adist(strings)
  lens <- sapply(strings, nchar)
  stats::as.dist(levs / outer(lens, lens, pmax))
}

#' Split strings into their constituent segments.
#'
#' Split strings into their constituent segments (and count them).
#'
#' @describeIn segment.string
#' Returns a list (of the same length as \code{x}), each item a vector of
#'   character vectors.
#' @param x one or more strings to be split (and, optionally, counted)
#' @param split the boundary character or sequence at which to segment the
#'   string(s). The default, \code{NULL}, splits the string after every character.
#' @examples
#' segment.string(c("asd", "fghj"))
#'
#' segment.string(c("la-dee-da", "lala-la"), "-")
#' @export
segment.string <- function(x, split=NULL)
  strsplit(x, split, fixed=TRUE)

#' @describeIn segment.string
#' Calculate the frequency of individual characters in one or more strings.
#' Returns a matrix with one row for every string in \code{x}.
#' @examples
#' segment.counts(c("asd", "aasd", "asdf"))
#' @export
segment.counts <- function(x, split=NULL) {
  chars <- segment.string(x, split)
  set <- unique(unlist(chars))
  structure(t(sapply(chars, function(s) {
    tabulate(match(s, set), nbins=length(set))
  })), dimnames=list(x, set))
}

#' Calculate the bag-of-characters similarity between strings.
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
  indices <- cbind(repmatrix(cbind(el1s, el2s), each.row=length(dist.columns)), 1:length(dist.columns))
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
#' Checks or fixes the given distance matrix specification and returns an
#' equivalent, symmetric \code{matrix} object with 0s in the diagonal.
#' 
#' If the argument is a matrix, check whether it is a valid specification of a
#' distance matrix and return it, making it symmetric if it isn't already.
#' 
#' If the argument is a list, calls \code{check.dist} on every of its elements
#' and returns a list of the results.
#' 
#' For all other object types, attempts to coerce the argument to a \code{dist}
#' object and return the corresponding distance matrix (see above).
#' 
#' @param x an object (or list of objects) specifying a distance matrix
#' @return a symmetric \code{matrix} object (or list of such objects) of the
#'   same dimension as \code{x}
#' @seealso \code{\link[stats]{dist}}
#' @export
check.dist <- function(x)
  UseMethod("check.dist")

#' @export
check.dist.dist <- function(x)
  as.matrix(x)

#' @export
check.dist.default <- function(x)
  as.matrix(tryCatch(stats::as.dist(x), warning=function(w) stop(w)))

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

#' @export
check.dist.list <- function(x)
  lapply(x, check.dist)

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
# https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/pMatrix-class.html
