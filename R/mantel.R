#' An R package for computing distance matrices and performing Mantel tests
#' @name mantel

if (suppressWarnings(require(gputools, quietly=TRUE))) {
  cor <- gpuCor
} else {
  message("If you have an nVidia graphics card and fancy faster calculations, consider installing the R package 'gputools'.")
}


#' Perform a Mantel test.
#' 
#' Performs a Mantel permutation test on two distance matrices.
#' 
#' @param m1 a (meaning) distance matrix
#' @param m2 a (string) distance matrix
#' @param maxtrials maximum number of permutations to be tested
#' @param conflate if \code{TRUE}, entries with 0 distance in \code{m1} will be
#'   conflated by averaging over the corresponding entries in \code{m2}
#' @param shuffle a function applied to \code{m2} to shuffle entries, taking a
#'   matrix and an optional permutation specification as its arguments - you
#'   shouldn't normally have to change this
#' @return a list specifying the results of the Mantel test
#' @examples
#' # small distance matrix, Mantel test run deterministically
#' mantel.test(dist(1:4), dist(1:4))
#' # smallest distance matrix using random permutations (with maxtrials=1000)
#' mantel.test(dist(1:7), dist(1:7))
#' @seealso \code{\link{hammingdists}}, \code{\link{normalisedlevenshteindists}}
#' @export
mantel.test <- function(m1, m2, maxtrials=1000, conflate=FALSE, shuffle=shuffle.locations) {
  m1 <- check.dist(m1)
  m2 <- check.dist(m2)
  d <- dim(m2)[1]
  size <- d/dim(m1)[1]
  if (size != 1) {
    if (size %% 1 == 0) {
      message("Dimensionality of second matrix is a multiple of the first, replicating the first")
      m1 <- rep.matrix(m1, times=size)
    } else {
      stop("The two distance matrices have incompatible dimensions")
    }
  }
  if (conflate) {
    # a row's unique referent is the number of the first column where
    # the matrix contains a zero in that row
    uniqueids <- apply(m1, 1, function(row)match(0,row))
    m1 <- conflate.rows(m1, uniqueids, TRUE)
    m2 <- conflate.rows(m2, uniqueids, FALSE)
  }
  indices <- which(lower.tri(m1))
  # extract values relevent for correlation computation
  m1 <- m1[indices]
  veridical <- cor(m1, m2[indices])
  # determine whether it would be more efficient to enumerate all permutations:
  # assuming we take 'maxtrials' samples, how many of the factorial(d) possible
  # values haven't been selected yet compared to the maximum number of trials?
  enumerate <- tryCatch(maxtrials/((factorial(d)*pgeom(maxtrials, 1/factorial(d), lower.tail=FALSE))) >= 1, warning=function(x)FALSE)
  if (enumerate && length(find.package("combinat", quiet=TRUE))) {
    message("Permutation space is small, enumerating all ", factorial(d), " possible permutations.")
    msample <- sapply(combinat::permn(d), function(perm)cor(m1, shuffle(m2, perm)[indices]))
  } else {
    msample <- replicate(maxtrials, cor(m1, shuffle(m2)[indices]))
  }
  c <- sum(msample >= veridical)
  mn <- mean(msample)
  s <- sd(msample)
  z <- (veridical-mn)/s
  return(list(mean=mn, sd=s, veridical=veridical, p=c/maxtrials, z=z, msample=msample))
}

#' Run Mantel tests on a set of consecutive data.
#' 
#' Performs multiple Mantel tests and returns their results in a matrix,
#' optionally visualising part of the data in a plot. This function can be
#' called with experimental data, distance matrix calculation is taken care
#' of internally.
#' 
#' @param meanings a matrix specifying all meaning combinations, as described in \link{\code{hammingdists}}
#' @param strings either a vector of strings, or a matrix of strings per generation
#' @param test.args a list of named arguments passed on to \code{mantel.test}
#' @param plot specifies which mantel test outcome should be plotted (if any):
#'   currently supported are \code{"r"} and \code{"z"}.
#' @param ... extra arguments passed on to the plotting function - here you
#'   might want to specify parameters like \code{ylim}..
#' @return a matrix of test results
#' @examples
#' mantel.development(allmeaningcombinations(c(2,2)), c("asd", "asdf", "", "f"), plot="r")
#' m <- matrix(c("sadasd", "iuerwh", "sdfgkj", "uofidsgf", "asd", "asdf", "", "f"), nrow=2, byrow=T)
#' mantel.development(allmeaningcombinations(c(2,2)), m, plot="r")
#' mantel.development(allmeaningcombinations(c(2,2)), m, plot="z")
#' @seealso \link{\code{hammingdists}}
#' @seealso \link{\code{mantel.test}}
#' @export
mantel.development <- function(meanings, strings, test.args=NULL, ..., plot=NULL) {
  if (is.vector(strings)) {
    strings <- t(strings)
  }
  d1 <- hammingdists(meanings)
  mantels <- do.call(rbind, apply(strings, 1, function(row)do.call(mantel.test, c(list(m1=d1, m2=normalisedlevenshteindists(row)), test.args))))
  rownames(mantels) <- 0:(nrow(mantels)-1)
  if (!is.null(plot)) {
    plot.mantels(mantels, plot, ...)
  }
  return(mantels)
}

#' Read data from a tab or comma-separated file and run Mantel tests on the data.
#' 
#' Reads the given .csv file and runs Mantel tests on it, optionally plotting
#' the results. Plotting is controlled solely through the \code{...} arguments
#' which are passed on to \link{\code{mantel.development}} and
#' \link{\code{plot.mantels}} respectively.
#' If none of the column-arguments are specified it is assumed that the first
#' column gives the strings and all remaining columns encode the different
#' meanings. See \code{mantelexample.csv} for a minimal example.
#' 
#' @param filename csv file to be read - not specifying a file will result in
#'   a file selection dialog being opened
#' @param header whether the csv file contains a header line or not
#' @param stringcolumns name(s) or index(es) of the columns containing the strings
#' @param meaningcolumns names or indices of the columns specifying the meanings
#' @param generationcolumn if the strings from consecutive generations are
#'   entered along rows rather than columns, this parameter indicates the column
#'   which specifies the generation to which a row belongs. When this parameter
#'   is specified, \code{stringcolumns} should be a single index/name.
#' @param ... extra arguments are passed on to \link{\code{mantel.development}}
#'   (and potentially further to \link{\code{mantel.test}} and \link{\code{plot.mantels}})
#' @examples
#' mantel.file(system.file("minimalexample.csv", package="mantel"), plot="r")
#' mantel.file(system.file("minimalexample.csv", package="mantel"), plot="r", test.args=list(maxtrials=10, conflate=TRUE))
#' mantel.file(system.file("generationsalongrows.csv", package="mantel"), plot="r", generationcolumn=4)
#' @seealso \link{\code{mantel.development}}
#' @seealso \link{\code{mantel.test}}
#' @seealso \link{\code{plot.mantels}}
#' @export
mantel.file <- function(filename=NULL, sep="\t", header=FALSE, stringcolumns=1, meaningcolumns=-c(stringcolumns,generationcolumn), generationcolumn=NULL, test.args=NULL, ...) {
  if (is.null(filename)) {
    filename <- file.choose()
  }
  data <- read.table(filename, stringsAsFactors=FALSE, sep=sep, header=header)
  if (!is.null(generationcolumn)) {
    data <- aggregate(subset(data[order(data[,generationcolumn]),], select=stringcolumns), by=do.call(list, data[,meaningcolumns]), FUN=c)
    meaningcolumns <- 1:length(meaningcolumns)
    stringcolumns <- length(meaningcolumns)+1
  }
  mantel.development(data[,meaningcolumns], t(data[,stringcolumns]), test.args=test.args, ...)
}
