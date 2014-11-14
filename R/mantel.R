#' An R package for computing distance matrices and performing Mantel tests
#' @name mantel

if (suppressWarnings(require(gputools, quietly=TRUE))) {
  cor <- gpuCor
} else {
  message("If you have an nVidia graphics card and fancy faster calculations, consider installing the R package 'gputools'.")
}

# mantel object constructor
mantel.new <- function(...) {
  res <- data.frame(..., stringsAsFactors=FALSE)
  class(res) <- c("mantel", class(res))
  return(res)
}

#' Perform a Mantel test.
#' 
#' Performs a Mantel permutation test on two distance matrices.
#' 
#' Given two distance matrices (or list of distance matrices) \code{m1} and
#' \code{m2}, performs a Mantel permutation test on the two matrices (or all
#' pairwise combinations of matrices if at least one of the parameters is a
#' list of matrices), and returns the outcome of the test(s) in a data frame.
#' 
#' The distance matrices can either be of type \link{\code{dist}}, plain R
#' matrices or any object that can be interpreted by \link{\code{check.dist}}.
#' 
#' If the number of possible permutations is reasonably close to
#' \code{maxtrials}, a deterministic enumeration of all the permutations will
#' be carried out instead of random sampling.
#' 
#' A warning will be produced if the randomly sampled \code{r} coefficients are
#' not normally distributed according to a Kolmogorov-Smirnoff test, at the
#' level specified by \code{ks.level}). Being normally distributed is a
#' condition for the calculation of a meaningful z score. Failing the test of
#' normality will also be signalled by \code{p.smoothed} being set to \code{NA}.
#'
#' @param m1 a (meaning) distance matrix, or list of matrices
#' @param m2 a (string) distance matrix, or list of matrices
#' @param maxtrials maximum number of permutations to be tested
#' @param conflate if \code{TRUE}, entries with 0 distance in \code{m1} will be
#'   conflated by averaging over the corresponding entries in \code{m2}
#' @param ks.level p-level below which a failure to pass 
#' @param shuffle a function applied to \code{m2} to shuffle entries, taking a
#'   matrix and an optional permutation specification as its arguments - you
#'   shouldn't normally have to change this
#' @return a dataframe specifying the results of the Mantel test(s)
#' @examples
#' # small distance matrix, Mantel test run deterministically
#' mantel.test(dist(1:4), dist(1:4))
#' # smallest distance matrix using random permutations (with maxtrials=1000)
#' mantel.test(dist(1:7), dist(1:7))
#' @seealso \code{|link{cor}}, \code{\link{hammingdists}},
#'   \code{\link{normalisedlevenshteindists}}
#' @export
mantel.test <- function(m1, m2, maxtrials=1000, conflate=FALSE, ks.level=0.05, shuffle=shuffle.locations) {
  m1 <- check.dist(m1)
  m2 <- check.dist(m2)
  d <- dim(m2)[1]
  size <- d/dim(m1)[1]
  if (size != 1) {
    if (size %% 1 == 0) {
      warning("Dimensionality of second matrix is a multiple of the first, replicating the first")
      m1 <- rep.matrix(m1, times=size)
    } else {
      stop("The two arguments have incompatible dimensions")
    }
  }
  if (conflate) {
    # a row's unique referent is the number of the first column where
    # the matrix contains a zero in that row
    uniqueids <- apply(m1, 1, function(row)match(0,row))
    if (length(unique(uniqueids)) == d) {
      warning("Specified conflate=TRUE, but there was nothing to conflate.")
    } else {
      m1 <- conflate.rows(m1, uniqueids, TRUE)
      m2 <- conflate.rows(m2, uniqueids, FALSE)
    }
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
  greater <- sum(msample >= veridical)
  p <- (greater+1)/(length(msample)+1)
  mn <- mean(msample)
  s <- sd(msample)
#  names(estimates) <- c("mean", "sd")
  z <- (veridical-mn)/s
#  names(z) <- "z"
  p.smoothed <- pnorm(z, lower.tail=FALSE)
  # suppress warnings about tied r's
  normality <- suppressWarnings(ks.test(msample, "pnorm", mn, s))
  if (normality$p.value < ks.level) {
    warning("The sample of randomised r's is not normally distributed, use that z score with a grain of salt")
    p.smoothed <- NA
  }
  # http://www.inside-r.org/node/219017
  # https://searchcode.com/codesearch/view/13555928/
#  if (length(p) > 1) {
#    combs <- expand.grid(1:length(m1), 1:length(m2))
#    data.name <- paste(combs[,1], "x", combs[,2], sep="")
#  } else {
#    data.name <- "m1, m2"
#  }
  mantel.new(veridical=veridical, mean=mn, sd=s, z=z, p.value=p, p.smoothed=p.smoothed, is.unique.max=(greater==0), sample.size=length(msample), msample=I(list(msample)))
}

#' Run Mantel tests on a set of consecutive data.
#' 
#' Performs multiple Mantel tests and returns their results in a matrix,
#' optionally visualising part of the data in a plot. This function can be
#' called with experimental data, distance matrix calculation is taken care
#' of internally.
#' 
#' @param meanings a matrix specifying all meaning combinations, as described in \link{\code{hammingdists}}
#' @param strings either a vector of strings, or a matrix of strings per
#'   generation, organised along rows. The result of \code{read.table} can be
#'   passed directly (or \code{t(read.table("..."))} if strings from the same
#'   generation are organised in the same column rather than the same row).
#' @param test.args a list of named arguments passed on to \link{\code{mantel.test}}
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
#' mantel.development(allmeaningcombinations(c(2,2)), read.table("stringsonlyexample.csv"), plot="r")
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
    plot.mantel(mantels, plot, ...)
  }
  return(mantels)
}

#' Read data from a tab or comma-separated file and run Mantel tests on the data.
#' 
#' Reads the given .csv file and runs Mantel tests on it, optionally plotting
#' the results. Plotting is controlled solely through the \code{...} arguments
#' which are passed on to \link{\code{mantel.development}} and
#' \link{\code{plot.mantel}} respectively.
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
#' @param test.args a list of named arguments passed on to \link{\code{mantel.test}}
#' @param ... extra arguments are passed on to \link{\code{mantel.development}}
#'   (and potentially further to \link{\code{mantel.test}} and \link{\code{plot.mantel}})
#' @examples
#' mantel.file(system.file("minimalexample.csv", package="mantel"), plot="r")
#' mantel.file(system.file("minimalexample.csv", package="mantel"), plot="r", test.args=list(maxtrials=10, conflate=TRUE))
#' mantel.file(system.file("generationsalongrows.csv", package="mantel"), plot="r", generationcolumn=4)
#' @seealso \link{\code{mantel.development}}
#' @seealso \link{\code{mantel.test}}
#' @seealso \link{\code{plot.mantel}}
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
