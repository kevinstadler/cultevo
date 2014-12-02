#' Object class for storing the results of one or more Mantel tests
#' 
#' A data frame storing the results of one or more Mantel tests, one test per
#' row. Multiple \code{mantel} objects can easily be combined by calling
#' \code{rbind(test1, test2, ...)}
#' 
#' The data frame contains the following columns:
#' \describe{
#'    \item{\code{method}}{Character string specifying the type of correlation
#'          coefficient used (typically "pearson")}
#'    \item{\code{veridical}}{The correlation coefficient between the distances
#'          in the two distance matrices}
#'    \item{\code{sample.size}}{The number of permutations carried out}
#'    \item{\code{msample}}{A list of correlation coefficients calculated from
#'          the \code{sample.size} permutations of the input matrices}
#'    \item{\code{mean}}{Average correlation coefficient produced by the permutations}
#'    \item{\code{sd}}{Standard deviation of the correlation coefficients}
#'    \item{\code{z}}{Standard score, \code{(veridical-mean)/sd}}
#'    \item{\code{ks.level}}{p-value of the Kolmogorov-Smirnov test of
#'          normality of \code{msample}. If this value is low, \code{sd} and
#'          \code{z} might not be very meaningful.}
#'    \item{\code{p.value}}{Empirical p-value computed from the Mantel test:
#'          let \code{ngreater} be the number of correlation coefficients in
#'          \code{msample} greater than or equal to \code{veridical}, then
#'          \code{p.value} is \code{(ngreater+1)(sample.size+1}}
#'    \item{\code{p.smoothed}}{The theoretical p-value that would correspond
#'          to the standard \code{z} score as calculated above. If the
#'          \code{ks.level} is less than was specified as an argument to
#'          \link{mantel.test}, \code{p.smoothed} will be set to \code{NA}.}
#'    \item{\code{is.unique.max}}{Logical, \code{TRUE} iff the veridical
#'          correlation coefficient is greater than any of the coefficients
#'          calculated for the permutations. If this is true, then
#'          \code{p.value == 1/(sample.size+1)}}
#'  }
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
#' @param ks.level p-level below which a failure to pass the Kolmogorov-Smirnov
#'   test will result in a warning, and \code{p.smoothed} be set to \code{NA}
#' @param shuffle a function applied to \code{m2} to shuffle entries, taking a
#'   matrix and an optional permutation specification as its arguments - you
#'   shouldn't normally have to change this
#' @return a dataframe of class \code{\link{mantel}} specifying the results of
#'   the Mantel test(s)
#' @examples
#' # small distance matrix, Mantel test run deterministically
#' mantel.test(dist(1:4), dist(1:4))
#' # smallest distance matrix using random permutations (with maxtrials=1000)
#' mantel.test(dist(1:7), dist(1:7))
#' @seealso \code{|link{mantel}}, \code{|link{cor}}, \code{\link{hammingdists}},
#'   \code{\link{normalisedlevenshteindists}}, \code{\link{ks.test}}
#' @export
mantel.test <- function(m1, m2, maxtrials=1000, method="pearson", conflate=FALSE, ks.level=0.05, shuffle=shuffle.locations) {
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
  veridical <- cor(m1, m2[indices], method=method)
  # determine whether it would be more efficient to enumerate all permutations:
  # assuming we take 'maxtrials' samples, how many of the factorial(d) possible
  # values haven't been selected yet compared to the maximum number of trials?
  enumerate <- tryCatch(maxtrials/((factorial(d)*pgeom(maxtrials, 1/factorial(d), lower.tail=FALSE))) >= 1, warning=function(x)FALSE)
  if (enumerate && length(find.package("combinat", quiet=TRUE))) {
    message("Permutation space is small, enumerating all ", factorial(d), " possible permutations.")
    msample <- sapply(combinat::permn(d), function(perm)cor(m1, shuffle(m2, perm)[indices], method=method))
  } else {
    msample <- replicate(maxtrials, cor(m1, shuffle(m2)[indices], method=method))
  }
  greater <- sum(msample >= veridical)
  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/
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
    warning("The sample of randomised r's does not look normally distributed, use that z score with a grain of salt")
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
  mantel.new(method=method, veridical=veridical, mean=mn, sd=s, z=z, p.value=p, ks.level=normality$p.value, p.smoothed=p.smoothed, is.unique.max=(greater==0), sample.size=length(msample), msample=I(list(msample)))
}

#' Run consecutive Mantel tests on a set of data.
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
#' @return a dataframe of class \code{\link{mantel}} specifying the results of
#'   the Mantel test(s)
#' @examples
#' mantel.development(allmeaningcombinations(c(2,2)), c("asd", "asdf", "", "f"), plot="msample")
#' m <- matrix(c("sadasd", "iuerwh", "sdfgkj", "uofidsgf", "asd", "asdf", "", "f"), nrow=2, byrow=T)
#' mantel.development(allmeaningcombinations(c(2,2)), m, plot="r")
#' mantel.development(allmeaningcombinations(c(2,2)), m, plot="z")
#' mantel.development(allmeaningcombinations(c(2,2)), read.table(system.file("stringsonlyexample.csv", package="mantel")), plot="r")
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
#' @return a dataframe of class \code{\link{mantel}} specifying the results of
#'   the Mantel test(s)
#' @examples
#' mantel.file(system.file("minimalexample.csv", package="mantel"), plot="msample")
#' mantel.file(system.file("minimalexample.csv", package="mantel"), plot="msample", test.args=list(method="kendall", maxtrials=500))
#' mantel.file(system.file("generationsalongrows.csv", package="mantel"), plot="msample", generationcolumn=4)
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
