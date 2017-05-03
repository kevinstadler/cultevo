#' Perform one or more Mantel permutation tests.
#' 
#' Perform correlation tests between pairs of matrices. The Mantel test is
#' different from classical correlation tests (such as
#' \code{\link[stats]{cor.test}}) in that the null distribution of the
#' correlation coefficient is determined empirically by shuffling the locations
#' in one of the matrices and calculating the resulting correlations to
#' generate an empirical null distribution for the given data set.
#' 
#' If the number of possible permutations of the matrices is reasonably close
#' to the number of permutations specified by the \code{trials} parameter, a
#' deterministic enumeration of all the permutations will be carried out
#' instead of random sampling.
#'
#' @param x a formula, distance matrix, or list of distance matrices (see below)
#' @param y a data frame, distance matrix, or list of matrices of the same length as \code{x}
#' @param plot logical: immediately produce a plot of the test results (using
#'   \code{plot()})
#' @param method correlation coefficient to be computed (passed on to
#'   \code{\link[stats]{cor}}, default: "pearson")
#' @param trials maximum number of permutations to be tested
#' @param omitzerodistances logical: if \code{TRUE}, the calculation of the
#'   correlation coefficient omits pairs of off-diagonal cells which contain a
#'   0 in the first distance matrix argument. (For the formula interface, this
#'   is the matrix which specifies the meaning distances.)
#' @param groups when \code{x} is a formula: column name by which the data in
#'   \code{y} is split into separate data sets to run several Mantel tests on
#' @param stringdistfun when \code{x} is a formula: edit distance function used
#'   to compute the distance matrix from the specified string column. Supports
#'   any edit distance function that returns a distance matrix from a vector or
#'   list of character strings. Default is Levenshtein distance (function
#'   \code{\link[utils]{adist}}), other options from this package include
#'   \code{\link{normalisedlevenshteindists}} and
#'   \code{\link{orderinsensitivedists}}.
#' @param meaningdistfun when \code{x} is a formula: meaning distance function
#'   used to compute the distance matrix from the specified meaning columns.
#'   Defaults to Hamming distances between meanings (function
#'   \code{\link{hammingdists}}), custom meaning functions can be created
#'   easily using \code{\link{wrap.meaningdistfunction}}.
#' @param ... further arguments which are passed on to the default method (in
#'   particular \code{plot}, \code{method}, \code{trials} and
#'   \code{omitzerodistances})
#' @return A dataframe of class \code{mantel}, with one row per Mantel
#' test carried out, containing the following columns:
#' \describe{
#'    \item{\code{method}}{Character string specifying the type of correlation
#'          coefficient used ("spearman", "kendall", or, inadvisable in the
#'          likely case of ties: "pearson")}
#'    \item{\code{statistic}}{The veridical correlation coefficient between
#'          the entries in the two distance matrices}
#'    \item{\code{rsample}}{A list of correlation coefficients calculated
#'          from the permutations of the input matrices}
#'    \item{\code{mean}}{Average correlation coefficient produced by the permutations}
#'    \item{\code{sd}}{Standard deviation of the sampled correlation coefficients}
#'    \item{\code{p.value}}{Empirical p-value computed from the Mantel
#'          test: let \code{ngreater} be the number of correlation coefficients
#'          in \code{rsample} greater than or equal to \code{statistic}, then
#'          \code{p.value} is \code{(ngreater+1)/(length(rsample)+1}}
#'    \item{\code{p.approx}}{The theoretical p-value that would correspond
#'          to the standard \code{z} score as calculated above.}
#'    \item{\code{is.unique.max}}{Logical, \code{TRUE} iff the veridical
#'          correlation coefficient is greater than any of the coefficients
#'          calculated for the permutations. If this is true, then
#'          \code{p.value == 1 / (length(rsample)+1)}}
#'  }
#' Multiple \code{mantel} objects can easily be combined by calling
#' \code{rbind(test1, test2, ...)}. 
#' @examples
#' # small distance matrix, Mantel test run deterministically
#' mantel.test(dist(1:7), dist(1:7))
#'
#' # smallest distance matrix using random permutations (based on default trials=9999)
#' plot(mantel.test(dist(1:8), dist(1:8), method="kendall"))
#'
#' mantel.test(hammingdists(enumerate.meaningcombinations(c(2, 2, 2, 2))),
#'   dist(1:16), plot=TRUE)
#' @seealso \code{\link[stats]{cor}},
#'   \code{\link[utils]{adist}}, \code{\link{hammingdists}},
#'   \code{\link{normalisedlevenshteindists}},
#'   \code{\link{orderinsensitivedists}}
#' @export
mantel.test <- function(x, y, ...)
  UseMethod("mantel.test")

#' @describeIn mantel.test Perform Mantel correlation test on two distance
#' matrices. The distance matrices can either be of type
#' \code{\link[stats]{dist}}, plain R matrices or any object that can be
#' interpreted by \code{\link{check.dist}}.
#' @importFrom stats cor
mantel.test.default <- function(x, y, plot=FALSE, method="pearson", trials=9999, omitzerodistances=FALSE, ...) {
  m1 <- check.dist(x)
  m2 <- check.dist(y)
  d <- dim(m1)[1]
  if (dim(m2)[1] != d)
    stop("The two distance matrices have incompatible dimensions")

  # figure out actual complexity by counting *unique* combinations in rows:
  # check if a row i is equivalent to a row j which has i+j swapped (in which
  # case the permutation of the two does not lead to different r)
#  swaps <- combinat::combn(d, 2)
#  inconsequentialswaps <- apply(swaps, 2, function(i) {
#    !any(m1[i[1],] != m2[i[2], c(seq_len(i[1]-1), i[2], i[1]+seq_len(i[2]-i[1]-1), i[1], i[2]+seq_len(d-i[2]))])
#    })
#  if (any(inconsequentialswaps))
#    print(swaps[,inconsequentialswaps])

  # extract values relevent for correlation computation
  indices <- which(lower.tri(m1))
  if (omitzerodistances)
    indices <- setdiff(indices, which(m1 == 0))
  m1 <- m1[indices]
  duration <- max(0.001, system.time(veridical <- cor(m1, m2[indices], method=method))[[3]])

  # determine whether it would be more efficient to enumerate all permutations:
  # assuming we take 'trials' samples, how many of the factorial(d) possible
  # values haven't been selected yet compared to the maximum number of trials?
  exact <- tryCatch(trials/((factorial(d)*stats::pgeom(trials, 1/factorial(d), lower.tail=FALSE))) >= 1,
    warning=function(x) FALSE)
  if (exact)
    message("Permutation space is small, enumerating all ", factorial(d), " possible permutations.")
  if (duration*trials > 30)
    message("Estimated time to evaluate all ", trials, " permutations is ", duration*trials, "s, go get yourself a biscuit!")

  if (exact) {
    rsample <- sapply(combinat::permn(d), function(perm) cor(m1, shuffle.locations(m2, perm)[indices], method=method))
  } else {
    rsample <- replicate(trials, cor(m1, shuffle.locations(m2)[indices], method=method))
    # calculate lower bound by sampling local permutations
    # all pairwise swaps
    switches <- combinat::combn(d, 2)
    localrs <- apply(switches, 2, function(i) {
      cor(m1,
        shuffle.locations(m2,
          c(seq_len(i[1]-1), i[2], i[1]+seq_len(i[2]-i[1]), i[1], i[2]+seq_len(d-i[2])))[indices],
        method=method)})
    p.lowerbound <- sum(localrs >= veridical) / factorial(d)
  }

  greater <- sum(rsample >= veridical)
  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/
  p.empirical <- (greater + 1) / (length(rsample) + 1)

  mn <- mean(rsample)
  s <- stats::sd(rsample)

  result <- data.frame(method=method, statistic=c(r=veridical),
    N=length(indices), mean=mn, sd=s, p.value=p.empirical,
    p.approx=ifelse(exact, NA,
      stats::pnorm((veridical - mn) / s, lower.tail=FALSE)),
    is.unique.max=(greater==0), alternative="greater",
    rsample=I(list(rsample)), stringsAsFactors=FALSE)
  if (plot)
    plot.mantel(result)
  structure(result, class=c("mantel", "htest", "data.frame"))
}

#' @export
mantel.test.dist <- mantel.test.default

#' @describeIn mantel.test This function can be called with experimental result
#' data frames, distance matrix calculation is taken care of internally.
#' \code{x} is a formula of the type \code{s ~ m1 + m2 + ...} where \code{s}
#' specifies a column of character strings in data frame or matrix \code{y},
#' while \code{m1} etc. are columns specifying the different meaning
#' dimensions. To calculate the respective distances, the function
#' \code{stringdistfun} is applied to the strings, \code{meaningdistfun} to the
#' meaning columns.
#' @examples
#' # using the formula interface in combination with a data frame
#' mantel.test(word ~ Var1 + Var2, cbind(word=c("aa", "ab", "ba", "bb"),
#'   enumerate.meaningcombinations(c(2, 2))))
#' @export
mantel.test.formula <- function(x, y, groups=NULL, stringdistfun=utils::adist,
  meaningdistfun=hammingdists, ...) {
  t <- stats::terms(x)
  fields <- rownames(attr(t, "factors"))
  lhs <- fields[attr(t, "response")]
  rhs <- attr(t, "term.labels")
  # TODO implement grouping/conditioning term parsing for formula interface
  if (is.null(groups)) {
    # pass meanings first
    mantel.test.default(meaningdistfun(y[,rhs]), stringdistfun(y[,lhs]), ...)
  } else {
    levels <- sort(unique(y[,groups]))
    mantel.test.list(sapply(levels, function(lvl)
          meaningdistfun(y[which(y[[groups]] == lvl), rhs]),
        simplify=FALSE),
      sapply(levels, function(lvl)
          stringdistfun(y[which(y[[groups]] == lvl), lhs]),
        simplify=FALSE), ...)
  }
}

#' @describeIn mantel.test When \code{x} is a list of distance matrices, and
#' \code{y} is either a single distance matrix or a list of distance matrices
#' the same length as \code{x}: runs a Mantel test for every pairwise
#' combination of distance matrices in \code{x} and \code{y} and returns a
#' \code{mantel} object with as many rows.
#' @examples
#'
#' # running tests on a list of distance matrices and plotting them
#' plot(mantel.test(list(dist(1:8), dist(sample(8:1)), dist(runif(8))),
#'   hammingdists(enumerate.meaningcombinations(c(2, 2, 2)))), xlab="group")
#' @export
mantel.test.list <- function(x, y, plot=FALSE, ...) {
  # catch plot argument so that individual tests don't all produce a plot
  result <- do.call(rbind, mapply(mantel.test.default, x,
    # save single distance matrices y from being iterated into with list(y)
    if (is.list(y)) y else list(y), plot=FALSE, MoreArgs=..., SIMPLIFY=FALSE))
  result <- cbind(group=if (is.null(names(x))) seq_along(x) else names(x), result)
  if (plot)
    plot.mantel(result)
  structure(result, class=c("mantel", "data.frame"))
}

#' @export
print.mantel <- function(x, ...) {
  for (i in 1:nrow(x))
    cat("Mantel permutation test (method: ", x$method[i], ")\nr = ",
      format(x$statistic[i], digits=3), ", N = ", x$N[i], "\n",
      length(x$rsample[[i]]), " permutations, mean = ",
      format(x$mean[i], digits=3), ", sd = ", format(x$sd[i], digits=3),
      "\np (", if(is.na(x$p.approx[i])) "exact" else "empirical",
      ") = ", format(x$p.value[i], digits=3),
      ifelse(x$is.unique.max[i], " (veridical correlation is highest found)",
        ""), "\n\n", sep="")
  # TODO print approximate p value iff not is.unique.max
  invisible(x)
}

# functions for visualising the outcome of Mantel tests
blue.if.true <- function(x, default="black")
  mapply(function(x, d) if (x) "blue" else d, x, default)

method.label <- list(pearson="r", kendall=expression(tau), spearman=expression(rho))

# grab the xaxt, xlab, xlim and ylim arguments to stop them from being passed on via ...
plotmantelsample <- function(mantels, nbins=25, main="", xaxt=NULL, xlab=NULL, xlim=NULL, ylim=NULL, ...) {
  d <- list()
  if (nrow(mantels) > 1)
    graphics::par(mfrow=c(1, nrow(mantels)), mar=c(5.1, 2.5, 2, 1))
    # layout(matrix(c(rep(1, nrow(mantels)), 2:(nrow(mantels)+1)), nrow=2, byrow=TRUE))

  maxdensity <- 0
  for (i in 1:nrow(mantels)) {
    d[[i]] <- graphics::hist(unlist(mantels$rsample[i]), breaks=nbins, plot=FALSE)
    maxdensity <- max(maxdensity, d[[i]]$density, stats::dnorm(0, sd=mantels$sd[i]))
  }
  for (i in 1:nrow(mantels)) {
    xlim <- range(mantels$rsample[[i]], mantels$statistic[i]) + c(-.05, .05)
    graphics::plot(d[[i]], freq=FALSE, yaxs="i", xlim=xlim, ylim=c(0, maxdensity), xlab=method.label[[mantels$method[i]]], ylab="Density", border="gray", main=paste(method.label[[mantels$method[i]]], "=", format(mantels$statistic[i], digits=3), ", N=", mantels$N[i], ", ", length(mantels$rsample[[i]]), " permutations", sep=""), ...)
    # add fit used for z score estimation
    graphics::curve({stats::dnorm(x, mean=mantels$mean[i], sd=mantels$sd[i])}, lty=3, add=TRUE)
    # mark veridical r
    col <- blue.if.true(mantels$is.unique.max[i])
    level <- stats::dnorm(mantels$statistic[i], mean=mantels$mean[i], sd=mantels$sd[i])
    graphics::segments(mantels$statistic[i], 0, y1=level, col=col, lty=2)
    graphics::points(mantels$statistic[i], y=0.15+level, pch=25, bg=col, col=col)
    graphics::text(mantels$statistic[i], 0.15+level, paste("p=", format(mantels$p.value[i], digits=3), sep=""), pos=3)
  }
}

#' @details
#' \cr\code{plot()} called on a data frame of class \code{mantel} plots a
#' visualisation of the test results (in particular, the distribution of
#' the permutated samples against the veridical correlation coefficient). If
#' the veridical correlation coefficient is plotted in blue it means
#' that it was higher than all other coefficients generated by random
#' permutations of the data. When the argument contains the result of more than
#' two Mantel tests, a side-by-side boxplot visualisation of the samples is
#' produced (see examples). Additional parameters \code{...} to \code{plot()}
#' are passed on to \code{\link[graphics]{plot.default}}.
#' @param xlab the x axis label used when plotting the result of several Mantel
#'   tests next to each other
#' @rdname mantel.test
#' @export
plot.mantel <- function(x, xlab="generation", ...) {
  if (nrow(x) <= 2) {
    plotmantelsample(x)
  } else {
    # TODO replace with boxplots
    Hmisc::errbar(1:nrow(x), x$mean, x$mean + x$sd, x$mean - x$sd, xlab=xlab,
      xaxt="n", ylab=paste(method.label[[x$method[1]]], " (mean+-sd)", sep=""),
      xlim=c(0.5, nrow(x) + 0.5),
      ylim=range(0.1, x$statistic, x$mean + x$sd, x$mean - x$sd)+c(-0.1, 0.1), 
      main=paste("N=", length(x$rsample[[1]]), sep=""), ...)
    # plot correlation coefficient reference points
    graphics::abline(h=-1:1, lty=c(3, 2, 3), col="grey")
    # blue points signify that no single larger r value has been sampled
    graphics::points(1:nrow(x), x$statistic, col=blue.if.true(x$is.unique.max))
    # label the veridical rs with their p values
    graphics::text(1:nrow(x), ifelse(x$statistic >= x$mean,
      pmax(x$statistic, x$mean + x$sd), pmin(x$statistic, x$mean - x$sd)),
      labels=paste("p", format(x$p.value, digits=2), sep="="), pos=2+sign(x$statistic))
    graphics::axis(1, at=1:nrow(x), labels=rownames(x))
  }
}
