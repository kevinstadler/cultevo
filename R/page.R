#' Calculate Page's L statistic for the given dataset.
#'
#' @param data a matrix or dataframe 
#' @param ties.method how to resolve tied ranks. Passed on to \code{\link{rank}}, so
#'   can be one of "average" (the default), "first", "last", "random", "max",
#'   "min"
#' @param verbose print the final rankings based on which the L statistic is
#'   computed? default: \code{TRUE}
#' @seealso \code{\link[base]{rank}}
#' @export
page.L <- function(data, ties.method="average", verbose=TRUE) {
  # perform row-wise rankings
  ranks <- apply(data, 1, function(r) rank(r, na.last="keep", ties.method=ties.method))
  if (verbose) {
    writeLines("Ranking used:")
    print(structure(t(ranks), dimnames=list(paste("N=", 1:ncol(ranks), sep=""))))
  }
  sum(seq(ncol(data)) * rowSums(ranks))
}

# mean L for which p(x<=L) == 0.5 (as is the case when all ranks are tied)
#page.L.mean <- function(k, N=1)
#  N*k*((k+1)/2)^2

# possible row-wise ls (equivalent to Spearman's Rho)
rho.null.distribution <- function(k) {
  ls <- sum((1:k)^2) : sum(1:k * k:1)
  nd <- pspearman:::spearman.list[[k]]
  # mirror the symmetric null distribution
  if (length(ls) %% 2 == 0)
    c(nd, rev(nd))
  else
    c(nd, rev(nd)[-1])
  # sanity check: sum(rho.null.distribution(k)) == factorial(k)
}

L.null.distribution <- function(k, N) {
  # possible row-wise ls
  ls <- sum((1:k)^2) : sum(1:k * k:1)
  convolutions <- list(rho.null.distribution(k))
  names(convolutions[[1]]) <- ls

  currentls <- ls
  i <- 1
  while (i < N) {
    currentls <- (currentls[1] + ls[1]) : (currentls[length(currentls)] + ls[length(ls)])
    # calculate all joint probabilities
    jointp <- outer(convolutions[[i]], convolutions[[1]])
    # pad matrix so that equal sums of L are aligned in rows
    # (resulting matrix will have nrow(jointp)+ncol(jointp)-1 rows)
    jointp <- sapply(1:ncol(jointp), function (col) {
      c(rep(0, col-1), jointp[,col], rep(0, max(0, ncol(jointp)-col)))
    })
    # it's that easy
    newconvolution <- rowSums(jointp)
    names(newconvolution) <- currentls

    i <- i+1
    convolutions[[i]] <- newconvolution
  }
  return(convolutions[[N]] / sum(convolutions[[N]]))
}
# don't recompute
if (requireNamespace("memoise", quietly = TRUE))
  L.null.distribution <- memoise::memoise(L.null.distribution)

#' Calculate exact significance levels of the Page L statistic.
#'
#' Returns the null probability that the Page statistic with the given k, N is >= L.
#'
#' @param k number of conditions/generations
#' @param N number of replications/chains
#' @param L value of the Page L statistic
#' @return the significance level of L for the given k, N
#' @examples
#' page.compute.exact(6, 4, 322)
#' @seealso \code{\link{page.test}}
#' @export
page.compute.exact <- function(k, N, L) {
  nd <- L.null.distribution(k, N)
  sum(nd[1:match(L, names(nd))])
}

page.compute.normal.approx <- function(k, N, L) {
  # Hollander & Wolfe (1999) formulation
  zs <- (L - (N*k*(k+1)^2)/4) / sqrt(N * k^2 * (k+1) * (k^2-1) / 144)
  sapply(zs, function(z) stats::pnorm(z, lower.tail=FALSE))
}
# the Normal approximation of crank::page.trend.test is demonstrably off
# ps <- replicate(50, {
#   k <- 12
#   N <- 4
#   d <- t(replicate(N, sample(k)))
#   c(cultevo=page.compute.normal.approx(k, N, page.L(d, verbose=FALSE)),
#     crank=crank::page.trend.test(d)$pZ)
# })
#hist(ps[2,], col="red")
#hist(ps[1,], col="blue", add=TRUE)

#' Perform the Page test of monotonicity.
#'
#' The Page test tests the given matrix for monotonically increasing ranks
#' across \code{k} linearly ordered conditions (along columns) based on
#' \code{N} replications (along rows). To test for monotonically
#' \emph{decreasing} trends you can either reverse the order of columns, or
#' simply invert the ranks (i.e. by prepending the dataset with a \code{-}).
#'
#' Exact p values are computed for \code{k} up to 22, using the pre-computed
#' null distributions from the \code{pspearman} package. For larger \code{k}, p
#' values are computed based on a Normal distribution approximation.
#'
#' @param data a matrix or dataframe with the different conditions along
#'   columns and the replications along rows. Conversion of the data to ranks
#'   is taken care of internally (see \code{\link{page.L}})
#' @param ... extra parameters which are passed on to \code{\link{page.L}}
#' @return a list of class \code{pagetest} (and \code{htest}) containing the following fields:
#' \describe{
#'    \item{\code{L}}{Value of the \code{L} statistic for the given data set}
#'    \item{\code{k}}{Number of conditions (columns of the data set)}
#'    \item{\code{N}}{Number of replications (rows of the data set)}
#'    \item{\code{p.value}}{Significance level}
#'    \item{\code{p.type}}{Whether the computed p value is "exact" or "approx"}
#' }
#' @examples
#' page.test(t(replicate(4, sample(4))))
#' page.test(t(replicate(4, sample(10))))
#' @export
page.test <- function(data, ...) {
  k <- ncol(data) # 'n' in Page (1963)
  N <- nrow(data) # 'm' in Page (1963)
  if (N < 2)
    stop("Not enough chains/replications, for individual samples use the Mann-Kendall test instead")
    # http://finzi.psych.upenn.edu/R/library/EnvStats/html/kendallSeasonalTrendTest.html
  if (k < 3)
    stop("Need at least 3 conditions/rank levels")
  # ranking is happening in here
  L <- page.L(data, ...)

  # if exact null distribution is known: compute exact
  if (k <= length(pspearman:::spearman.list)) {
    p <- page.compute.exact(k, N, L)
    p.type <- "exact"
  } else {
    warning("Exact p value is not available for given k, using Normal approximation which is unreliable for small N")
    p <- page.compute.normal.approx(k, N, L)
    p.type <- "approx"
  }
  structure(list(statistic=L, k=k, N=N, p.value=p, p.type=p.type,
      alternative="", method="Page test of monotonicity"),
    class=c("htest", "pagetest"))
}

#' @export
print.pagetest <- function(x, ...)
  invisible(cat("Page test of monotonicity: L=", x$statistic, ", k=", x$k,
    ", N=", x$N, "\np = ",
    format(x$p.value, digits=4), " (", x$p.type, ")\n", sep=""))

# =============================================================================
# Below here are functions used to pre-compute the exact distribution of the
# Page L statistic for individual rows (rankings). This is currently tractable
# for k up to 14
# =============================================================================

# calculate the multinomial coefficient (number of permutations of a multiset,
# i.e. the number of possible orders in which n elements can be drawn when some
# of the individual elements are identical. might overflow for large n/k.

# ks == vector of positive integers specifying the multiplicities of the
# individual elements. 1 elements can be omitted.
mcombn <- function(n, ks)
#  if (!all.equal(c(n, ks), as.integer(c(n, ks))) || sum(ks) > n)
#    stop("All arguments must be integer with sum(ks) == n")
  factorial(n) / prod(sapply(ks, factorial))

# Calculates the probabilities of one replication (row) producing a certain L
# under the null hypothesis. Loops through all factorial(k) possible rankings,
# calculates their row-wise Ls and returns a table of frequencies (ordered by
# descending L)
rowwise.ls <- function(k) {
  ncombinations <- factorial(k)
#  ls <- numeric(ncombinations)
  comb <- 1:k
  maxl <- sum(comb * comb)
  minl <- sum(comb * rev(comb))
  ls <- numeric(1+maxl-minl)
  ls[1] <- 1

  # helper variables to detect opportunity for half-way interruption
  meanl <- (minl+maxl)/2
  highls <- 1
  meanls <- 0
  # Heap's non-recursive algorithm for enumerating all rank permutations
  # cf. http://www.cs.princeton.edu/~rs/talks/perms.pdf
  cs <- rep(1, k)
  n <- 1
  while (n <= k) {
    if (cs[n] < n) {
      if (n %% 2) {
        comb[c(1,n)] <- comb[c(n,1)] # odd recursion level
      } else {
        comb[c(cs[n],n)] <- comb[c(n,cs[n])] # even recursion level
      }
      cs[n] <- cs[n]+1
      n <- 1
      l <- sum(1:k * comb)
      ls[1+maxl-l] <- ls[1+maxl-l] + 1

      # probability space is symmetric, so can fold around half-way point
      if (l > meanl) {
        highls <- highls + 1
      } else if (l == meanl) {
        meanls <- meanls + 1
      }
      if (2*highls+meanls == ncombinations) {
        halfpoint <- length(ls)/2
        ls[length(ls):ceiling(halfpoint+1)] <- ls[1:floor(halfpoint)]
        break
      }
    } else {
      cs[n] <- 1
      n <- n+1
    }
  }
  # sanity check that Heap's worked: sum(ls) == factorial(k)
  # divide by total count to get probabilities
  names(ls) <- maxl:minl
  ls/ncombinations
}
