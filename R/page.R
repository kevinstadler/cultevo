#' Page test for monotonicity of ranks.
#'
#' Given \code{N} replications of \code{k} different treatments/conditions,
#' tests whether the \emph{median ordinal ranks} \eqn{m_i} of the treatments
#' are identical \deqn{m_1 = m_2 = \ldots = m_k} against the alternative
#' hypothesis \deqn{m_1 \leq m_2 \leq \ldots \leq m_k} where \emph{at least
#' one} of the inequalities is a strict inequality (Siegel and Castellan
#' 1988, p.184). Given that even a single point change in the distribution of
#' ranks across conditions represents evidence against the null hypothesis,
#' the Page test is simply a test for \emph{some ordered differences in
#' ranks}, but not a 'trend test' in any meaningful way (see also the
#' [Page test tutorial](https://kevinstadler.github.io/cultevo/articles/page.test.html)).
#'
#' Tests the given matrix for monotonically *increasing* ranks across `k`
#' linearly ordered conditions (along columns) based on `N` replications
#' (along rows). To test for monotonically *decreasing* ranks, either reverse
#' the order of columns, or simply invert the rank ordering by calling `-` on
#' the entire dataset.
#'
#' Exact p-values are computed for `k` up to 22, using the pre-computed null
#' distributions from the
#' [`pspearman`](https://CRAN.R-project.org/package=pspearman) package. For
#' larger `k`, p-values are computed based on a Normal distribution
#' approximation (Siegel and Castellan, 1988).
#'
#' @param data a matrix with the different conditions along its \code{k}
#'   columns and the \code{N} replications along rows. Conversion of the data
#'   to ordinal ranks is taken care of internally.
#' @param verbose whether to print the final rankings based on which the L
#'   statistic is computed
#' @return \code{page.test} returns a list of class \code{pagetest} (and
#'   \code{htest}) containing the following elements:
#' \describe{
#'    \item{\code{statistic}}{value of the L statistic for the data set}
#'    \item{\code{parameter}}{a named vector specifying the number of
#'      conditions (k) and replications (N) of the data (which is the number of
#'      columns and rows of the data set, respectively)}
#'    \item{\code{p.value}}{significance level}
#'    \item{\code{p.type}}{whether the computed p-value is `"exact"` or
#'      `"approximate"`}
#' }
#' @examples
#' # exact p value computation for N=4, k=4
#' page.test(t(replicate(4, sample(4))))
#' 
#' # exact p value computation for N=4, k=10
#' page.test(t(replicate(4, sample(10))))
#' 
#' # approximate p value computation for N=4, k=23
#' result <- page.test(t(replicate(4, sample(23))), verbose = FALSE)
#' 
#' print(result)
#' 
#' @references Siegel, S., and N. J. Castellan, Jr. (1988). Nonparametric
#' Statistics for the Behavioral Sciences. McGraw-Hill.
#' @describeIn page.test See above.
#' @export
page.test <- function(data, verbose=TRUE) {
  k <- ncol(data) # 'n' in Page (1963)
  N <- nrow(data) # 'm' in Page (1963)
  if (N < 2)
    stop("Not enough chains/replications, for individual samples use the Mann-Kendall test instead")
    # http://finzi.psych.upenn.edu/R/library/EnvStats/html/kendallSeasonalTrendTest.html
  # ranking is happening in here
  L <- page.L(data, verbose=verbose)

  # if exact null distribution is known: compute exact
  if (k <= 22) {
    p <- page.compute.exact(k, N, L)
    p.type <- "exact"
  } else {
    message("Exact p value is not available for k=", k, ", using Normal approximation.")
    p <- page.compute.normal.approx(k, N, L)
    p.type <- "approximate"
  }
  structure(list(statistic=c(L=L), parameter=c(k=k, N=N), p.value=p,
      p.type=p.type, data.name=paste(deparse(substitute(data))),
      method="Page test of monotonicity of ranks",
      alternative="at least one difference in ranks"),
    class=c("htest", "pagetest"))
}

#' @describeIn page.test
#' Calculate Page's L statistic for the given dataset.
#'
#' @param ties.method how to resolve tied ranks. Passed on to
#' \code{\link[base]{rank}}, should be left on "average" (the default).
#' @seealso \code{\link[base]{rank}}, [Page test tutorial](https://kevinstadler.github.io/cultevo/articles/page.test.html)
#' @export
page.L <- function(data, verbose=TRUE, ties.method="average") {
  # perform row-wise rankings
  ranks <- apply(data, 1, function(r) rank(r, na.last="keep", ties.method=ties.method))
  if (verbose) {
    writeLines("Ranking used:")
    print(structure(t(ranks), dimnames=list(paste("N=", 1:ncol(ranks), sep=""))))
  }
  sum(seq(ncol(data)) * rowSums(ranks))
}

# mean L for which p(x<=L) == 0.5 (as is the case when all ranks are tied)
page.L.mean <- function(k, N)
  N*k*((k+1)/2)^2

# possible row-wise ls (equivalent to Spearman's Rho)
rho.null.distribution <- function(k) {
  # access pspearman:::spearman.list null distribution directly
#  ls <- sum((1:k)^2) : sum(1:k * k:1)
#  nd <- pspearman:::spearman.list[[k]]
  # mirror the symmetric null distribution
#  if (length(ls) %% 2 == 0)
#    c(nd, rev(nd))
#  else
#    c(nd, rev(nd)[-1])
  # sanity check: sum(rho.null.distribution(k)) == factorial(k)

  # reconstruct null distribution from calls to pspearman()
  nls <- sum((1:k)^2) - sum(1:k * k:1)
  ps <- sapply(2*0:floor(nls/2), function(x) pspearman::pspearman(x, k))
  ps <- c(ps[1], diff(ps))
  if (nls %% 2 == 0)
    c(ps, rev(ps)[-1])
  else
    c(ps, rev(ps))
}
# don't recompute
if (requireNamespace("memoise", quietly = TRUE))
  rho.null.distribution <- memoise::memoise(rho.null.distribution)

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

#' @describeIn page.test Calculate exact significance levels of the Page L
#' statistic. Returns a single numeric indicating the null probability of
#' the Page statistic with the given `k`, `N` being greater or equal than the
#' given `L`.
#'
#' @param k number of conditions/generations
#' @param N number of replications/chains
#' @param L value of the Page L statistic
#' @examples
#' # raw calculation of the significance levels
#' page.compute.exact(6, 4, 322)
#' @export
page.compute.exact <- function(k, N, L) {
  if (k < 2)
    stop("Need at least 2 conditions/rank levels")
  range <- c(sum(seq(k) * N*k:1), sum(seq(k) * N*1:k))
  if (L < range[1] || L > range[2])
    stop("Valid range of L for k=", k, ", N=", N, ": [",
      paste(range, collapse=", "), "]")
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

#' @export
print.pagetest <- function(x, ...)
  invisible(cat("Page test of monotonicity: L=", x$statistic, ", k=", x$k,
    ", N=", x$N, "\np = ",
    format(x$p.value, digits=4), " (", x$p.type, ")\n", sep=""))
