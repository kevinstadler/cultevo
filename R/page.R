#' Pre-computed lookup table for p values of Page's test for monotonicity.
#'
#' Functions for computing and assembling the table can be found in the file
#' \code{page.compute.table.R}
#' @seealso \code{\link{page.test}}
"Ltable"

#' Calculate Page's L statistic for the given dataset.
#'
#' @param data a matrix or dataframe 
#' @param ties how to resolve tied ranks. Passed on to \code{\link{rank}}, so
#'   can be one of "average" (the default), "first", "last", "random", "max",
#'   "min"
#' @param verbose print the final rankings based on which the L statistic is
#'   computed? default: \code{TRUE}
#' @seealso \code{\link{rank}}
#' @export
pages.L <- function(data, ties.method="average", verbose=TRUE) {
  # perform row-wise rankings
  ranks <- apply(data, 1, function(r) rank(r, na.last="keep", ties.method=ties.method))
  if (verbose) {
    print("Ranking used:")
    print(ranks)
  }
  sum(1:ncol(data) * rowSums(ranks))
}

# mean L for which p(x<=L) == 0.5 (as is the case when all ranks are tied)
pages.L.mean <- function(k, N=1)
  N*k*((k+1)/2)^2

page.test.normal.approx <- function(L, k, N) {
  # Hollander & Wolfe (1999) formulation
  zs <- (L - (N*k*(k+1)^2)/4) / sqrt(N * k^2 * (k+1) * (k^2-1) / 144)
  sapply(zs, function(z) pnorm(z, lower.tail=FALSE))
}

#' Perform the Page test of monotonicity on the given data set.
#'
#' The Page test tests the given matrix for monotonically increasing ranks
#' across \code{k} linearly ordered conditions (along columns) based on
#' \code{N} replications (along rows). To test for monotonically
#' \emph{decreasing} trends you can either reverse the order of columns, or
#' simply invert the ranks (i.e. by prepending the dataset with a \code{-}).
#'
#' Exact p values taken from a precalculated table are provided where possible,
#' for combinations of large \code{k} and \code{N} p values are computed based
#' on a Normal distribution approximation.
#'
#' @param data a matrix or dataframe with the different conditions along
#'   columns and the replications along rows. Conversion of the data to ranks
#'   is taken care of internally.
#' @param ... extra parameters which are passed on to \code{\link{pages.L}}
#' @return a list of class \code{pagetest}
#' @examples
#' page.test(t(replicate(4, sample(4))))
#' page.test(t(replicate(4, sample(10))))
#' @export
page.test <- function(data, ...) {
  k <- ncol(data) # 'n' in Page (1963)
  N <- nrow(data) # 'm' in Page (1963)
  if (N < 2)
    stop("Not enough chains/replications, for individual samples use the Mann-Kendall test instead")
  if (k < 3)
    stop("Need at least 3 conditions/rank levels")
  # ranking is happening in here
  L <- pages.L(data, ...)

  p <- NA
  p.str <- NULL
  # table lookup
  if (all(dim(Ltable) >= c(N,k)) && !is.null(Ltable[N,k])) {
    p <- match(TRUE, L == Ltable[N,k][[1]])
    if (is.na(p)) {
      # no exact p value available
      p.str <- ifelse(L <= pages.L.mean(k, N), " >= .5", " > .05")
    } else {
      p <- sum(Ltable[N,k][[1]]$p.ind[1:p])
    }
  } else if (k == 3) {
    # do exact calculation for k=3 (this is generally tractable)
    p <- page.combinations.bruteforce(3, N)
    p <- sum(p[1:match(TRUE, L == names(p))])
  }
  if (!is.na(p)) {
    p.str <- paste("=", format(p, scientific=TRUE), sep="")
  } else {
    warning("Exact p value is not available for given k, N, using Normal approximation which is unreliable for small N")
    approx <- page.test.normal.approx(L, k, N)
  }
  result <- list(L=L, k=k, N=N, p.exact=p, p.formatted=p.str, p.approx=approx)
  class(result) <- "pagetest"
  return(result)
}

#' @export
print.pagetest <- function(x, ...) {
  cat("Page test of monotonicity: L=", x$L, ", k=", x$k, ", N=", x$N, ", p", sep="")
  approx <- format(p$p.approx, scientific=TRUE)
  if (is.null(x$p.formatted)) {
    cat(" (approx) = ", approx, sep="")
  } else {
    cat(x$p.formatted)
    if (is.na(x$p.exact))
      cat(", p (approx) =", approx)
    else
      cat(" (exact)")
  }
  cat("\n")
  invisible(NULL)
}

#' Look up critical L values for the given significance levels from the
#' precomputed table
criticalLs <- function(Lcell, p.lvls)
  unlist(lapply(p.lvls, function(p.lvl) Lcell$L[which(cumsum(Lcell$p.ind) > p.lvl)[1]-1]))

#' Plot the goodness of the normal approximation of Page's test.
#'
#' Plots the goodness of the normal approximation of Page's test for the given
#' critical levels against the exactly computed critical Ls, where available.
#'
#' Red cells mean that the normal approximation is too liberal for the given
#' combination of \code{k}, \code{N}, blue cells mean the normal approximation
#' is too conservative. In general, we find that for low values of N the
#' approximation is slightly too \emph{liberal} at the .05 level, but too
#' \emph{conservative} for lower p-levels.
#'
#' @param p.lvls the significance levels for which the normal approximation
#'   should be compared to the exact values.
#' @example
#'par(mfrow=c(1, 2))
#'page.test.normal.approx.goodness(.05)
#'page.test.normal.approx.goodness(.01)
page.test.normal.approx.goodness <- function(p.lvls=c(0.05, 0.01, 0.001)) {
  goodness <- array(NA, c(length(p.lvls), dim(Ltable)))
  for (k in 3:ncol(Ltable)) {
    for (N in 2:nrow(Ltable)) {
      if (is.null(Ltable[N,k][[1]]))
        next
      # there might be NAs here
      critlvls <- criticalLs(Ltable[N,k][[1]], p.lvls)
      if (length(critlvls) > 0)
        # calculate difference: approximated critical L - true critical L
        suppressWarnings(goodness[,N,k] <- ceiling(qnorm(p.lvls, lower.tail=FALSE)*sqrt(N*k^2*(k+1)*(k^2-1)/144) + (N*k*(k+1)^2)/4)
          - critlvls)
    }
  }
  par(bg="lightgrey")
  image(x=1:ncol(Ltable), y=(1:(length(p.lvls)*nrow(Ltable))) / length(p.lvls),
    z=t(apply(goodness, 3, c)), xlab="k", ylab="N", main=paste("Goodness of fit of Normal approximation for p-level =", paste(p.lvls, collapse=", ")),
    xlim=c(2.5, ncol(Ltable)+0.5), ylim=c(nrow(Ltable)+0.5, 1),
    breaks=c(-10, -2, -1, 0, 1, 10), # [-Inf, -2] (-2 -1] (-1 0] (0 1] (1 Inf]
    col=rev(temperature.colors(5)))
  legend("bottomright", fill=temperature.colors(5), title="critical L derived from Normal approximation is:",
    legend=c("too conservative by >1", "too conservative by 1", "correct", "too liberal by 1", "too liberal by >1"))
}
#page.test.normal.approx.goodness(.05)

page.test.normal.approx <- function(N, k) {
  plot(Ltable[[N, k]][[1]], cumsum(Ltable[[N, k]][[2]]), type="p", pch=4,
    col="red", xlab="L", ylab="P(X>=L)", ylim=0:1, yaxs="i")
  curve(pnorm((x - (N*k*(k+1)^2)/4) / sqrt(N * k^2 * (k+1) * (k^2-1) / 144), lower.tail=FALSE), add=TRUE)
}
#page.test.normal.approx(5, 5)
