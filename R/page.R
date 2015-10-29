.onLoad <- function(libname, pkgname) {
  data("Ltable", package=pkgname, envir=parent.env(environment()))
}

#' Calculate Page's L statistic for the given dataset.
#'
#' @param data a matrix or dataframe 
#' @param ties how to resolve tied ranks. Passed on to \code{\link{rank}}, so can
#'   be one of "average" (the default), "first", "last", "random", "max", "min"
#' @seealso \code{\link{rank}}
#' @export
pagesL <- function(data, ties.method="average")
  sum(1:ncol(data) * rowSums(apply(data, 1, function(r)rank(r, na.last="keep", ties.method=ties.method))))

# mean L for which p(x<=L)=0.5 (same as having all ranks tied)
meanL <- function(k, N=1)
  N*k*((k+1)/2)^2

pages.test.normal.approx <- function(L, k, N) {
  if (N <= 10)
    warning("Normal approximation is only reliable for N>10")
  # Hollander & Wolfe (1999) formulation
  zs <- (L-(N*k*(k+1)^2)/4) / sqrt(N*k^2*(k+1)*(k^2-1)/144)
  sapply(zs, function(z)c(p=pnorm(z, lower.tail=FALSE), z=z))
}

#' Perform Page's test on the given data set, providing exact p values whre possible.
#'
#' Page's test tests for monotonically increasing ranks across linearly ordered
#' conditions. To test for monotonically \emph{decreasing} trends you can
#' either reverse the order of columns, or simply invert the ranks (i.e. by
#' prepending the dataset with a \code{-}).
#'
#' @param data a matrix or dataframe with the different conditions along
#'   columns and the replications along rows. Conversion of the data to ranks
#'   is taken care of internally.
#' @param ... extra parameters which are passed on to \code{\link{pagesL}}
#' @return a list of class \code{PagesTest}
#' @examples
#' pages.test(t(replicate(4, sample(4))))
#' pages.test(t(replicate(4, sample(10))))
#' @export
pages.test <- function(data, ...) {
  k <- ncol(data) # 'n' in Page (1963)
  N <- nrow(data) # 'm' in Page (1963)
  if (N==1)
    stop("Only one sample/block/replication, you should be running a Mann-Kendall test instead, e.g.: Kendall::MannKendall(data)")
  if (k<3)
    stop("Need at least 3 conditions/rank levels")
  # ranking is happening in here
  L <- pagesL(data, ...)
  # table lookup
  p <- NA
  p.str <- NULL
  if (all(dim(Ltable) >= c(N,k)) && !is.null(Ltable[N,k])) {
    p <- match(TRUE, L == Ltable[N,k][[1]])
    if (is.na(p)) {
      # no exact p value available
      p.str <- ifelse(L <= meanL(k,N), " >= .5", " > .05")
    } else {
      p <- sum(Ltable[N,k][[1]]$p.ind[1:p])
    }
  } else if (k == 3) {
    # do exact calculation for k=3
    p <- page.combinations.bruteforce(3, N)
    p <- sum(p[1:match(TRUE, L == names(p))])
  }
  if (!is.na(p)) {
    p.str <- paste("=", format(p, scientific=TRUE), sep="")
  }
  approx <- suppressWarnings(pages.test.normal.approx(L, k, N))
  result <- list(L=L, k=k, N=N, p.exact=p, p.str=p.str, p.approx=approx[1], z.L=approx[2])
  class(result) <- "PagesTest"
  return(result)
}

print.PagesTest <- function(p) {
  cat("Page's test: L=", p$L, ", k=", p$k, ", N=", p$N, ", p", sep="")
  approx <- format(p$p.approx, scientific=TRUE)
  if (is.null(p$p.str)) {
    cat(" (approx) = ", approx, sep="")
  } else {
    cat(p$p.str)
    if (is.na(p$p.exact))
      cat(", p (approx) =", approx)
    else
      cat(" (exact)")
  }
  cat("\n")
}
