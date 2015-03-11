pagesL <- function(data, ties="average")
  sum(1:ncol(data) * rowSums(apply(data, 1, function(r)rank(r, na.last="keep", ties.method=ties))))

normal.approx.p <- function(L, k, N) {
  if (N <= 10)
    warning("Normal approximation is only reliable for N>10")
  # Hollander & Wolfe (1999) formulation
  zs <- (L-(N*k*(k+1)^2)/4) / sqrt(N*k^2*(k+1)*(k^2-1)/144)
  sapply(zs, function(z)c(p=pnorm(z, lower.tail=FALSE), z=z))
}

.onLoad <- function(libname, pkgname) {
  data("Ltable", package=pkgname, envir=parent.env(environment()))
}

#' Perform Page's test on the given data set.
#'
#' Page's test tests for monotonically increasing ranks across linearly ordered
#' conditions. To test for monotonically \emph{decreasing} trends you can
#' either reverse the order of columns, or simply invert the ranks (i.e. by
#' prepending a \code{-}).
#'
#' @param data a matrix or dataframe with the different conditions along
#'   columns and the replications along rows. Conversion to ranks is taken care
#'   of internally.
#' @examples
#' pages.test(t(replicate(4, sample(4))))
#' pages.test(t(replicate(4, sample(10))))
#' @export
pages.test <- function(data, ...) {
  k <- ncol(data) # 'n' in Page (1963)
  N <- nrow(data) # 'm' in Page (1963)
  if (N==1)
    stop("Only one sample/block, you should be running a Mann-Kendall test instead, e.g.: Kendall::MannKendall(data)")
  if (k<3)
    stop("Need at least 3 rank levels")
  # ranking is happening in here
  L <- pagesL(data, ...)
  # 1. lookup table
  p <- NA
  p.str <- NULL
  if (all(dim(Ltable) >= c(N,k)) && !is.null(Ltable[N,k])) {
    p <- match(TRUE, L == Ltable[N,k][[1]])
    if (!is.na(p)) {
      p <- Ltable[N,k][[1]]$p[p]
      p.str <- paste("=", p, sep="")
    } else {
      p.str <- if (L <= N*sum(1:k)^2/k) " > .5" else " > .05"
    }
  }
  approx <- suppressWarnings(normal.approx.p(L,k,N))
  result <- list(L=L, k=k, N=N, p=p, p.str=p.str, p.approx=approx[1], z.L=approx[2])
  class(result) <- "pages.test"
  return(result)
}

#' @export
print.pages.test <- function(p) {
  cat("Page's test: L=", p$L, ", k=", p$k, ", N=", p$N, ", p", sep="")
  if (p$p.approx <= 0.0005)
    approx <- p$p.approx
  else
    approx <- round(p$p.approx, digits=4)
  if (is.na(p$p.str))
    cat(" (approx) = ", approx, sep="")
  else
    cat(p$p.str)
    if (is.na(p$p))
      cat(", p (approx) =", approx)
    else
      cat(" (exact)")
  cat("\n")
}
