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
    message("Permutation space is small, enumerating all possible permutations.")
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
