#' Enumerate all substrings of a string.
#'
#' @param string a character string
#' @return a vector containing all substrings of the string (including duplicates)
#' @examples
#' enumerate.substrings("abccc")
#' @export
enumerate.substrings <- function(string) {
  unlist(sapply(1 : nchar(string),
    function(start) {
      sapply(start : nchar(string),
        function(end) {
          substr(string, start, end)
        })
    }))
}

#' Count occurence of all substrings in a string.
#'
#' Returns a matrix with the original strings along rows and all substrings
#' of those strings along columns. The cell values indicate whether (and how
#' many times) the substring is contained in each of the strings.
#'
#' @param strings a list or vector of character sequences
#' @return a matrix
#' @examples
#' count.substring.occurrences(c("asd", "asdd"))
#' @export
count.substring.occurrences <- function(strings) {
  # list of all subsequences per input string
  allstrings <- sapply(strings, enumerate.substrings, simplify=FALSE)
  substrings <- unique(unlist(allstrings))
  sapply(substrings, function(all)
    sapply(allstrings, function(sub)
      length(which(sub == all))))
}

#' Spike's measure of additive compositionality.
#'
#' Implementation of the Spike-Montague information-theoretic measure of
#' additive compositionality (Spike 2016), which finds the most predictive
#' association between substrings and categorical meaning features. Additive
#' means that it does not take the ordering of elements into account, i.e.
#' \code{GREEN DOG = GREEN + DOG = DOG + GREEN}.
#'
#' The measure really captures the degree to which a synonymy-free signalling
#' system exists at the level of semantic \emph{features}, rather than looking
#' for complex meanings per se. The resulting segmentations may therefore
#' overlap.
#'
#' The segmentation algorithm scans through all sub-strings found in
#' \code{strings} to find the most bijective mapping onto all the meaning
#' features present in \code{meanings}, i.e. the pairings of sub-strings and
#' meaning features that are \emph{most predictive of each other}.
#' Mathematically, for every meaning feature \eqn{f\in M}, it tries to find
#' the sub-string \eqn{s_{ij}} from the set of strings \eqn{S} that maximises
#' \deqn{comp(f,S) = \max_{s_{ij}\in S} P(f|s_{ij}) \cdot P(s_{ij}|f)}.
#'
#' @section Meaning data format:
#' The \code{meanings} can be a matrix or data frame in one of two formats. If
#' it is a matrix of logicals (TRUE/FALSE values), then the columns are
#' assumed to refer to meaning \emph{features}, with individual cells
#' indicating whether the meaning feature is present or absent in the message
#' indicated by that row (see \code{\link{binaryfeaturematrix}} for an
#' explanation).
#' If \code{meanings} is a data frame or matrix of any other type, it is
#' assumed that the columns specify different meaning dimensions, with the
#' cell values showing the levels with which the different dimensions can be
#' realised. This dimension-based representation is automatically converted to
#' a feature-based one using \code{\link{binaryfeaturematrix}}. As a
#' consequence, whatever the actual types of the columns in the meaning
#' matrix, \emph{they will be treated as categorical factors} in order to
#' represent them as atomic meaning features.
#'
#' @describeIn sm.compositionality
#' Calculates the mean predictability of all meaning features,
#' \eqn{C = \frac{1}{|M|} \sum_{f\in M} comp(m,S)}. Returns a vector of three
#' elements: \code{N}, the number of signal-meaning pairings on which the
#' computation was based, \code{comp}, the compositionality measure, and
#' \code{M}, the number of distinct meaning features over which the measure
#' was computed. (When \code{groups} is not \code{NULL}, returns a matrix with
#' the same elements along columns, with one row for every group.)
#' @references Spike, M. (2016). Minimal requirements for the cultural
#'   evolution of language. PhD thesis, The University of Edinburgh.
#' @param strings a list or vector of character sequences
#' @param meanings a matrix or data frame with as many rows as there are
#'   strings (see below)
#' @param groups a list or vector with as many items as strings, used to split
#'   \code{strings} and \code{meanings} into data sets for which
#'   compositionality measures are computed separately.
#' @examples
#' # perfect communication system
#' sm.compositionality(c("a", "b", "ab"),
#'   cbind(a=c(T, F, T), b=c(F, T, T)))
#' sm.segmentation(c("a", "b", "ab"),
#'   cbind(a=c(T, F, T), b=c(F, T, T)))
#'
#' sm.compositionality(c("as", "bas", "basf"),
#'   cbind(a=c(T, F, T), b=c(F, T, T)))
#' sm.segmentation(c("as", "bas", "basf"),
#'   cbind(a=c(T, F, T), b=c(F, T, T)))
#'
#' # the function also accepts meaning-dimension based matrix definitions:
#' enumerate.meaningcombinations(c(animal=2, colour=2))
#' sm.segmentation(c("greendog", "bluedog", "greencat", "bluecat"),
#'   enumerate.meaningcombinations(c(animal=2, colour=2)))
#' @seealso \code{\link{binaryfeaturematrix}}
#' @export
sm.compositionality <- function(strings, meanings, groups=NULL) {
  if (is.null(groups)) {
    segmentation <- sm.segmentation(strings, meanings)
    c(N=length(strings), comp=mean(segmentation$p), M=nrow(segmentation))
  } else {
    do.call(rbind, lapply(unique(groups), function(grp) {
      c(group=grp, sm.compositionality(strings[groups==grp],
                                       meanings[groups==grp,]))
    }))
  }
}
# TODO implement formula interface

# calculate P(s_{ij}|m)
substring.predictability <- function(substringmatrix, meanings)
  sapply(colnames(meanings),
    function(m) {
      # condition on meaning
      relevantrows <- as.logical(meanings[,m])
      # for each substring, count occurrences when that meaning was present
      apply(substringmatrix[relevantrows,,drop=FALSE], 2,
        function(s) length(which(as.logical(s))) / sum(relevantrows))
    })

# Calculate the respective predictability of meaning features \eqn{f} in
# \code{meanings} given their co-occurrence with substrings, i.e. calculate
# \eqn{P(m|s_{ij})} for every \eqn{f\in M} and \eqn{s_{ij}\in S}.
meaning.predictability <- function(substringmatrix, meanings)
  t(sapply(colnames(substringmatrix),
    function(s) {
      # condition on string
      relevantrows <- as.logical(substringmatrix[,s])
      # for each substring, count occurrences when that meaning was present
      apply(meanings[relevantrows,,drop=FALSE], 2,
        function(m) length(which(as.logical(m))) / sum(relevantrows))
    }))

# meanings already have to be in binary representation!
mutualpredictability <- function(strings, meanings) {
  substringmatrix <- count.substring.occurrences(strings)
  sgivenm <- substring.predictability(substringmatrix, meanings)
  mgivens <- meaning.predictability(substringmatrix, meanings)
  sgivenm * mgivens # mutual information: I(s;m) = H(m) - H(m|s)
}

#' @describeIn sm.compositionality
#' Finds the most predictive segment(s) for every meaning feature, i.e. the
#' substrings \eqn{s_{ij}} that maximise \eqn{P(m|s_{ij}) \cdot P(s_{ij}|m)}.
#' Returns a matrix with one row for every meaning feature, in descending order
#' of their predictability from (and to) their corresponding string segments.
#' @export
sm.segmentation <- function(strings, meanings) {
  # make sure meaning matrix is in ultra-long binary format
  meanings <- binaryfeaturematrix(meanings)

  x <- mutualpredictability(strings, meanings)
  highestps <- apply(x, 2, max)
  mostpredictive <- lapply(1:ncol(x), function(i) which(x[,i] == highestps[i]))
  segmentation <- data.frame(#meaning=colnames(meanings),
    N=colSums(meanings), p=highestps,
    segments=I(lapply(mostpredictive, function(x)
      # put longest sequences first
      names(x)[order(nchar(names(x)), decreasing=TRUE)])),
    tiedsegments=sapply(mostpredictive, length))
  # put most predictable meaning/sequence pairings first
  structure(segmentation[order(-segmentation$p),], class=c("smcomp", "data.frame"))
}
