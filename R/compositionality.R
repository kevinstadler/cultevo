#' An implementation of Spike's measure of compositionality.
#' @references Spike, M. 2016 \emph{...}. PhD thesis, The University of Edinburgh.
#' @param meaningmatrix a
#' @param strings a
#' @return a list
#' @export
compositionality <- function(meaningmatrix, strings) {
  
}
#compositionality.formula <- function() {
#
#}

#' Enumerate all substrings of a string.
#'
#' @param string a character string
#' @return a vector of containing all substrings of the string (including duplicates)
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

# turn strings into a substring-present-matrix
# substring counts can be greater than 1!
count.substring.occurrences <- function(strings) {
  # list of all subsequences per input string
  allstrings <- sapply(strings, enumerate.substrings, simplify=FALSE)
  substrings <- unique(unlist(allstrings))
  sapply(substrings, function(all)
    sapply(allstrings, function(sub)
      length(which(sub == all))))
}

# create a meaning-by-substring co-occurrence count matrix
# meaningmatrix has to be in 'long' binary (present/non-present) format
joint.occurrences <- function(meaningmatrix, strings) {
  if (length(strings) != nrow(meaningmatrix))
    stop("Number of strings doesn't match number of rows of the meaning matrix")

  substringmatrix <- count.substring.occurrences(strings)
  substrings <- colnames(substringmatrix)

  meanings <- colnames(meaningmatrix)
  # P(string,meaning) = P(string|meaning) * P(meaning)
  sapply(meanings,
    function(m) {
      # condition on meaning
      relevantrows <- as.logical(meaningmatrix[,m])
      # for each substring, count occurrences when that meaning was present
      apply(subset(substringmatrix, relevantrows), 2,
        function(s) {
          length(which(as.logical(s)))
        }) / sum(relevantrows)
    })
}
#joint.occurrences(cbind(a=c(0,0,1,1), b=c(1,1,0,0)), c("asd", "basd", "asdf", "basdf"))

# find the segment that maximises p(m|s) * p(s|m)
mostpredictivesegments <- function(meaningmatrix, strings,
    meaningprobabilities = colSums(meaningmatrix) / sum(meaningmatrix)) {
  jointps <- joint.occurrences(meaningmatrix, strings)
  sgivenm <- jointps / meaningprobabilities
  print(jointps)
  print(meaningprobabilities)
  print(rowSums(sgivenm))
  print(colSums(sgivenm))
  # hard
#  mgivens <- sgivenm * p(m)
}
#mostpredictivesegments(cbind(a=c(0,0,1,1), b=c(1,1,0,0)), c("asd", "basd", "asdf", "basdf"))
