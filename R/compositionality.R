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

#' Count occurences of all possible substrings in one more strings.
#'
#' @param strings a list or vector of character sequences
#' @param sortbylength logical indicating whether the substring columns should
#'   be ordered according to the (decreasing) length of the substrings. Default
#'   is to leave them in the original order in which they occur in the given
#'   strings.
#' @return A matrix with the original strings along rows and all substrings
#'   of those strings along columns. The cell values indicate whether (and how
#'   many times) the substring is contained in each of the strings.
#' @examples
#' count.substring.occurrences(c("asd", "asdd", "foo"))
#' @export
count.substring.occurrences <- function(strings, sortbylength=FALSE) {
  # list of all subsequences per input string
  allstrings <- sapply(strings, enumerate.substrings, simplify=FALSE)
  substrings <- unique(unlist(allstrings))
  if (sortbylength)
    substrings <- substrings[order(sapply(substrings, nchar), decreasing=TRUE)]
  sapply(substrings, function(all)
    sapply(allstrings, function(sub)
      length(which(sub == all))))
}

#' @importFrom stats weighted.mean
weighted.mean.mp <- function(segmentation)
  weighted.mean(segmentation$mp, segmentation$N)

#' Spike's segmentation and measure of additive compositionality.
#'
#' Implementation of the Spike-Montague segmentation and measure of additive
#' compositionality (Spike 2016), which finds the most predictive associations
#' between meaning features and substrings. Computation is deterministic and
#' fast.
#'
#' The algorithm works on compositional meanings that can be expressed as sets
#' of categorical meaning features (see below), and does not take the order
#' of elements into account. Rather than looking directly at how complex
#' meanings are expressed, the measure really captures the degree to which a
#' homonymy- and synonymy-free signalling system exists at the level of
#' *individual semantic features*.
#'
#' The segmentation algorithm provided by `sm.segmentation()` scans through
#' all sub-strings found in `strings` to find the pairings of meaning features
#' and sub-strings whose respective presence is *most predictive of each
#' other*. Mathematically, for every meaning feature \eqn{f\in M}, it finds
#' the sub-string \eqn{s_{ij}} from the set of strings \eqn{S} that yields the
#' highest *mutual predictability* across all signals,
#' \deqn{mp(f,S) = \max_{s_{ij}\in S}\ P(f|s_{ij}) \cdot P(s_{ij}|f)\;.}
#'
#' Based on the mutual predictability levels obtained for the individual
#' meaning features, `sm.compositionality` then computes the mean mutual
#' predictability weighted by the individual features' relative frequencies of
#' attestation, i.e.
#' \deqn{mp(M,S) = \sum_{f\in M} freq_f \cdot mp(f,S)\;,}
#' as a measure of the overall compositionality of the signalling system.
#'
#' Since mutual predictability is determined seperately for every meaning
#' feature, the most predictive sub-strings posited for different meaning
#' features as returned by `sm.segmentation()` can overlap, and even coincide
#' completely. Such results are generally indicative of either limited data
#' (in particular frequent co-occurrence of the meaning features in question),
#' or spurious results in the absence of a consistent signalling system. The
#' latter will also be indicated by the significance level of the given mutual
#' predictability.
#'
#' @section Null distribution and p-value calculation:
#' A perfectly unambiguous mapping between a meaning feature to a specific
#' string segment will always yield a mutual predictability of `1`. In the
#' absence of such a regular mapping, on the other hand, chance co-occurrences
#' of strings and meanings will in most cases stop the mutual predictability
#' from going all the way down to `0`. In order to help distinguish chance
#' co-occurrence levels from significant signal-meaning associations,
#' `sm.segmentation()` provides significance levels for the mutual
#' predictability levels obtained for each meaning feature.
#'
#' What is the baseline level of association between a meaning feature and a
#' set of sub-strings that we would expect to be due to chance co-occurrences?
#' This depends on several factors, from the number of data points on which the
#' analysis is based to the frequency of the meaning feature in question and,
#' perhaps most importantly, the overall makeup of the different substrings
#' that are present in the signals. Since every substring attested in the data
#' is a candidate for signalling the presence of a meaning feature, the
#' absolute number of different substrings greatly affects the likelihood of
#' chance signal-meaning associations. (Diversity of the set of substrings is
#' in turn heavily influenced by the size of the underlying alphabet, a factor
#' which is often not appreciated.)
#'
#' For every candidate substring, the degree of association with a specific
#' meaning feature that we would expect by chance is again dependent on the
#' absolute number of signals in which the substring is attested.
#'
#' Starting from the simplest case, take a meaning that is featured in \eqn{m}
#' of the total \eqn{n} signals (where \eqn{0 < m \leq n}). Assume next that
#' there is a string segment that is attested in \eqn{s} of these signals
#' (where again \eqn{0 < s \leq n}). The degree of association between the
#' meaning feature and string segment is dependent on the number of times that
#' they co-occur, which can be no more than \eqn{c_{max} = min(m, s)} times.
#' The null probability of getting a given number of co-occurrences can be
#' obtained by considering all possible reshufflings of the meaning feature in
#' question across all signals: if \eqn{s} signals contain a given substring,
#' how many of \eqn{s} randomly drawn signals from the pool of \eqn{n} signals
#' would contain the meaning feature if a total of \eqn{m} signals in the pool
#' did? Approached from this angle, the likelihood of the number of
#' co-occurrences follows the
#' [hypergeometric distribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution),
#' with \eqn{c} being the number of successes when taking \eqn{s} draws without
#' replacement from a population of size \eqn{n} with fixed number of successes
#' \eqn{m}.
#'
#' For every number of co-occurrences \eqn{c \in [0, c_{max}]}, one can
#' compute the corresponding mutual probability level as
#' \eqn{p(c|s) \cdot p(c|m)} to obtain the null distribution of mutual
#' predictability levels between a meaning feature and *one* substring of a
#' particular frequency \eqn{s}:
#' \deqn{Pr(mp = p(c|s) \cdot p(c|m)) = f(k=c; N=n, K=m, n=s)}
#'
#' From this, we can now derive the null distribution for the entire set of
#' attested substrings as follows: making the simplifying assumption that the
#' occurrences of different substrings are independent of each other, we first
#' aggregate over the null distributions of all the individual substrings to
#' obtain the mean probability \eqn{p=Pr(X\ge mp)} of finding a given mutual
#' predictability level at least as high as \eqn{mp} for one randomly drawn
#' string from the entire population of substrings. Assuming the total number
#' of candidate substrings is \eqn{|S|}, the overall null probability that at
#' least one of them would yield a mutual predictability at least as high is
#' \deqn{Pr(X\ge 0), X \equiv B(n=|S|, p=p)\;.}
#'
#' Note that, since the null distribution also depends on the frequency with
#' which the meaning feature is attested, the significance levels corresponding
#' to a given mutual predictability level are not necessarily identical for
#' all meaning features, even within one analysis.
#' 
#' (In theory, one can also compute an overall p-value of the weighted mean
#' mutual predictability as calculated by `sm.compositionality`. However, the
#' significance levels for the individual meaning features are much more
#' insightful and should therefore be consulted directly.)
#'
#' @section Meaning data format:
#' The `meanings` argument can be a matrix or data frame in one of two formats.
#' If it is a matrix of logicals (`TRUE`/`FALSE` values), then the columns are
#' assumed to refer to meaning *features*, with individual cells indicating
#' whether the meaning feature is present or absent in the signal represented
#' by that row (see [binaryfeaturematrix()] for an explanation). If `meanings`
#' is a data frame or matrix of any other type, it is assumed that the columns
#' specify different meaning dimensions, with the cell values showing the
#' levels with which the different dimensions can be realised. This
#' dimension-based representation is automatically converted to a
#' feature-based one by calling [binaryfeaturematrix()]. As a consequence,
#' whatever the actual types of the columns in the meaning matrix, *they will
#' be treated as categorical factors* for the purpose of this algorithm, also
#' discarding any explicit knowledge of which 'meaning dimension' they might
#' belong to.
#'
#' @references Spike, M. 2016 \emph{Minimal requirements for the cultural
#'   evolution of language}. PhD thesis, The University of Edinburgh.
#'   <http://hdl.handle.net/1842/25930>.
#' @param x a list or vector of character sequences specifying the signals to
#'   be analysed. Alternatively, \code{x} can also be a formula of the format
#'   \code{s ~ m1 + m2 + ...}, where \code{s} and \code{m1}, \code{m2}, etc.
#'   specify the column names of the signals and meaning features found in the
#'   data frame that is passed as the second argument.
#' @param y a matrix or data frame with as many rows as there are signals,
#'   indicating the presence/value of the different meaning dimensions along
#'   columns (see section Meaning data format). If \code{x} is a formula, the
#'   \code{y} data frame can contain any number of columns, but only the ones
#'   whose column name is specified in the formula will be considered.
#' @param groups a list or vector with as many items as strings, used to split
#'   \code{strings} and \code{meanings} into data sets for which
#'   compositionality measures are computed separately.
#' @param strict logical: if \code{TRUE}, perform additional filtering of
#'   candidate segments. In particular, it removes combinations of segments
#'   (across meanings) which overlap in at least one of the strings where they
#'   co-occur. For convenience, it also removes segments which are shorter
#'   substrings of longer candidates (for the same meaning feature).
#' @return \code{sm.segmentation} provides detailed information about the most
#'   predictably co-occurring segments for every meaning feature. It returns
#'   a data frame with one row for every meaning feature, in descending order
#'   of the mutual predictability from (and to) their corresponding string
#'   segments. The data frame has the following columns:
#'   \describe{
#'     \item{\code{N}}{The number of signal-meaning pairings in which this
#'       meaning feature was attested.}
#'     \item{\code{mp}}{The highest mutual predictability between this
#'       meaning feature and one (or more) segments that was found.}
#'     \item{\code{p}}{Significance levels of the given mutual predictability,
#'       i.e. the probability that the given mutual predictability level could
#'       be reached by chance. The calculation depends on the frequency of the
#'       meaning feature as well as the number and relative frequency of all
#'       substrings across all signals (see below).}
#'     \item{\code{ties}}{The number of substrings found in \code{strings}
#'       which have this same level of mutual predictability with the meaning
#'       feature.}
#'     \item{\code{segments}}{For \code{strict=FALSE}: a list containing the
#'       \code{ties} substrings in descending order of their length (the
#'       ordering is for convenience only and not inherently meaningful). When
#'       \code{strict=TRUE}, the lists of segments for each meaning feature
#'       are all of the same length, with a meaningful relationship of the
#'       order of segments across the different rows: every set of segments
#'       which are found in the same position for each of the different
#'       meaning features constitute a valid segmentation where the segments
#'       occurrences in the actual signals do not overlap.}
#'   }
#'
#'   \code{sm.compositionality} calculates the weighted average of the
#'   mutual predictability of all meaning features and their most predictably
#'   co-occurring strings, as computed by \code{sm.segmentation}. The function
#'   returns a data frame of three columns:
#'   `N` is the total number of signals (utterances) on which the computation
#'   was based, `M` the number of distinct meaning features attested across
#'   all signals, and `meanmp` the mean mutual predictability across all these
#'   features, weighted by the features' relative frequency. When `groups` is
#'   not `NULL`, the data frame contains one row for every group.
#' @examples
#' # perfect communication system for two meaning features (which are marked
#' # as either present or absent)
#' sm.compositionality(c("a", "b", "ab"),
#'   cbind(a=c(TRUE, FALSE, TRUE), b=c(FALSE, TRUE, TRUE)))
#' sm.segmentation(c("a", "b", "ab"),
#'   cbind(a=c(TRUE, FALSE, TRUE), b=c(FALSE, TRUE, TRUE)))
#'
#' # not quite perfect communication system
#' sm.compositionality(c("as", "bas", "basf"),
#'   cbind(a=c(TRUE, FALSE, TRUE), b=c(FALSE, TRUE, TRUE)))
#' sm.segmentation(c("as", "bas", "basf"),
#'   cbind(a=c(TRUE, FALSE, TRUE), b=c(FALSE, TRUE, TRUE)))
#'
#' # same communication system, but force candidate segments to be non-overlapping
#' # via the 'strict' option
#' sm.segmentation(c("as", "bas", "basf"),
#'   cbind(a=c(TRUE, FALSE, TRUE), b=c(FALSE, TRUE, TRUE)), strict=TRUE)
#'
#'
#' # the function also accepts meaning-dimension based matrix definitions:
#' print(twobytwoanimals <- enumerate.meaningcombinations(c(animal=2, colour=2)))
#'
#' # note how there are many more candidate segments than just the full length
#' # ones. the less data we have, the more likely it is that shorter substrings
#' # will be just as predictable as the full segments that contain them.
#' sm.segmentation(c("greendog", "bluedog", "greencat", "bluecat"), twobytwoanimals)
#'
#' # perform the same analysis, but using the formula interface
#' print(twobytwosignalingsystem <- cbind(twobytwoanimals,
#'   signal=c("greendog", "bluedog", "greencat", "bluecat")))
#'
#' sm.segmentation(signal ~ colour + animal, twobytwosignalingsystem)
#'
#' # since there is no overlap in the constituent characters of the identified
#' # 'morphemes', they are all tied in their mutual predictiveness with the
#' # (shorter) substrings they contain
#' #
#' # to reduce the pool of candidate segments to those which are
#' # non-overlapping and of maximal length, again use the 'strict=TRUE' option:
#' 
#' sm.segmentation(signal ~ colour + animal, twobytwosignalingsystem, strict=TRUE)
#' 
#' @seealso [binaryfeaturematrix()], [ssm.compositionality()]
#' @importFrom stats terms
#' @export
sm.compositionality <- function(x, y, groups=NULL, strict=FALSE)
  UseMethod("sm.compositionality")

#' @export
sm.compositionality.formula <- function(x, y, groups=NULL, strict=FALSE)
  runfromformula(sm.compositionality.default, x, y, groups, strict)
#   t <- stats::terms(x)
#   fields <- rownames(attr(t, "factors"))
#   strings <- y[ , fields[attr(t, "response")]]
#   if (!is.null(groups))
#     groups <- y[ , groups]
#   meanings <- y[ , attr(t, "term.labels")]
#   sm.compositionality.default(strings, meanings, groups, strict)
# }

#' @export
sm.compositionality.default <- function(x, y, groups=NULL, strict=FALSE) {
  # TODO make sure it's a data frame and data/columns are well-formed?
  if (is.null(groups)) {
    segmentation <- sm.segmentation(x, y, strict=strict)
    data.frame(N=length(x), M=nrow(segmentation), meanmp=weighted.mean.mp(segmentation))
    # TODO if (strict) cbind(rate=mean(segmentation$rate)) # m=sum(), mr=sum(), 
  } else {
    do.call(rbind, lapply(unique(groups), function(grp) {
      # FIXME preserve type of group (numeric/factor)?
      cbind(group=grp, sm.compositionality(x[groups == grp],
        y[groups == grp, ], strict=strict))
    }))
  }
}

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

# take a mutual predictability matrix and select most predictive segments,
# ordering meaning features by decreasing predictability
topsegments <- function(x) {
  # find top mutual predictiveness per meaning
  mps <- apply(x, 2, max)
  # get indices of segments tied for top predictiveness
  mostpredictive <- lapply(1:ncol(x), function(i) which(x[,i] == mps[i]))

  # get tied segments and order by length
  segments <- lapply(mostpredictive, function(is)
    # could use rownames(x)[is] instead?
    names(is)[order(nchar(names(is)), decreasing=TRUE)])
  data.frame(mp=mps, segments=I(segments))
}

# we have a selection of candidate ranges where all relevant segments are.
# check that there's at least one selection of ranges without overlaps
has.nonoverlapping.ranges <- function(pos)
  is.null(find.selections(pos, function(boundaries) {
      # order by first (start) row
      boundaries <- boundaries[, order(boundaries[1, ])]
      # do any next-segment starts precede the previous' ends?
      return(any(boundaries[1, -1] <= boundaries[2, -ncol(boundaries)]))
    }, first=TRUE, noptions=sapply(pos, nrow),
    selectionfun=function(p, i) p[i, ]))

# check if the given segmentation produces overlap in any of the strings
# meanings is a logical matrix of |strings| rows and |segments| columns
#' @importFrom stringi stri_locate_all_fixed
is.overlapping <- function(strings, meanings, segments) {
  for (i in seq_along(strings)) {
    # find segment positions in string
    pos <- stri_locate_all_fixed(strings[i], segments[meanings[i,]])
    # segments might be unattested in this string
    pos <- pos[!sapply(pos, function(p) is.na(p[1]))]

    # just one meaning? trivially non-overlapping
    if (length(pos) < 2) # sum(meanings[i,])
      next
    if (!has.nonoverlapping.ranges(pos))
      return(TRUE)
  }
  FALSE
}

#' @rdname sm.compositionality
#' @export
sm.segmentation <- function(x, y, strict=FALSE)
  UseMethod("sm.segmentation")

#' @export
sm.segmentation.formula <- function(x, y, ...)
  runfromformula(sm.segmentation.default, x, y, ...)

#' @importFrom stats terms
#' @export
sm.segmentation.default <- function(x, y, strict=FALSE) {

  # make sure meaning matrix is in ultra-long binary format
  meanings <- binaryfeaturematrix(y, x)
  meaningfrequencies <- colSums(meanings)
  substringmatrix <- count.substring.occurrences(x)

  sgivenm <- substring.predictability(substringmatrix, meanings)
  mgivens <- meaning.predictability(substringmatrix, meanings)
  mp <- sgivenm * mgivens

  meaninglabels <- colnames(meanings)
#  meanings <- unname(meanings)

  top <- topsegments(mp)
  segments <- top$segments

  # determine significance levels
  p <- mp.plvls(length(x), substringmatrix, meaningfrequencies, top$mp)

  if (strict) {
    # find selections of segments which don't overlap in individual strings.
    message("Applying strict selection, checking ",
      prod(sapply(top$segments, length)), " segment combinations for overlap")
    nonoverlapping <- find.selections(top$segments,
      function(elements) !is.overlapping(x, meanings, elements))
    # if we simply consider overlap as non-occurrence of one of the segments we
    # would have to re-calculate the mutual predictability for non-overlapping
    # selections of segments (see mutualpredictability.strict) which might
    # lower the score massively. the mutual predictability approach by itself
    # is also problematic when we further try to improve the score by merging
    # currently non-predictive features together since greedily maximising the
    # overall *mean* mutual predictability at every step will not stop until
    # all categories are merged.
    segments <- mapply('[', top$segments, lapply(seq(ncol(nonoverlapping)),
      function(i) nonoverlapping[,i]))

    if (is.null(segments)) {
      stop("No combination of non-overlapping per-meaning segments found among most predictive segment set. Note that this does not necessarily mean that no such segmentation exists. Try running with strict=FALSE.")
    }
    # more than one candidate segmentation
    if (is.matrix(segments)) {
      # limit to longest combination of candidate segments across all slots
      # (just a heuristic! doesn't mean it's actually maximising string
      # coverage across signals...)
      segmentlengths <- apply(segments, 1, function(ss) sum(nchar(ss)))
      # also make into list appropriate for merging into data frame below
      longest <- which(segmentlengths == max(segmentlengths))
      if (length(longest) == 1) {
        segments <- segments[longest,]
      } else {
        segments <- lapply(which(segmentlengths == max(segmentlengths)),
          function(row) segments[row,])
      }
    }

    # nonoverlapping <- find.selections(replicate(ncol(mp),
    #     rownames(x), simplify=FALSE),
    #   function(elements) ! is.overlapping(x, meanings, elements))
  }

  # TODO if meanings was in wide format: include dimensions+values
  structure(data.frame(row.names=meaninglabels,
    N=meaningfrequencies, mp=top$mp, p=p,
    ties=sapply(segments, length),
    segments=I(segments))[order(-top$mp, p), ], class=c("sm", "data.frame"),
    nsignals=length(x), ntotalchars=sum(sapply(x, nchar)),
    signalentropy=entropy(meaningfrequencies/length(x)),
    # mutual information: I(S;M) = H(M) - H(M|S) where
    # H(M|S) = - SUM[S] p(s) * SUM[M] p(m|s)*log(p(m|s))
    mi=entropy(meaningfrequencies/length(x)) +
      weighted.mean(colSums(sgivenm*log2(sgivenm)),
        meaningfrequencies/length(x)))
}

entropy <- function(p)
  -sum(p * log2(p))

# n = number of signals
# m = number of signals in which meaning feature is attested
# stringfreqs = density = number of segments with frequency of index (1 and up)
#' @importFrom stats aggregate dhyper
#' @importFrom utils head
mp.null.distribution <- function(n, m, stringfreqs) {
  # calculate weighted probability distribution over mp levels
  ps <- do.call(rbind, lapply(seq_along(stringfreqs), function(s) {
    # s = number of signals in which substring is attested.
    # between 0 and min(m,s) of the meaning-relevant sigs can contain the segment
    hits <- 0:min(m, s)
    # calculate mp values for each (trivially, mp = 0 when they never co-occur)
    mps <- (hits / s) * (hits / m) # p(m|s) * p(s|m)
    # the likelihood of each number of 'hits' occurring follows a
    # hypergeometric distribution (sampling from a population with fixed
    # number of successes)
    # this should be P(mp < X) so that we can multiply them all later on
    cbind(s=s, mp=mps, p=dhyper(hits, m, n-m, s))
  }))
  # to create appropriate weighted probabilities over all string frequencies,
  # multiply every element with its stringfreq, then divide by sum(stringfreqs)
  ps <- cbind(ps, weightedp=ps[,'p'] * stringfreqs[ ps[,'s'] ] / sum(stringfreqs))
  # aggregate (by addition) based on mp level to get global null distribution
  ps <- aggregate(weightedp ~ mp, ps, sum)
  # return P(mp < X)
  data.frame(mp=ps$mp, p=c(0, cumsum(head(ps$weightedp, -1))))
}
#mp.null.distribution(10, 2, 5:1)

# don't recompute
#if (requireNamespace("memoise", quietly = TRUE))
#  mp.null.distribution <- memoise::memoise(mp.null.distribution)

#' @importFrom stats dbinom
mp.null.probability <- function(n, m, stringfreqs, mp) {
  nulldist <- mp.null.distribution(n, m, stringfreqs)
  # strict matching, safer but might mess up because floating point
  i <- match(mp, nulldist$mp)
  # based on this P(mp < X) for one string, what's the probability that *all*
  # available strings produce a mp<X? then take complement of that for P(mp>=X)
  1 - dbinom(sum(stringfreqs), sum(stringfreqs), nulldist$p[i])
}
#mp.null.probability(10, 2, 5:1, 1/4)

# convenience function: meaningfreqs and mps need to be of the same length
mp.plvls <- function(nsignals, substringmatrix, meaningfrequencies, mps)
  sapply(seq_along(meaningfrequencies), function(i)
    mp.null.probability(nsignals, meaningfrequencies[i],
      tabulate(colSums(substringmatrix > 0)), mps[i]))

#' Find a segmentation that maximises the overall string coverage across all signals.
#'
#' This algorithm builds on Spike's measure of compositionality (see
#' \code{\link{sm.compositionality}}), except instead of simply determining
#' which segment(s) have the highest mutual predictability for each
#' meaning feature separately, it attempts to find a combination of
#' non-overlapping segments for each feature that maximises the overall string
#' coverage over all signals. In other words, it tries to find a segmentation
#' which can account for (or 'explain') as much of the string material in the
#' signals as possible.
#'
#' For large data sets and long strings, this computation can get very slow.
#' If the attested signals are such that no perfect segmentation is possible,
#' this algorithm is not guaranteed to find any segmentation (as no such
#' segmentation might exist).
#'
#' @param x a list or vector of character sequences
#' @param y a matrix or data frame with as many rows as there are
#'   strings (see section Meaning data format)
#' @param groups a list or vector with as many items as strings, used to split
#'   the signals and meanings into data sets for which the compositionality
#'   measures are computed separately.
#' @param mergefeatures logical: if \code{TRUE}, \code{ssm.segmentation} will
#'   try to improve on the initial solution by incrementally merging pairs of
#'   meaning features as long as doing so improves the overall string coverage
#'   of the segmentation.
#' @param verbose logical: if \code{TRUE}, messages detailed information about
#'   the number of segment combinations considered for every coverage computed.
#' @examples
#' ssm.segmentation(c("as", "bas", "basf"),
#'   cbind(a=c(TRUE, FALSE, TRUE), b=c(FALSE, TRUE, TRUE)))
#'
#'
#' # signaling system where one meaning distinction is not encoded in the signals
#' print(threebytwoanimals <- enumerate.meaningcombinations(list(animal=c("dog", "cat", "tiger"),
#'   colour=c("col1", "col2"))))
#'
#' ssm.segmentation(c("greendog", "bluedog", "greenfeline", "bluefeline", "greenfeline", "bluefeline"),
#'   threebytwoanimals)
#'
#' # the same analysis again, but allow merging of features
#' ssm.segmentation(c("greendog", "bluedog", "greenfeline", "bluefeline", "greenfeline", "bluefeline"),
#'   threebytwoanimals, mergefeatures=TRUE)
#' @seealso \code{\link{sm.compositionality}}
#' @export
ssm.compositionality <- function(x, y, groups=NULL)
  UseMethod("ssm.compositionality")

#' @export
ssm.compositionality.formula <- function(x, y, groups=NULL)
  runfromformula(ssm.compositionality.default, x, y, groups)

#' @export
ssm.compositionality.default <- function(x, y, groups=NULL) {
  if (is.null(groups)) {
    segmentation <- ssm.segmentation.default(x, y)
    data.frame(N=length(x), M=ncol(y),
      chars=sum(sapply(x, nchar)),
      covered=0, # TODO
      attr(segmentation, "weightedsignalcoverage"),
      attr(segmentation, "meansignalcoverage"))
  } else {
    do.call(rbind, lapply(unique(groups), function(grp) {
      cbind(group=grp, ssm.compositionality.default(x[groups == grp],
                                                    y[groups == grp,]))
    }))
  }
}

#' @rdname ssm.compositionality
#' @export
ssm.segmentation <- function(x, y, mergefeatures=FALSE, verbose=FALSE)
  UseMethod("ssm.segmentation")

#' @export
ssm.segmentation.formula <- function(x, y, mergefeatures=FALSE, verbose=FALSE)
  runfromformula(ssm.segmentation.default, x, y, mergefeatures, verbose)

#' @importFrom utils tail
#' @export
ssm.segmentation.default <- function(x, y, mergefeatures=FALSE, verbose=FALSE) {
  totalchars <- sum(sapply(x, nchar))
  # make sure meaning matrix is in ultra-long binary format
  meanings <- binaryfeaturematrix(y)

  coverage <- maximise.stringcoverage(x, meanings)
  message("Initial segmentation covers ", coverage$charscovered, " of ", totalchars, " characters, mean mp ", round(coverage$meanmp, 3))

  # based on this first selection of segments we can also consider (always
  # trying to maximise coverage$meanstringcoverage):
  # 1. less predictive segments of the current segmentations (in the limit
  #   considering all substrings as candidates for all meanings).

  if (mergefeatures) {
    # 2. merge currently non-predictive features together.
    # a heuristic based on incrementally merging the least predictive
    # category and greedily maximising the overall *mean* mutual predictability
    # at every step will not stop until all categories are merged.
    # instead: find combination of segments that maximises string coverage

    # as long as we still distinguish *some* features...
    # TODO fix substring predictability to work with one feature, then do >= 2
    while (ncol(meanings) > 2) {
      mergers <- as.matrix(combinat::combn(ncol(meanings), 2))
      if (verbose)
        message("Trying ", ncol(mergers),
          " different combinations of meaning feature mergers")
      # TODO parallel::mclapply() here, plus limiting to within-dim mergers
      mergers <- apply(mergers, 2, function(merge) {
        mergeddata <- cbind(meanings[, -merge], apply(meanings[, merge], 1, any))
        colnames(mergeddata) <- c(colnames(meanings)[-merge],
          paste(colnames(meanings)[merge], collapse="|"))
        list(meanings=mergeddata,
          segmentation=maximise.stringcoverage(x, mergeddata, coverage$charscovered, verbose=verbose))
      })
      # don't consider merges which didn't at least match in string coverage
      mergers <- mergers[sapply(mergers, function(y) !is.null(y$segmentation))]
      # don't consider merges which didn't improve mean mutual predictability
      if (length(mergers) == 0)
        break

      mergers <- mergers[sapply(mergers, function(y) y$segmentation[["meanmp"]] > coverage$meanmp)]
      # maybe also discard mergers which would reduce the meanrate?

      # greedily select merger which maximises string coverage (doesn't work
      # as well when going for highest mean mutual predictability)
      highest <- which.max(sapply(mergers, function(merger)
        merger$segmentation[["charscovered"]])) # charscovered or meanrate
      if (length(highest)) {
        meanings <- mergers[[highest]]$meanings
        coverage <- mergers[[highest]]$segmentation
        message("Merging ", tail(colnames(meanings), n=1),
          " improved coverage to ", coverage$charscovered, " out of ",
          totalchars, ", mean mp ", round(coverage$meanmp, 3))
      } else {
        break
      }
    }
  }

  meaninglabels <- colnames(meanings)
  meaningfrequencies <- colSums(meanings)
  meaningrealisations <- colSums(coverage$meanings, na.rm=TRUE)
  p <- mp.plvls(length(x), coverage$substringmatrix, meaningfrequencies, coverage$mps)

  # TODO if meanings was in wide format: include dimensions+values
  structure(data.frame(row.names=meaninglabels,
      N=meaningfrequencies, matches=meaningrealisations,
      matchrate=meaningrealisations / meaningfrequencies,
      mp=coverage$mps, #ties=sapply(segments, length),
      p=p,
      segments=I(coverage$segments)) #reorder
        [order(-coverage$mps, -meaningrealisations / meaningfrequencies), ],
    class=c("ssm", "data.frame"),
    nsignals=length(x),
    ntotalchars=totalchars,
    charscovered=sum(coverage$charscoveredperstring),
    meansignalcoverage=coverage$meancoverage,
    weightedsignalcoverage=coverage$weightedcoverage)
}

# calculate the mean string coverage of the given segments for the given
# strings, discounting overlaps
#' @importFrom stringi stri_locate_all_fixed
stringcoverage <- function(strings, meanings, segments) {
  realisedsegments <- sapply(seq_along(strings), function(i) {
    pos <- stri_locate_all_fixed(strings[i], segments[meanings[i,]])
    # mark segments that are altogether missing/unattested
    missing <- sapply(pos, function(p) is.na(p[1]))
    if (any(missing))
      meanings[i, which(meanings[i, ])[missing] ] <- NA

    # next check for overlap in attested segments
    if (sum(!missing) >= 2 && !has.nonoverlapping.ranges(pos[!missing])) {
      # TODO only set minimum number necessary to avoid overlap to NA
      # strict (allaround) punishment of any overlap
      meanings[i, which(meanings[i, ])] <- NA
    }

    meanings[i,]
  })

  # calculate coverage (in absolute characters)
  lengths <- sapply(segments, nchar)
  charscovered <- sapply(seq_along(strings), function(i)
    # strings/signals are along columns this time!
    sum(lengths[which(realisedsegments[ , i])]))
  list(meaningscoveredperstring=colSums(realisedsegments, na.rm=TRUE),
    charscoveredperstring=charscovered,
    meanings=I(t(realisedsegments)))
}

#' @importFrom stats weighted.mean
#' @importFrom utils head
# @export
maximise.stringcoverage <- function(strings, meanings, charscovered=-1,
    maxcombinations=500, verbose=TRUE) {
  substringmatrix <- count.substring.occurrences(strings, TRUE)
  segments <- colnames(substringmatrix)
  # FIXME this breaks if meanings is unnamed earlier...
  sgivenm <- substring.predictability(substringmatrix, meanings)
  meanings <- unname(meanings)
  mgivens <- meaning.predictability(substringmatrix, meanings)
  x <- sgivenm * mgivens

  # find combination of segments that maximises string coverage (similar task
  # to finding the highest strict mutual predictability)
  predorder <- apply(-x, 2, function(xs) rank(xs, ties.method="min"))
  # only look at top-scoring (in mutual predictability) segments..
  segmentids <- lapply(1:ncol(x), function(i) which(predorder[,i] == 1))

  meaningfrequencies <- colSums(meanings)
  segmentlengths <- sapply(segments, nchar)

  # figure out if we can do full enumeration...
  if (prod(sapply(segmentids, length)) > maxcombinations) {
    # make a clever selection based on which has a) most ties b) lowest mp?
    # just try first/longest segment per feature
    ids <- rbind(rep(1, ncol(meanings)),
      # also try any tied top segments which have the highest sgivenm
      sapply(seq_along(segmentids),
        function(i) which.max(sgivenm[segmentids[[i]],i])))
    segmentids <- lapply(seq_along(segmentids),
      function(i) segmentids[[i]][unique(ids[,i])])
  }
  if (verbose)
    message("Checking ", prod(sapply(segmentids, length)), " segment combinations for overlaps...")
  sels <- unname(as.matrix(do.call(expand.grid, segmentids)))

  # ways to optimise:
  # 1. pre-compute upper limit
  # upper bound for coverage = sum_m( nchar(s) * p(s|m) * p(m) )
  upperbounds <- apply(sels, 1, function(sel)
    sum(segmentlengths[sel] *
      sgivenm[cbind(sel, seq_along(sel))] *
      meaningfrequencies))

  # 2. start with longest (theoretical) coverage strings
  selectionid <- 0
  coverage <- NULL
  ord <- order(upperbounds, decreasing=TRUE)
  j <- 0
  for (i in ord) {
    j <- j+1
    # accept equal coverage too, mutual predictability might be higher
    if (upperbounds[i] < charscovered) {
#      message("Interrupting after ", j)
      break
    }
  # 3. memoise some of the overlap functions?
    nextcoverage <- stringcoverage(strings, meanings, segments[sels[i, ]])
    nextcharscovered <- sum(nextcoverage$charscoveredperstring)
    if (nextcharscovered >= charscovered) {
      selectionid <- i
      coverage <- nextcoverage
      charscovered <- nextcharscovered
    }
  }

  if (!is.null(coverage)) {
    coverage$substringmatrix <- substringmatrix
    coverage$sels <- sels[selectionid, ]
    coverage$segments <- segments[coverage$sels]
    coverage$charscovered <- charscovered
    coverage$mps <- sapply(seq_along(coverage$sels), function(i)
      x[coverage$sels[i], i])
    coverage$meanmp <- weighted.mean(coverage$mps, colSums(meanings))
    # feature-frequency-weighted mean rate == overall hit rate
    coverage$meanrate <-
      sum(coverage$meanings, na.rm=TRUE) / sum(meaningfrequencies)
    perstringcoverage <- coverage$charscoveredperstring/sapply(strings, nchar)
    coverage$meancoverage <- mean(perstringcoverage)
    # weight by number of features encoded (this is often uniform anyway)
    coverage$weightedcoverage <- weighted.mean(perstringcoverage, rowSums(meanings))
    coverage
  }
}

#' @export
print.sm <- function(x, digits=3, ...) {
  if (!is.null(digits)) {
    d <- getOption("digits")
    options(digits=digits)
  }
  # format p values nicely
  if (!is.null(x$p))
    x$p <- pvalue.str(x$p)
  print.data.frame(x)
  cat("\nMean feature-wise mutual predictability, weighted by feature frequency:", weighted.mean.mp(x), "\n")
  if (!is.null(digits))
    options(digits=d)
  invisible(x)
}

#' @export
print.ssm <- function(x, digits=3, ...) {
  print.sm(x, digits, ...)
#  cat("Mean signal-wise character coverage, unweighted:",
#    attr(x, "meansignalcoverage"), "\n")
  if (!is.null(digits)) {
    d <- getOption("digits")
    options(digits=digits)
  }
  cat("Mean signal-wise character coverage, weighted by features per signal:",
    attr(x, "weightedsignalcoverage"), "\n\nSegmentation is based on",
    attr(x, "nsignals"), "signals totalling", attr(x, "ntotalchars"),
    "characters.\nDiscounting overlaps, the segmentation above accounts for",
    attr(x, "charscovered"),
    "of those characters.\nTotal character coverage rate:",
    attr(x, "charscovered") / attr(x, "ntotalchars"), "\n")
  if (!is.null(digits))
    options(digits=d)
  invisible(x)
}

# helper functions

runfromformula <- function(FUN, x, y, ...) {
  t <- stats::terms(x)
  fields <- rownames(attr(t, "factors"))
  lhs <- fields[attr(t, "response")]
  rhs <- attr(t, "term.labels")
  # TODO replace with attach()
  FUN(y[,lhs], y[,rhs,drop=FALSE], ...)
}

# selectfrom is a list of N elements, each specifying one or more elements that
# could be selected for that slot (by default determined via length(element)).
# given a function testfun(selected-elements),
# returns either the first or a matrix of all selections of elements from
# selectfrom which satisfy testfun, with selection indices along columns.
# if mapply(selectionfun, selectfrom, is) is not supposed to simplify the
# elements to a vector or matrix, protect the result with a list()
find.selections <- function(selectfrom, testfun, first=FALSE,
    noptions=sapply(selectfrom, length), selectionfun="[[",
    include.value=FALSE) { # , incremental=FALSE
  selections <- NULL
  is <- rep(1, length(selectfrom))
  # increment range selections from back to front
  nextinc <- length(selectfrom)
  i <- 0
  while (nextinc > 0) {
    val <- testfun(mapply(selectionfun, selectfrom, is))
    if (!is.null(val) && (length(val) > 1 || val)) {
      if (first)
        return(is) # TODO include value here too?
      else
        selections <- rbind(selections, c(list(sels=is), val))
    }

    # move on to next selection combination
    while (nextinc > 0) {
      if (is[nextinc] == noptions[nextinc]) {
        # overflow: reset and carry over to next segment
        is[nextinc] <- 1
        nextinc <- nextinc - 1
      } else {
        # try next selection for current field
        is[nextinc] <- is[nextinc] + 1
        # if incremental==TRUE, at least one selection has to be on index 1
#        if (!incremental || any(is == 1)) {
          nextinc <- length(selectfrom)
          break
#        }
      }
    }
  }
  if (include.value)
    selections
  else if (!is.null(selections))
    t(sapply(selections[,1], I))
}
