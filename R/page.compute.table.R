# This file contains functions used to calculate the exact p values for the
# Page L statistic. This is computationally expensive for large N and k, so
# these functions are used for pre-computation

#' Display the range of exact p values for the Page test provided by this package.
#'
#' Displays a plot of exact p values available to \code{page.test()} for
#' various k, N.
#'
#' @export
page.test.exact.values <- function() {
  image(x=1:ncol(Ltable), y=1:nrow(Ltable), z=t(apply(Ltable, c(1,2),
    function(cell) sum(cell[[1]]$p.ind))),
    main="Range of exact p values available to pages.test()", xlab="k",
    ylab="N", col=c("green", "darkgreen"), breaks=c(0.05, 0.5, 1.1))
  legend("topright", legend=c("exact up to p=1", "exact up to p=.05"), fill=c("darkgreen", "green"))
}

#' Calculate the multiset coefficient.
#'
#' Calculates the number of possible ways to select k out of n elements with
#' replacement.
#' @param n
#' @param k
#' @seealso \url{http://en.wikipedia.org/wiki/Multiset#Counting_multisets}
#' @export
multichoose <- function(n, k)
  choose(n+k-1, k) # prod(n+0:(k-1))/factorial(k)}

#' Calculate the multinomial coefficient (number of permutations of a multiset)
#'
#' Calculates the number of possible orders in which n elements can be drawn
#' when some of the individual elements are identical. Only uses R core
#' functions so should be fast but might overflow for large n/k. In such cases,
#' have a look at \code{multinom} in the \code{multicool} or \code{multichoose}
#' in the \code{iterpc} package.
#'
#' @param n
#' @param ks vector of positive integers specifying the multiplicities of the
#'   individual elements. \code{1} elements can be omitted.
#' @examples
#' # only one order for the same sample drawn five times
#' mcombn(5, 5)
#' # the full 5! possible orders of drawing five unique samples
#' mcombn(5, rep(1, 5))
#' @export
mcombn <- function(n, ks) {
#  if (!all.equal(c(n, ks), as.integer(c(n, ks))) || sum(ks) > n)
#    stop("All arguments must be integer with sum(ks) == n")
  factorial(n) / prod(sapply(ks, factorial))
}

# Calculates the probabilities of one replication (row) producing a certain L
# under the null hypothesis. Loops through all factorial(k) possible rankings,
# calculates their row-wise Ls and returns a table of frequencies (ordered by
# descending L)
rowwise.ls <- function(k) {rowwise.ls
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
#  i <- 2
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
#      i <- i+1
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

dropselection <- function(selection)
  vapply(which(head(selection, -1) > 0),
    function(dropindex) {
      selection[c(dropindex, dropindex+1)] <- selection[c(dropindex, dropindex+1)]+c(-1,1)
      return(selection)
    },
    numeric(length(selection)))

ramused <- NULL
dropselections.apply <- function(selections) {
  x <- matrix(unlist(apply(selections, 2, dropselection)), nrow=nrow(selections))
  ramused <<- max(ramused, object.size(selections) + object.size(x))
  unique(x, MARGIN=2)
}

contains.column <- function(mat, column) {
  for (i in 1:ncol(mat)) {
    if (all(column == mat[,i]))
      return(TRUE)
  }
  return(FALSE)
}

# about ten times slower than the apply variant but possibly less memory-hungry
dropselections.loop <- function(selections) {
  newselections <- matrix(0, nrow=nrow(selections), ncol=ncol(selections))
  s <- 0
  for (i in 1:ncol(selections)) {
    dropped <- dropselection(selections[,i])
    for (j in 1:ncol(dropped)) {
      if (s > 0 && contains.column(as.matrix(newselections[,1:s]), dropped[,j]))
        next
      s <- s+1
      if (s > ncol(newselections))
        newselections <- cbind(newselections, dropped[,j])
      else
        newselections[,s] <- dropped[,j]
    }
  }
  return(as.matrix(newselections[,1:s]))
}

# non-clever bruteforce enumeration - slow but necessary to compute values for
# k=3 because it doesn't produce all row-wise Ls (10, 11, 13 and 14 but not 12)
page.combinations.bruteforce <- function(k, N, lfreq=rowwise.ls(k)) {
  ls <- as.integer(names(lfreq))
  lps <- numeric(1 + N*max(ls) - N*min(ls))
  names(lps) <- (N*max(ls)):(N*min(ls))
  selections <- matrix(c(N, rep(0, length(lfreq)-1)), ncol=1)
  lps[[as.character(max(ls)*N)]] <- lfreq[[1]]^N
  while (length(selections <- dropselections.apply(selections))) {
    if (ncol(selections) == 0)
      break
    sls <- colSums(apply(selections, 2, function(r) r*ls))
    for (i in 1:ncol(selections))
      lps[[as.character(sls[i])]] <- lps[[as.character(sls[i])]] +
        mcombn(N, selections[which(selections[,i] > 0), i]) *
        prod(sapply(which(selections[,i]>0), function(j) prod(rep(lfreq[[j]], selections[j,i]))))
  }
  return(lps)
}

# helper function for incremental computation status writeout
append.csv.row <- function(filename, vector)
  write.table(t(vector), file=filename, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)

#' Calculate exact significance levels of the Page \emph{L} statistic.
#'
#' Calculate critical values of \emph{L} for the given significance levels.
#'
#' If only one p value is supplied, the function will return the full table of
#' p values for all \emph{L} up to that p value, rather than just the critical
#' \emph{L} for that particular p.
#'
#' @param k number of conditions/generations
#' @param N number of replications/chains
#' @param p.lvls vector of p levels for which the exact critical L value is sought
#' @param message.progress logical, whether to \code{message()} information on
#'   the progress of the computation to the output
#' @param L the lowest L for which the p value should be computed
#' @examples
#' page.test.exact(6, 4, p.lvls=1)
#' @seealso \code{\link{page.test}}
page.compute.exact <- function(k, N, p.lvls=c(0.05, 0.01, 0.001), message.progress=(k*N >= 20), L=0, csvfile=NULL, lfreq=rowwise.ls(k)) {
  if (k == 3) {
    lps <- page.combinations.bruteforce(k, N, lfreq=lfreq)
  } else {
    if (file.exists(csvfile))
      stop("Output file ", csvfile, " already exists, won't overwrite.")
    if (message.progress) {
      message("==========")
      message("k=", k, ", N=", N)
      message("==========")
    }
    # row-wise l frequencies are in descending order
    ls <- as.integer(names(lfreq))
    ndistinctrowls <- length(ls)
    highestL <- max(ls)*N
    if (message.progress) {
      message("The ", factorial(k), " possible rank orderings per replication for k=",
        k, " only yield ", ndistinctrowls, " distinct row-wise contributions to the overall L, ranging from ",
        min(ls), " to ", max(ls), ", so there can really only be ",
        max(ls)*N-min(ls)*N+1, " attested values of L, ranging from ", min(ls)*N, " to ", highestL)
      message("To cover the entire probability space exactly we'd have to look at all ",
        round(multichoose(length(lfreq), N)), " combinations of the ", max(ls)*N-min(ls)*N+1, " possible per-row samples..")
    }
    message("Starting enumeration, terminating once we reach p level of ", max(p.lvls))
    selections <- matrix(c(N, rep(0, ndistinctrowls-1)), ncol=1)

    # TODO cachedmcombn <- memoise::memoise(mcombn)
    lps <- list()
    lps[[as.character(highestL)]] <- lfreq[[1]]^N
    p <- lps[[1]]
    downto <- max(L, ceiling(pages.L.mean(k,N)))
    for (l in (highestL-1):downto) {
      if (p >= max(p.lvls)) {
        break
      }
      proc <- proc.time()
      selections <- dropselections.apply(selections)
      # a combination with replacement can occur in several orders (i.e. if we
      # permute the identical elements) so we have to multiply each selection
      # with its appropriate multinomial coefficient
      lps[[as.character(l)]] <- sum(apply(selections, 2,
        function(selection) {
          mcombn(N, selection[which(selection > 0)]) *
            prod(sapply(which(selection > 0), function(i) lfreq[[i]]^selection[i]))
        }))
      p <- p + lps[[as.character(l)]]
      if (!is.null(csvfile))
        # write: L, p.ind, p, MB RAM used, proc.time() (user, system, real)
        append.csv.row(csvfile, c(l, lps[[as.character(l)]], p, ramused/1048576, (proc.time()-proc)[1:3]))
      if (message.progress)
        message("k=", k, ", N=", N, ", L=", l, ", cumulative p=", round(p, 6), " (", ncol(selections), " possible selections, ", round(ramused/1048576, 2), "MB)")
      ramused <<- NULL
    }
    if (message.progress)
      message("Done!")
    if (downto == ceiling(pages.L.mean(k, N))) # get remaining values for free
      lps[as.character((as.integer(tail(names(lps),1))-1):(N*min(ls)))] <- lps[(length(lps)-ifelse(pages.L.mean(k, N)%%1,0,1)):1]
  }
  if (length(p.lvls) == 1) {
    x <- rbind(L=as.numeric(names(lps)), p.ind=unlist(lps), p=cumsum(unlist(lps)))
    return(x)
  } else {
    p <- cumsum(lps)
    return(sapply(p.lvls,
      function (lvl) {
        i <- which(p > lvl)[1] - 1
        c(L=as.integer(names(p)[i]), p=if(i != 0) p[[i]])
      }))
  }
}

# faster method by incremental convolution of row-wise distributions
page.compute.exact <- function(k, maxN, p.lvl=1, lfreq=rowwise.ls(k)) {
  ls <- as.integer(names(lfreq))

  currentls <- ls
  convolutions <- list(lfreq)

  N <- 1
  while (N < maxN) {
    currentls <- (currentls[1] + ls[1]) : (currentls[length(currentls)] + ls[length(ls)])
    # calculate all joint probabilities
    jointp <- outer(convolutions[[N]], lfreq)
    # pad matrix so that equal sums of L are aligned in rows
    # (resulting matrix will have nrow(jointp)+ncol(jointp)-1 rows)
    jointp <- sapply(1:ncol(jointp), function (col) {
      c(rep(0, col-1), jointp[,col], rep(0, max(0, ncol(jointp)-col)))
    })
    # it's that easy
    newconvolution <- rowSums(jointp)
    names(newconvolution) <- currentls

    N <- N+1
    convolutions[[N]] <- newconvolution
  }
  return(convolutions)
}

page.csv.filename <- function(k, N)
  paste("page-k", k, "N", N, ".csv", sep="")

# calculate table of exact p values and store it in a csv file
write.page.csv <- function(k, N, p.lvl=1.0) {
  csvfile <- page.csv.filename(k, N)
  lfreqfile <- paste("lfreq-k", k, ".rda", sep="")
  if (file.exists(lfreqfile)) {
    load(lfreqfile)
  } else {
    message("Calculating null distribution of the ", factorial(args[1]), " possible row-wise Ls..")
    pretime <- system.time(lfreq <- rowwise.ls(args[1]))
    save(lfreq, file=lfreqfile)
    append.csv.row("page-resources.csv", c(k, NA, object.size(lfreq)/1048576, pretime[1:3]))
  }

  totaltime <- system.time(result <- t(page.compute.exact(k, N, p.lvl, csvfile=csvfile, lfreq=lfreq)))
  if (k != 3) {
    while (!file.exists(csvfile))
     Sys.sleep(1)
    stats <- read.csv(csvfile, col.names=c("L", "p.ind", "p", "RAM", "user", "system", "real"))
    # write resource usage summary
    append.csv.row("page-resources.csv", c(k, N, max(stats$RAM), totaltime[1:3]))
    message("Complete enumeration took ", paste(totaltime[1:3], collapse=" "),
      ", peak RAM usage was ", round(max(stats$RAM), 2), "MB")
    # rewrite bare results without resource stats
  }
  write.csv(result, csvfile, row.names=FALSE)
}

if ("--file=page.compute.table.R" %in% commandArgs()) {
  args <- as.numeric(commandArgs(trailingOnly=TRUE))
  if (is.na(args[1]))
    stop("Usage: Rscript page.compute.table.R <k> [[minN] maxN]")
  if (!is.na(args[3])) {
    Ns <- args[2]:args[3]
  } else if (!is.na(args[2])) {
    Ns <- 2:args[2]
  } else {
    Ns <- 2:12
  }
  source("page.R")
  for (N in Ns)
    write.page.csv(args[1], N)
}

plot.resource.usage <- function() {
  library(lattice)
  data <- read.csv("page-resources.csv", header=FALSE, col.names=c("k", "N", "RAM", "user", "system", "real"))
  compN <- xyplot(user ~ N, data=data, group=k, type=c("l", "p"),
    auto.key=list(title="k", space="right", reverse.rows=TRUE))
  memoryN <- xyplot(RAM ~ N, data=data, group=k, type=c("l", "p"))

  # add baseline
  data[is.na(data$N), 'N'] <- 0
  compk <- xyplot(user ~ k, data=data, group=N, type=c("l", "p"),
    auto.key=list(title="N", space="right", reverse.rows=TRUE))
  memoryk <- xyplot(RAM ~ k, data=data, group=N, type=c("l", "p"))
  gridExtra::grid.arrange(compN, memoryN, compk, memoryk, ncol=2, nrow=2)
}

# once you've run write.page.csv() for a number of k and N, use this function
# to load the pre-computed p values from the csv files, collocate them in a big
# list and save the resulting lookup table as an R object
save.ltable <- function(ks=3:14, Ns=2:20) {
  Ltable <- array(vector("list"), dim=c(max(Ns), max(ks)))
  for (k in ks) {
    for (N in Ns) {
      if (file.exists(page.csv.filename(k, N))) {
        data <- read.csv(page.csv.filename(k, N))
        # don't store redundant p value contributions (symmetric around mean L)
#        maxindex <- match(TRUE, data$p >= 0.5)
        Ltable[N,k] <- list(data[c("L", "p.ind")])
      }
    }
  }
  print(Ltable)
  save(Ltable, file="Ltable.rda")
}
