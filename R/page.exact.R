#' Calculate the multiset coefficient.
#'
#' Calculates the number of possible ways to select k out of n elements with
#' replacement.
#' @seealso \url{http://en.wikipedia.org/wiki/Multiset#Counting_multisets}
#' @export
multichoose <- function(n,k)
  choose(n+k-1, k) # prod(n+0:(k-1))/factorial(k)}

#' Calculate the multinomial coefficient (number of permutations of a multiset)
#'
#' Calculates the number of possible orders in which n elements can be drawn
#' when some of the individual elements are identical. Pure R implementation
#' which is fast enough but might fail for large n/k. In this case you probably
#' want to have a look at \code{multinom} in the \code{multicool} package.
#'
#' @param n
#' @param ks vector of positive integers specifying the multiplicities of the
#'   individual elements. '1' elements can be omitted.
#' @examples
#' # only one order for the same sample drawn five times
#' mcombn(5,5)
#' # the full 5! possible orders of drawing five unique samples
#' mcombn(5,rep(1,5))
#' # sum(ks) should never exceed n, so this will give funny results:
#' mcombn(5,c(2,4))
#' @export
mcombn <- function(n,ks)
  factorial(n)/prod(sapply(ks, factorial))

# Calculates the probabilities of one replication (row) producing a certain L
# under the null hypothesis. Loops through all factorial(k) possible rankings,
# calculates their row-wise Ls and returns a table of frequencies (ordered by
# descending L)
rowwise.ls <- function(k) {
  ls <- numeric(factorial(k))
  # Heap's non-recursive algorithm cf. http://www.cs.princeton.edu/~rs/talks/perms.pdf
  comb <- 1:k
  ls[1] <- sum(comb^2)
  i <- 2
  cs <- rep(1,k)
  n <- 1
  while (n <= k) {
    if (cs[n] < n) {
      if (n%%2) {
        comb[c(1,n)] <- comb[c(n,1)] # odd recursion level
      } else {
        comb[c(cs[n],n)] <- comb[c(n,cs[n])] # even recursion level
      }
      cs[n] <- cs[n]+1
      n <- 1
      ls[i] <- sum(1:k * comb)
      i <- i+1
    } else {
      cs[n] <- 1
      n <- n+1
    }
  }
  # sanity check that Heap's worked: factorial(k) == length(which(ls>0))
  # create tabulation and divide by total count to get probabilities
  rev(table(ls))/factorial(k)
}

dropselection <- function(selection)
  vapply(which(head(selection,-1) > 0),
    function(dropindex) {
      selection[c(dropindex,dropindex+1)] <- selection[c(dropindex,dropindex+1)]+c(-1,1)
      return(selection)
    },
    numeric(length(selection)))

dropselections.apply <- function(selections)
  unique(matrix(unlist(apply(selections, 2, dropselection)), nrow=nrow(selections)), MARGIN=2)

contains.column <- function(mat, column) {
  for (i in 1:ncol(mat)) {
    if (all(column == mat[,i]))
      return(T)
  }
  return(F)
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
page.combinations.bruteforce <- function(k, N) {
  lfreq <- rowwise.ls(k)
  ls <- as.integer(names(lfreq))
  lps <- numeric(1 + N*max(ls) - N*min(ls))
  names(lps) <- (N*max(ls)):(N*min(ls))
  selections <- matrix(c(N, rep(0, length(lfreq)-1)), ncol=1)
  lps[[as.character(max(ls)*N)]] <- lfreq[[1]]^N
  while (length(selections <- dropselections.apply(selections))) {
    if (ncol(selections) == 0)
      break
    sls <- colSums(apply(selections, 2, function(r)r*ls))
    for (i in 1:ncol(selections))
      lps[[as.character(sls[i])]] <- lps[[as.character(sls[i])]] + mcombn(N,selections[which(selections[,i]>0),i])*prod(sapply(which(selections[,i]>0), function(j)prod(rep(lfreq[[j]], selections[j,i]))))
  }
  return(lps)
}

#' Calculate exact significance levels for Page's L.
#'
#' Calculate critical values of Page's L for the given significance levels.
#'
#' If only one p value is supplied, the function will return the full table of
#' p values for all Ls computed up to that p value, rather than just the
#' critical L for that particular p.
#'
#' @param k number of conditions/generations
#' @param N number of replications/chains
#' @param p.lvls vector of p levels for which the exact critical L value is sought
#' @param message.progress logical, whether to \code{message()} information on
#'   the progress of the computation to the output
#' @param L the lowest L for which the p value should be computed
#' @examples
#' pages.test.exact(6, 4, p.lvls=1)
#' @seealso page.test
#' @export
pages.test.exact <- function(k, N, p.lvls=c(0.05, 0.01, 0.001), message.progress=k*N>=30, L=0) {
  if (k==3) {
    lps <- page.combinations.bruteforce(k, N)
  } else {
    if (k+N>=11)
      message("k=", k, ", N=", N, " might take a while, I suggest you ", switch(as.character(k+N), "11"="go grab a cup of tea", "12"="go get a sandwich", "13"="go out for lunch", "14"="go on a holiday", "give up and go home"))
    if (message.progress)
      message("Calculating the distribution of the ", factorial(k), " possible row-wise L contributions..")
    lfreq <- rowwise.ls(k)
    # row-wise l frequencies are in descending order
    ls <- as.integer(names(lfreq))
    ndistinctrowls <- length(ls)
    highestL <- max(ls)*N
    if (message.progress) {
      message("The ", factorial(k), " possible rank orderings per replication for k=", k, " only yield ", ndistinctrowls, " distinct row-wise contributions to the overall L, ranging from ", min(ls), " to ", max(ls))
      message("So for k=", k, " and N=", N, " there can really only be ", max(ls)*N-min(ls)*N+1, " attested values of L, ranging from ", min(ls)*N, " to ", highestL)
      message("To cover the entire probability space exactly we'd have to look at all ", round(multichoose(length(lfreq), N)), " combinations of the possible per-row samples..")
    }
    message("Starting enumeration, terminating once we reach p level of ", max(p.lvls))
    selections <- matrix(c(N, rep(0, ndistinctrowls-1)), ncol=1)
    lps <- list()
    lps[[as.character(highestL)]] <- lfreq[[1]]^N
    p <- lps[[1]]
    downto <- max(L, ceiling(meanL(k,N)))
    for (l in (highestL-1):downto) {
      if (p >= max(p.lvls)) {
        break
      }
      selections <- dropselections.apply(selections)
      # a combination with replacement can occur in several orders (i.e. if we
      # permute the identical elements) so we have to multiply each selection
      # with its appropriate multinomial coefficient
      lps[[as.character(l)]] <- sum(apply(selections, 2, function(selection)mcombn(N,selection[which(selection>0)])*prod(sapply(which(selection>0), function(i)lfreq[[i]]^selection[i]))))
      p <- p + lps[[as.character(l)]]

      if (message.progress)
        message("L=", l, " generated by ", ncol(selections), " possible sample selections: ")
      if (message.progress)
        message("cumulative p=", p)
    }
    if (downto == ceiling(meanL(k,N))) # get remaining values for free
      lps[as.character((as.integer(tail(names(lps),1))-1):(N*min(ls)))] <- lps[(length(lps)-ifelse(meanL(k,N)%%1,0,1)):1]
  }
  if (length(p.lvls) == 1) {
    x <- rbind(L=as.numeric(names(lps)), p.ind=unlist(lps), p=cumsum(unlist(lps)))
    return(x)
  } else {
    p <- cumsum(lps)
    return(sapply(p.lvls, function(lvl) {i<-which(p>lvl)[1]-1;c(L=as.integer(names(p)[i]), p=if(i!=0)p[[i]])}))
  }
}

# calculate table of exact p values and store it in a csv file
write.page.csv <- function(k, N, p.lvl = ifelse(N+k<=10, 1.0, 0.05)) {
  if (file.exists(paste("page-k", k, "N", N, ".csv", sep="")))
    return(F)
  print(system.time(write.csv(t(pages.test.exact(k,N,p.lvl,T)), paste("page-k", k, "N", N, ".csv", sep=""), row.names=F)))
}

# load exact p values from csv files, collocate them in a big list and save it
# as an R object
save.ltable <- function(ks=3:11, Ns=2:12) {
  Ltable <- array(vector("list"), dim=c(max(Ns), max(ks)))
  for (k in ks) {
    for (N in Ns) {
      if (file.exists(paste("page-k", k, "N", N, ".csv", sep=""))) {
        data <- read.csv(paste("page-k", k, "N", N, ".csv", sep=""))
        # don't store redundant p value contributions (symmetric around mean L)
        maxindex <- match(TRUE, data$p >= 0.5)
        Ltable[N,k] <- list(data[c("L","p.ind")])
      }
    }
  }
  save(Ltable, file="Ltable.rda")
}

#' Display the range of exact p values for Page's test provided by this package.
#'
#' Displays a plot of exact p values available to \code{pages.test()} for various k, N
#'
#' @export
pages.test.exact.values <- function() {
  image(x=1:ncol(Ltable), y=1:nrow(Ltable), z=t(apply(Ltable, c(1,2), function(cell) sum(cell[[1]]$p.ind))), main="Range of exact p values available to pages.test()", xlab="k", ylab="N", col=c("green", "darkgreen"), breaks=c(0.05,0.5,1.1))
  legend("topright", legend=c("exact up to p=1", "exact up to p=.05"), fill=c("darkgreen", "green"))
}

# calculate and save the exact p values to .csv files, making your way up the
# cline of intractability. tractability 11 shouldn't take longer than a minute,
# 12 shouldn't take more than ten minutes, 13 up to 1 1/2 hours, and at 14
# memory starts to be an issue
#for (tractability in 6:14) {
#  print(tractability)
#  for (k in 4:(tractability-2)) {
#    write.page.csv(k,tractability-k)
#  }
#}
#save.ltable()

criticalLs <- function(Lcell, p.lvls=c(.05, .01, .001)) {
  unlist(lapply(p.lvls, function(p.lvl) Lcell$L[which(cumsum(Lcell$p.ind)>p.lvl)[1]-1]))
}

#' Plot the goodness of the normal approximation
#' @param p.lvls the significance levels for which the normal approximation
#'   should be compared to the exact values.
pages.test.normal.approx.goodness <- function(p.lvls=c(0.05, 0.01, 0.001)) {
  goodness <- array(NA, c(length(p.lvls), dim(Ltable)))
  for (k in 3:ncol(Ltable)) {
    for (N in 2:nrow(Ltable)) {
      if (is.null(Ltable[N,k][[1]]))
        next
      suppressWarnings(goodness[,N,k] <- ceiling(qnorm(p.lvls, lower.tail=FALSE)*sqrt(N*k^2*(k+1)*(k^2-1)/144) + (N*k*(k+1)^2)/4) - criticalLs(Ltable[N,k][[1]], p.lvls))
    }
  }
  image(x=1:ncol(Ltable), y=(1:(length(p.lvls)*nrow(Ltable)))/length(p.lvls), z=t(apply(goodness, 3, c)), xlab="k", ylab="N (.05, .01, .001)", col=c("red", "grey", "green", "darkgreen"), breaks=c(-2,-1,0,1,15))
}