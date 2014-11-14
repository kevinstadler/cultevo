# meaning spaces and utility functions for some specific experimental designs

#' Construct the meaning space matrix for the 'optional arguments' experiment.
#'
#' Constructs a matrix with all meaning combinations for the optional arguments
#' experiment. the globally shared part of this meaning space design is 2x2x2x2,
#' but then the upper half of the meaning space this is again combined with
#' all combinations of an additional 2x2 space (the ground) while the
#' corresponding lower half is filled with dummy values for the non-specified
#' ground according to the unspecifieddistance argument which can be 1, 2 or 3.
#' @examples
#' meaningspace.optionalarguments(3)
#' @export
meaningspace.optionalarguments <- function(unspecifieddistance) {
  objects <- c("ci", "sq")
  number <- c("sg", "pl")
  actions <- c("bn", "sl")
  perf <- c("pf", "cn")

  # 64 transitive meanings
  transitives <- allmeaningcombinations(list(objects, objects, number, number, actions, perf))

  # the three ways of encoding the optional arguments in the intransitive
  if (unspecifieddistance == 2) {
    # encode the dummy arguments as 0s but without specifying transitivity
    # with its own flag, leading to a distance of 2.
    rbind(transitives, allmeaningcombinations(list(objects, "00", number, "00", actions, perf)))
  } else {
    # prepend explicit transitivity flag
    transitives <- cbind("tr", transitives)
    if (unspecifieddistance == 1) {
      # specify the unrealised arguments as NA which leads to a
      # distance of 1 to the corresponding transitive sentence
      unmarked <- NA
    } else {
      # specifying the dummy arguments as 0 leads to a distance of 3
      unmarked <- "00"
    }
    intransitives <- cbind("in", allmeaningcombinations(list(objects, unmarked, number, unmarked, actions, perf)))
    rbind(transitives, intransitives)
  }
}

# Read optional arguments experiment data and plot Mantel tests.
# 
# Runs and visualises Mantel tests for the initial language and 8 generations
# of the optional arguments experiment, with data read from the given file.
#' @export
experiment.optionalarguments <- function(filename, unspecifieddistance, ...) {
  strings <- read.table(filename, sep="\t")[,1:9]
  strings <- apply(strings, c(1,2), toString)
  meaningspace <- meaningspace.optionalarguments(unspecifieddistance)
  mantel.development(meaningspace, t(strings), ...)
}

# Run and plot Mantel tests for all generations read from the local files
# "chain1.tsv" to "chain4.tsv"
#' @export
experiment.optionalarguments.allchains <- function(unspecifieddistance=3, plot="r", ...) {
  # plot all 4 chains next to each other
  if (!is.null(plot)) {
    par(mfrow=c(1,4))
  }
  for (i in 1:4) {
    print(paste("Chain", i))
    filename <- paste("chain", i, ".tsv", sep="")
    print(experiment.optionalarguments(filename, unspecifieddistance, plot=plot, ylim=c(-0.1, 0.8), ...))
  }
}

# visualise the two different measures
#experiment.optionalarguments.allchains(1, plot="r", ylim=c(-0.1, 0.7))
#experiment.optionalarguments.allchains(1, plot="z", ylim=c(-1, 30))
# original plot with the old distance measure (intransitives have distance 3)
#experiment.optionalarguments.allchains(3, plot="z", ylim=c(-1, 30))
# look at the underlying r and mantel test variance
#experiment.optionalarguments.allchains(3, plot="r", ylim=c(-0.1, 0.7))

