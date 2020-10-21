## Create children in different ways
recombineTwoSequences <- function(ntSequence1,
                                  ntSequence2,
                                  preferenceVector) {

  MI1 <- strsplit(ntSequence1, "")[[1]]
  MI1 <- paste0(MI1[c(TRUE, FALSE, FALSE)],
                MI1[c(FALSE, TRUE, FALSE)],
                MI1[c(FALSE, FALSE, TRUE)])

  MI2 <- strsplit(ntSequence2, "")[[1]]
  MI2 <- paste0(MI2[c(TRUE, FALSE, FALSE)],
                MI2[c(FALSE, TRUE, FALSE)],
                MI2[c(FALSE, FALSE, TRUE)])

  mode <- sample(c(1, 2, 3), prob = preferenceVector, size = 1)

  ## Mode 1 is random from parent one and two
  if (mode == 1) {

    child = ""
    for (i in seq_len(length(MI1))) {

      if (sample(c(0, 1), size = 1) == 0) {
        child <- paste0(child, MI1[[i]])
      } else {
        child <- paste0(child, MI2[[i]])
      }

    }
  }

  ## Mode 2 similating crossover
  if (mode == 2) {

    child = ""

    # generate crossover point
    crossoverpoint <- sample(c(2:(length(MI1) - 1)), 1)
    child <- paste(MI1[seq_len(crossoverpoint)],
                   collapse = "")
    child <- paste0(child, (paste(MI2[c((crossoverpoint + 1):
                                          length(MI1))],
                                  collapse = "")))

  }

  ## Mode 3 is random from parent one and two
  if (mode == 3) {

    child = MI1

    ## generate crossover point
    crossoverpointSeq <- sample(c(2:(length(MI1) - 1)), 1)
    crossoverpointSeq2 <- sample(c(2:(length(MI1) - 1)), 1)
    min <- min(crossoverpointSeq, crossoverpointSeq2)
    max <- max(crossoverpointSeq, crossoverpointSeq2)

    child[c(min:max)] <- MI2[c(min:max)]
    child <- paste(child, collapse = "")

  }

  child

}




