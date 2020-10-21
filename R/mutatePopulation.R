## Introduce mutations to sequenceVector
mutatePopulation <- function(sequenceVector,
                             codonReplacementChance) {

  ## Check if codonsInsert table is present
  if (!any(ls() == "codonsInsert"))
    codonsInsert <- Codons

  ## For every parent
  for (i in seq_len(length(sequenceVector))) {

    ## Split sequence into codons
    coSeq <- strsplit(sequenceVector[i], "")[[1]]
    coSeq <- paste0(coSeq[c(TRUE, FALSE, FALSE)],
                    coSeq[c(FALSE, TRUE, FALSE)],
                    coSeq[c(FALSE, FALSE, TRUE)])

    ## For every Codon
    for (changePos in seq_len(length(coSeq))) {

      ## By codonReplacementChance exchange codon
      aa <- Codons$AA[Codons$seq == coSeq[changePos]]
      aa <- codonsInsert$seq[codonsInsert$AA == aa]
      aaProb <- aa
      aaProb[aaProb != coSeq[changePos]] <- codonReplacementChance
      aaProb[aaProb == coSeq[changePos]] <- (1 - codonReplacementChance)

      coSeq[changePos] <- sample(aa, 1, prob = aaProb)


    }
    sequenceVector[i] <- paste(coSeq, collapse = "")
  }

  ## Return results
  return(sequenceVector)
}
