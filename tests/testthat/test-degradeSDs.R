
test_that("degradeSDs produces character output", {
  sdMaximalHBS <- 10
  increaseHZEI <- TRUE
  # Initiaing the Codons matrix plus corresponding amino acids
  ntSequence <- 'TTTTCGATCGGGATTAGCCTCCAGGTAAGTATCTATCGATCTATGCGATAG'
  # Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(ntSequence, '')[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  fan <- matrix(nrow=2, ncol=nchar(ntSequence)/3)
  fan[1, ] <- sst
  fan[2, ] <- sst
  x <- degradeSDs(fan, maxhbs=sdMaximalHBS, increaseHZEI=increaseHZEI)
  expect_is(x, "character")
})



test_that("degradeSDs keeps AA sequence intact", {
  sdMaximalHBS <- 10
  increaseHZEI <- TRUE
  # Initiaing the Codons matrix plus corresponding amino acids
  ntSequence <- 'TTTTCGATCGGGATTAGCCTCCAGGTAAGTATCTATCGATCTATGCGATAG'
  # Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(ntSequence, '')[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  fan <- matrix(nrow=2, ncol=nchar(ntSequence)/3)
  fan[1, ] <- sst
  fan[2, ] <- sst
  x <- degradeSDs(fan, maxhbs=sdMaximalHBS, increaseHZEI=increaseHZEI)

    ## Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(x, "")[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  codonSeqFinal <- Codons$AA[match(sst, Codons$seq)]

  ## Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(ntSequence, "")[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  codonSeqCDS <- Codons$AA[match(sst, Codons$seq)]

  expect_equal(codonSeqFinal, codonSeqCDS)

})





test_that("degradeSDs successfully degrades HBS", {
  sdMaximalHBS <- 10
  increaseHZEI <- TRUE
  # Initiaing the Codons matrix plus corresponding amino acids
  ntSequence <- 'TTTTCGATCGGGATTAGCCTCCAGGTAAGTATCTATCGATCTATGCGATAG'
  # Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(ntSequence, '')[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  fan <- matrix(nrow=2, ncol=nchar(ntSequence)/3)
  fan[1, ] <- sst
  fan[2, ] <- sst
  x <- degradeSDs(fan, maxhbs=sdMaximalHBS, increaseHZEI=increaseHZEI)
  xseq <- substr(x, 22, 32)
  xseq <- hbg$hbs[hbg$seq == xseq]

  xseqold <- substr(ntSequence, 22, 32)
  xseqold <- hbg$hbs[hbg$seq == xseqold]

  test <- xseqold > xseq
  expect_equal(test, TRUE)
})
