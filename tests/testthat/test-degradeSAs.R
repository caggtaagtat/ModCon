
test_that("degradeSAs produces character output", {
  sdMaximalHBS <- 10
  acMaximalMaxent <- 4
  increaseHZEI <- TRUE
  # Initiaing the Codons matrix plus corresponding amino acids
  ntSequence <- 'TTTTGTCTTTTTCTGTGTGGCAGTGGGATTAGCCTCCTATCGATCTATGCGATA'
  # Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(ntSequence, '')[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  fan <- matrix(nrow=2, ncol=nchar(ntSequence)/3)
  fan[1, ] <- sst
  fan[2, ] <- sst
  x <- degradeSAs(fan, maxhbs=sdMaximalHBS, maxME=acMaximalMaxent, increaseHZEI=increaseHZEI)
  expect_is(x, "character")
})



test_that("degradeSAs keeps AA sequence intact", {
  sdMaximalHBS <- 10
  acMaximalMaxent <- 4
  increaseHZEI <- TRUE
  # Initiaing the Codons matrix plus corresponding amino acids
  ntSequence <- 'TTTTGTCTTTTTCTGTGTGGCAGTGGGATTAGCCTCCTATCGATCTATGCGATA'
  # Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(ntSequence, '')[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  fan <- matrix(nrow=2, ncol=nchar(ntSequence)/3)
  fan[1, ] <- sst
  fan[2, ] <- sst
  x <- degradeSAs(fan, maxhbs=sdMaximalHBS, maxME=acMaximalMaxent, increaseHZEI=increaseHZEI)

    ## Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(x, "")[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  codonSeqFinal <- Codons$AA[match(sst, Codons$seq)]

  ## Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(ntSequence, "")[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  codonSeqCDS <- Codons$AA[match(sst, Codons$seq)]

  expect_equal(codonSeqFinal,codonSeqCDS)

})





test_that("degradeSDs successfully degrades HBS", {
  sdMaximalHBS <- 10
  acMaximalMaxent <- 4
  increaseHZEI <- TRUE
  # Initiaing the Codons matrix plus corresponding amino acids
  ntSequence <- 'TTTTGTCTTTTTCTGTGTGGCAGTGGGATTAGCCTCCTATCGATCTATGCGATA'
  # Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(ntSequence, '')[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  fan <- matrix(nrow=2, ncol=nchar(ntSequence)/3)
  fan[1, ] <- sst
  fan[2, ] <- sst
  x <- degradeSAs(fan, maxhbs=sdMaximalHBS, maxME=acMaximalMaxent, increaseHZEI=increaseHZEI)
  xseq <- substr(x, 4, 26)
  xseq <- calculateMaxEntScanScore(xseq, 3)

  xseqold <- substr(ntSequence, 4, 26)
  xseqold <- calculateMaxEntScanScore(xseqold, 3)

  test <- xseqold > xseq
  expect_equal(test, TRUE)
})
