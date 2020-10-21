
test_that("slidingWindowHZEImanipulation outputs character value", {
  inSeq <- 'ATGGAAGACGCCAAAAACATAAAGAAAGGCAGGCTAAGCCTAGCTTGCCATTGCCCGGCGCCATTCTATCCGCTGGAAGATGGAATT'
  x <- slidingWindowHZEImanipulation(inSeq, increaseHZEI=TRUE)
  y <- slidingWindowHZEImanipulation(inSeq, increaseHZEI=FALSE)

  x <- is.character(x)
  y <- is.character(y)

  x <- (x == TRUE) & (y == TRUE)

  expect_equal(x, TRUE)
})


test_that("slidingWindowHZEImanipulation works for example sequence", {
  inSeq <- 'ATGGAAGACGCCAAAAACATAAAGAAAGGCAGGCTAAGCCTAGCTTGCCATTGCCCGGCGCCATTCTATCCGCTGGAAGATGGAATT'
  x <- as.character(slidingWindowHZEImanipulation(inSeq, increaseHZEI=TRUE))
  y <- as.character(slidingWindowHZEImanipulation(inSeq, increaseHZEI=FALSE))

  ## Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(x, "")[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  codonSeqFinal <- Codons$AA[match(sst, Codons$seq)]

  ## Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(inSeq, "")[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)], sst[c(FALSE, FALSE, TRUE)])
  codonSeqCDS <- Codons$AA[match(sst, Codons$seq)]

  expect_equal(codonSeqFinal,codonSeqCDS)
})
