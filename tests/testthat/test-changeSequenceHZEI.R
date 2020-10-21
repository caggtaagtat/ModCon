test_that("changeSequenceHZEI is able to increase HZEI", {

  inSeq <- 'ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGAACC'
  res <- changeSequenceHZEI(inSeq, increaseHZEI=TRUE)
  x <- calculateHZEIint(res[[3]])
  xold <- calculateHZEIint(inSeq)
  check <- x>xold

  expect_equal(check, TRUE)
})


test_that("changeSequenceHZEI is able to decrease HZEI", {

  inSeq <- 'ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGAACC'
  res <- changeSequenceHZEI(inSeq, increaseHZEI=FALSE)
  x <- calculateHZEIint(res[[3]])
  xold <- calculateHZEIint(inSeq)
  check <- x<xold

  expect_equal(check, TRUE)
})



test_that("changeSequenceHZEI is able to decrease HZEI using the genetic algorithm", {

  inSeq <- 'ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGAACC'
  res <- changeSequenceHZEI(inSeq, increaseHZEI=FALSE, optiRate=50, nCores=1)
  x <- calculateHZEIint(res[[3]])
  xold <- calculateHZEIint(inSeq)
  check <- x<xold

  expect_equal(check, TRUE)
})



test_that("changeSequenceHZEI is able to increase HZEI using the genetic algorithm", {

  inSeq <- 'ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGAACC'
  res <- changeSequenceHZEI(inSeq, increaseHZEI=TRUE, optiRate=50, nCores=1)
  x <- calculateHZEIint(res[[3]])
  xold <- calculateHZEIint(inSeq)
  check <- x>xold

  expect_equal(check, TRUE)
})

