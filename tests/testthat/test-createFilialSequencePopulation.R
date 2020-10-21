test_that("createFilialSequencePopulation works for example sequences", {
  x <- createFilialSequencePopulation(c('AAABBBCCCDDDEEEFFF','GGGHHHIIIJJJKKKLLL'), 3)
  expect_is(x, "character")
})

