test_that("Test findDEgenes parallel",{
  cat("Checking finDEgenes parallel not changed\n")
  load("ORdensity_parallel.Rdata")
  load("findDEgenes_parallel.Rdata")
  expect_equal(findDEgenes(ORdensity_parallel), findDEgenes_parallel)
})