test_that("Test findDEgenes sequential",{
  cat("Checking finDEgenes sequential not changed\n")
  load("ORdensity_sequential.Rdata")
  load("findDEgenes_sequential.Rdata")
  expect_equal(findDEgenes(ORdensity_sequential), findDEgenes_sequential)
})
