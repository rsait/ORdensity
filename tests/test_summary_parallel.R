test_that("Test summary parallel",{
  cat("Checking summary parallel not changed\n")
  load("ORdensity_parallel.Rdata")
  load("summary_parallel.Rdata")
  expect_equal(summary(ORdensity_parallel), summary_parallel)
})
