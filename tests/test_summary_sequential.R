test_that("Test summary sequential",{
  cat("Checking summary sequential not changed\n")
  load("ORdensity_sequential.Rdata")
  load("summary_sequential.Rdata")
  expect_equal(summary(ORdensity_sequential), summary_sequential)
})
