test_that("Test preclusteredData parallel",{
  cat("Checking preclusteredData parallel not changed\n")
  load("ORdensity_parallel.Rdata")
  load("preclusteredData_parallel.Rdata")
  expect_equal(preclusteredData(ORdensity_parallel), preclusteredData_parallel)
})
