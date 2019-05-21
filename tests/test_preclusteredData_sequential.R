test_that("Test preclusteredData sequential",{
  cat("Checking preclusteredData sequential not changed\n")
  load("ORdensity_sequential.Rdata")
  load("preclusteredData_sequential.Rdata")
  expect_equal(preclusteredData(ORdensity_sequential), preclusteredData_sequential)
})
