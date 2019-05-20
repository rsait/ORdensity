test_that("Test simexpr",{
  cat("Checking simexpr not changed\n")
  load("simexpr_original.Rdata")
  expect_equal(simexpr, simexpr_original)
})
