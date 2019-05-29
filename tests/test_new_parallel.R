test_that("Test new parallel",{
  cat("Checking new ORdensity parallel from simexpr not changed\n")
  myORdensity <- new("ORdensity", Exp_cond_1 = EXC.1, Exp_cond_2 = EXC.2, parallel = TRUE)
  load("ORdensity_parallel.Rdata")
  expect_equal(myORdensity, ORdensity_parallel)
})
