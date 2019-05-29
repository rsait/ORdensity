test_that("Test new sequential",{
  cat("Checking new ORdensity sequential from simexpr not changed\n")
  myORdensity <- new("ORdensity", Exp_cond_1 = EXC.1, Exp_cond_2 = EXC.2)
  load("ORdensity_sequential.Rdata")
  expect_equal(myORdensity, ORdensity_sequential)
})
