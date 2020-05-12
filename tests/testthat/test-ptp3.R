test_that("ptp3 internal log consistency", {
  expect_equal(log(ptp3(0.1, 0, 1 ,0, param = "eps", FUN = plogis, log.p = F)),
               ptp3(0.1, 0, 1, 0, param = "eps", FUN = plogis, log.p = T))
})

test_that("ptp3 conistency log with plogis", {
  expect_equal(plogis(0.1, location = 0, scale = 1, log.p = T),
               ptp3(0.1, 0, 1, 0, param = "eps", FUN = plogis, log.p = T))
})

test_that("ptp3 conistency with function", {
  expect_equal(log(plogis(0.1, location = 0, scale = 1)),
               ptp3(0.1, 0, 1, 0, param = "eps", FUN = plogis, log.p = T))
})

