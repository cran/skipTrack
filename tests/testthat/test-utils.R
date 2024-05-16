test_that("Utility functions do not throw errors", {
  set.seed(1)
  tst <- skipTrack.simulate(10, skipProb = c(.7,.2,.1))
  foo <- skipTrack.fit(Y = tst$Y, cluster = tst$cluster, X = tst$X, Z = tst$Z,
                       chains = 2, reps = 200, useParallel = FALSE)

  expect_null(plot(foo))
  expect_null(print(foo))
  expect_no_error(summary(foo, burnIn = 20))
  expect_null(str(foo))
})
