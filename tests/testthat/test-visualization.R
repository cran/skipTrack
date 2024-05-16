test_that("skipTrack.visualize does not cause errors", {
  set.seed(1)
  tst <- skipTrack.simulate(10, skipProb = c(.7,.2,.1))
  foo <- skipTrack.fit(Y = tst$Y, cluster = tst$cluster, X = tst$X, Z = tst$Z,
                       chains = 2, reps = 100, useParallel = FALSE)

  expect_no_error(skipTrack.visualize(foo))
})
