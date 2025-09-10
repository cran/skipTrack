test_that("skipTrack.diagnostics runs on all options", {
  set.seed(1)
  tst <- skipTrack.simulate(20, skipProb = c(.9, .1))
  foo <- skipTrack.fit(Y = tst$Y, cluster = tst$cluster, X = tst$X, Z = tst$Z,
                       chains = 2, reps = 10, useParallel = FALSE)

  expect_no_error(skipTrack.diagnostics(foo, 'rho'))
  expect_no_error(skipTrack.diagnostics(foo, 'phi'))
  expect_no_error(skipTrack.diagnostics(foo, 'Betas'))
  expect_no_error(skipTrack.diagnostics(foo, 'Gammas'))
  expect_no_error(skipTrack.diagnostics(foo, 'muijs'))
  expect_no_error(skipTrack.diagnostics(foo, 'tauis'))
  expect_no_error(skipTrack.diagnostics(foo, 'cijs'))
})
