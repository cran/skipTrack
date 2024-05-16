test_that("results function runs", {
  set.seed(1)
  tst <- skipTrack.simulate(10, skipProb = c(.7,.2,.1))
  foo <- skipTrack.fit(Y = tst$Y, cluster = tst$cluster, X = tst$X, Z = tst$Z,
                       chains = 2, reps = 100, useParallel = FALSE)
  res <- skipTrack.results(foo, burnIn = 10)

  expect_equal(length(res), 4)
})

#test_that('on easiest example we get at least decent coverage', {
#  set.seed(1)
#  tst <- skipTrack.simulate(50, skipProb = c(.7,.2,.1))
#  foo <- skipTrack.fit(Y = tst$Y, cluster = tst$cluster, X = tst$X, Z = tst$Z,
#                       chains = 2, reps = 100, useParallel = FALSE)
#  res <- skipTrack.results(foo, tst, burnIn = 10)#
#
#  expect_gt(mean(res$cijs$Coverage), .5)
#  expect_true(res$Betas$Coverage)
#  expect_true(res$Gammas$Coverage)
#})
