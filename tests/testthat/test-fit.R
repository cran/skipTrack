test_that("skipTrack.fit runs", {
  tst <- skipTrack.simulate(5, skipProb = c(.5,.5), trueBetas = c(.01, 0),
                               trueGammas = c(.2, 0))
  foo <- skipTrack.fit(Y = tst$Y, cluster = tst$cluster, X = tst$X, Z = tst$Z,
                       chains = 2, reps = 10, useParallel = FALSE)

  expect_equal(names(foo), c('fit', 'data', 'useParallel', 'model'))
  expect_equal(length(foo$fit), 2)
  expect_equal(length(foo$fit[[1]]), 11)
})
