test_that('All simulation methods output necessary components', {
  st.tst <- skipTrack.simulate(5, model = 'skipTrack', skipProb = c(.5, .5))
  li.tst <- skipTrack.simulate(5, model = 'li', maxCycles = 2)
  mix.tst <- skipTrack.simulate(5, model = 'mixture', skipProb = c(.5, .5))

  nms <- c('Y', 'cluster', 'X', 'Z', 'Beta', 'Gamma', 'NumTrue', 'Underlying')

  expect_equal(names(st.tst), nms)
  expect_equal(names(li.tst), nms)
  expect_equal(names(mix.tst), nms)
})
