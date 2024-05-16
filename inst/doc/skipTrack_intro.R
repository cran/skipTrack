## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

set.seed(1)

## -----------------------------------------------------------------------------
library(skipTrack)

## -----------------------------------------------------------------------------
#Simulate data
dat <- skipTrack.simulate(n = 100, model = 'skipTrack', skipProb = c(.75, .2, .05))

names(dat)

## ----fig.align='center', fig.width = 7----------------------------------------
#Histogram of observed outcomes
hist(dat$Y, breaks = 10:150)

## -----------------------------------------------------------------------------
ft <- skipTrack.fit(Y = dat$Y, cluster = dat$cluster,
                    reps = 1000, chains = 4, useParallel = FALSE)

## ----fig.align='center', fig.width = 7, fig.height=7--------------------------
skipTrack.diagnostics(ft, param = 'cijs')

## ----fig.align='center', fig.width = 7, fig.height=7--------------------------
plot(ft)

## -----------------------------------------------------------------------------
summary(ft)

summary(ft, burnIn = 500)

