#This document holds all of the functions that provide random draws from full conditional
#posteriors for the skipTrack algorithm.

#' Draw from Posterior Distribution for Beta Parameters
#'
#' In our model, \code{mui} follows a normal distribution with mean \eqn{X_{i}^{T} \beta + b_{i}} and precision \eqn{tau_i}.
#' Additionally, we assume that \code{beta} follows a multivariate normal prior with mean 0 and precision \eqn{rhoBeta \cdot I}.
#' This function draws from the posterior distribution of \code{beta} under these assumptions.
#'
#' @param rhoBeta A scalar representing the prior precision parameter for beta.
#' @param X A matrix of covariates, where each row represents a cycle and each column represents a covariate.
#' @param b A vector where each element is the random effect/intercept for individual \code{i}.
#' @param m A vector of observed means (\code{mui}) for each cycle.
#' @param tau A vector where each element is the precision for individual \code{i} (length = number of individuals).
#' @param indFirst An logical vector (length = number of individuals); each entry is TRUE if this is the first cycle for that individual in the vector of observations. Used to identify submatrices of X and m.
#'
#' @return A numeric vector representing a draw from the posterior distribution of beta parameters.
#'
#' @details
#' For each individual, the function extracts the relevant rows of \code{X} and \code{m} using \code{indFirst}, and multiplies by the individual's precision \code{tau[i]}. It then computes the updated posterior precision and mean for \code{beta} and returns a sample from the resulting multivariate normal distribution. Requires the \code{mvtnorm} package.
#'
postBeta <- function(rhoBeta = .01, X, b, m, tau, indFirst){
  #Assuming X is (num cycles)x(dimension of beta) matrix of covariates
  ind <- cumsum(indFirst) #Variables to identify individuals for sub matrices

  #Extract Xi matrices from X
  Xi <- lapply(1:max(ind), function(i){
    Xi <- X[ind == i,,drop = F]
  })

  #Build XTX matrices
  XTXs <- lapply(1:max(ind), function(i){
    tau[i]*t(Xi[[i]]) %*% Xi[[i]]
  })

  #Build other sum in mean
  meanSums <- lapply(1:max(ind), function(i){
    tau[i] * t(Xi[[i]]) %*% (m[ind == i] - rep(b[i], sum(ind == i)))
  })

  postPre <- (rhoBeta)*diag(1, dim(XTXs[[1]])[1]) +  Reduce('+', XTXs)
  postVar <- solve(postPre)
  postMean <- postVar %*% Reduce('+', meanSums)

  return(mvtnorm::rmvnorm(1, mean = postMean, postVar))
}

#' Perform a Metropolis-Hastings Step for Drawing a New Gamma
#'
#' Our model assumes that tau_i ~ Gamma(alpha_i, phi) where alpha_i/phi = theta_i and
#' g(theta_i) = Z_i^T Gamma, where g is the log-link function.
#' Because of this GLM formulation we cannot simply draw from a posterior here and instead use a
#' Metropolis-Hastings step with proposal distribution propGamma ~ MVNormal(currentGamma, rhoGamma*I)
#'
#' @param taui A vector of length num_Individuals representing precision parameters for individuals.
#' @param Z A matrix of covariates, where each row represents an individual and each column represents a covariate.
#' @param currentGamma A vector (or matrix with 1 row) representing the current Gamma value.
#' @param phi A scalar, the prior rate for tau_i.
#' @param rhoGamma A scalar representing the proposal distribution precision parameter.
#'
#' @details This draw step assumes a log-link function for the Gamma GLM that we are fitting.
#'
#' @return A list containing the new Gamma value and the corresponding thetai values.
#'
postGamma <- function(taui, Z, currentGamma, phi = 1, rhoGamma = 1000){
  #make sure things are formatted correctly
  currentGamma <- matrix(currentGamma, nrow = 1)

  #Set proposal covariance matrix
  sig <- diag(1/rhoGamma, ncol(currentGamma))

  #Assuming Zi is (num Individuals)x(dimension of gamma) matrix of covariates

  #Get proposal draw
  propGamma <- mvtnorm::rmvnorm(1, mean = currentGamma,
                                sigma = sig)

  #Calculate thetais under proposal and current, note we are using the log-link, not the canonical
  propThetas <- Z %*% t(propGamma)
  propThetas <- exp(propThetas)

  currentThetas <- Z %*% t(currentGamma)
  currentThetas <- exp(currentThetas)

  #Calculate qs
  qOld <- mvtnorm::dmvnorm(currentGamma, mean = propGamma, sigma = sig)
  qNew <- mvtnorm::dmvnorm(propGamma, mean = currentGamma, sigma = sig)

  #Calculate ps
  pOld <- (sum(log(dgamma(taui, shape = currentThetas*phi, rate = phi))))
  pNew <- (sum(log(dgamma(taui, shape = propThetas*phi, rate = phi))))

  #Calculate a
  a <- max(0, min(1, exp(pNew-pOld)*(qOld/qNew)))

  #Flip a coin and return
  if(as.logical(rbinom(1,1,a))){
    return(list(Gamma = propGamma,
                thetai = propThetas))
  }else{
    return(list(Gamma = currentGamma,
                thetai = currentThetas))
  }
}

#' Metropolis-Hastings step to draw a new value for phi.
#'
#' In our model the data are drawn from LogN(mu_i + log(c_ij), tau_i). The prior for tau_i
#' is given as Gamma(thetai*phi, phi). This function uses a MH step to draw a new sample of phi.
#' Proposal distribution is Gamma(currentPhi*rhoPhi, rhoPhi).
#' Note that we parameterize with RATE, not SCALE.
#'
#' @param taui Numeric vector, individuals precisions.
#' @param thetai Numeric vector. individuals precisions means (estimate)
#' @param currentPhi Previous draw of phi
#' @param rhoPhi Proposal rate for gamma distribution that draws proposal for phi, default is 1000.
#'
#' @return Numeric, new draw of phi
postPhi <- function(taui, thetai, currentPhi, rhoPhi = 1000){
  #Draw proposal phi
  propPhi <- rgamma(1, currentPhi*rhoPhi, rate = rhoPhi)

  #Calculate nuis
  currentNus <- thetai*currentPhi
  propNus <- thetai*propPhi

  #Calculate qs
  qOld <- dnorm(currentPhi, mean = propPhi, sd = rhoPhi)
  qNew <- dnorm(propPhi, mean = currentPhi, sd = rhoPhi)

  #Calculate ps
  pOld <- sum(log(dgamma(taui, shape = currentNus, rate = currentPhi)))
  pNew <- sum(log(dgamma(taui, shape = propNus, rate = propPhi)))

  #Calculate a
  a <- max(0, min(1, exp(pNew-pOld)*(qOld/qNew)))

  #Flip a coin and return
  if(as.logical(rbinom(1,1,a))){
    return(propPhi)
  }else{
    return(currentPhi)
  }
}

#' Sample a value from the full conditional posterior of rho
#'
#' In our model the data are drawn from LogN(mu_ij + log(c_ij), tau_i). mu_ij = Xi^T %*%Beta + bi, where the prior
#' for bi is given as N(0, rho). This function draws from the conditional posterior of rho, given
#' that the prior on rho is a uniform prior on the standard deviation.
#'
#' @param b Numeric vector, Random intercepts for individuals
#'
#' @return Numeric
#'
postRho <- function(b){
  #n is the length of b
  n <- length(b)

  #Set posterior parameters
  postA <- (n-1)/2
  postB <- sum(b^2)/2

  #Draw value for rho and return
  return(rgamma(1, shape = postA, rate = postB))
}

#' Sample a value from the full conditional posterior of b_i
#'
#' @param taui Numeric, tau_i for individual i
#' @param mi Numeric vector, length equal to the number of cycles contributed by individual i
#' @param XiBeta Numeric vector, equal to t(X_i) %*% Beta
#' @param rho Numeric, rho value for prior precision of b_i
#'
#' @return Numeric value, single draw of bi
postB <- function(taui, mi, XiBeta, rho){
  #Number of cycles for individual i
  ni <- length(mi)

  #Calculate posterior parameters
  postPrec <- (rho + ni*taui)
  postMean <- taui*sum(mi - XiBeta)/postPrec

  return(rnorm(1, postMean, sqrt(1/postPrec)))
}

#' Sample a value from the full conditional posterior of tau_i
#'
#' In our model the data are drawn from LogN(mu_i + log(c_ij), tau_i). The prior for tau_i
#' is given as Gamma(thetai*phi, phi). This function draws from the conditional
#' posterior of tau_i. Note that we parameterize with RATE, not SCALE.
#'
#' Additionally, note that in order to vectorize the remainder of the MCMC algorithm
#' this function returns the sampled value repeated for length(yij)
#'
#' @param yij Numeric vector, cycle lengths for a single individual
#' @param cij Positive Integer vector, a sampled vector of length(yij) where the corresponding
#'  values in cij indicate a sampled number of TRUE cycles in each cycle length given by yij
#' @param mui Numeric, log of sampled mean of this individual's yijs
#' @param thetai Numeric, mean of prior (gamma) distribution on taui
#' @param phi Numeric, rate for Taui prior
#'
#' @return Numeric vector, repeated sampled value of length(yij)
postTaui <- function(yij, cij, mui, thetai, phi = 1){
  #Ni is the length of yij
  Ni <- length(yij)

  #Set posterior parameters
  postA <- thetai*phi + Ni/2
  postB <- phi + sum((log(yij/cij) - mui)^2)/2

  #Draw from posterior and return
  dr <- rgamma(1, shape = postA, rate = postB)
  #return(rep(dr, Ni))
  return(dr)
}

#' Sample a vector of values from the full conditional posterior of the c_ij vector
#'
#' In our model the data are drawn from LogN(mu_i + log(c_ij), tau_i) The prior for
#'  c_ij is a categorical prior with category probabilities pi1, ..., pik, and c_ij can
#'  take values 1, ..., k where k is the length of pi. This function samples from the
#'  full conditional posterior of all c_ijs, given vectors of equal length yijs, muis, tauis
#'
#' @param yijs Numeric Vector, cycle lengths
#' @param pi Numeric vector, must sum to 1. Sampled probabilities for c_ijs
#' @param muis Numeric vector, log of sampled mean for all individuals yijs
#' @param tauis Numeric vector > 0, sampled precision for all individuals yijs
#'
#' @return Integer vector
#'
postCij <- function(yijs, pi, muis, tauis){
  probs <- sapply(1:length(pi), function(j){
    #likelihood
    lik <- dlnorm(yijs, meanlog = muis + log(j), sdlog = sqrt(1/tauis))
    return(lik*pi[j])
  })

  #Set minimum at machine 0. Prevents the situation where all probabilities are 0 and we get
  #an error
  probs <- pmax(probs, rep(5e-324,length(probs)))

  return(glmnet::rmult(probs))
}

#' Sample a value from the full conditional posterior of pi
#'
#' In our model the data are drawn from LogN(mu_i + log(c_ij), tau_i) The prior for
#'  c_ij is a categorical prior with category probabilities pi1, ..., pik, and c_ij can
#'  take values 1, ..., k where k is the length of pi. This function samples from the
#'  posterior of pi = pi1, ..., pik, assuming that pi follows Dirichlet(priorAlphas)
#'
#' @param ci Integer vector, all of the sampled cij values for all individuals
#' @param priorAlphas Numeric vector, prior dirichlet parameters for pi
#'
#' @return Numeric vector
#'
postPi <- function(ci, priorAlphas){
  #Get posterior dirichlet parameters
  postAlphas <- sapply(1:length(priorAlphas), function(k){
    return(priorAlphas[k] + sum(ci == k))
  })

  #Sample from dirichlet and return
  return(as.numeric(LaplacesDemon::rdirichlet(1, alpha = postAlphas)))
}
