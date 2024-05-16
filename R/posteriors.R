#This document holds all of the functions that provide random draws from full conditional
#posteriors for the 'duttweiler' algorithm.

#' Draw from Posterior Distribution for Beta Parameters
#'
#' In our model mui follows a normal distribution with mean Xi^T %*% beta and precision rho. Additionally
#' we assume that beta follows a mvnormal prior with mean 0 and precision (rho_Beta) * I.This function draws
#' from the posterior distribution of beta under these assumptions.
#'
#' @param rhoBeta A scalar representing the prior precision parameter for beta.
#' @param rho A scalar representing the precision parameter.
#' @param Xi A matrix of covariates, where each row represents an individual and each column represents a covariate.
#' @param muI A vector where each element is the mean for individual i.
#' @return A vector representing a draw from the posterior distribution of beta parameters.
#'
#' @details This function assumes that \code{Xi} is a (num Individuals) x (dimension of beta) matrix of covariates.
#'
postBeta <- function(rhoBeta = .01, rho, Xi, muI){
  #Assuming Xi is (num Individuals)x(dimension of beta) matrix of covariates
  XTX <- t(Xi) %*% Xi
  Xmu <- t(Xi) %*% muI

  postPre <- (rhoBeta)*diag(1, dim(XTX)[1]) + rho*XTX
  postVar <- solve(postPre)
  postMean <- postVar %*% (rho*Xmu)

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
#' @param Zi A matrix of covariates, where each row represents an individual and each column represents a covariate.
#' @param currentGamma A vector (or matrix with 1 row) representing the current Gamma value.
#' @param phi A scalar, the prior rate for tau_i.
#' @param rhoGamma A scalar representing the proposal distribution precision parameter.
#'
#' @details This draw step assumes a log-link function for the Gamma GLM that we are fitting.
#'
#' @return A list containing the new Gamma value and the corresponding thetai values.
#'
postGamma <- function(taui, Zi, currentGamma, phi = 1, rhoGamma = 1000){
  #make sure things are formatted correctly
  currentGamma <- matrix(currentGamma, nrow = 1)

  #Set proposal covariance matrix
  sig <- diag(1/rhoGamma, ncol(currentGamma))

  #Assuming Zi is (num Individuals)x(dimension of gamma) matrix of covariates

  #Get proposal draw
  propGamma <- mvtnorm::rmvnorm(1, mean = currentGamma,
                                sigma = sig)

  #Calculate thetais under proposal and current, note we are using the log-link, not the canonical
  propThetas <- Zi %*% t(propGamma)
  propThetas <- exp(propThetas)

  currentThetas <- Zi %*% t(currentGamma)
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
#' In our model the data are drawn from LogN(mu_i + log(c_ij), tau_i). The prior for mu_i
#' is given as N(mu, rho). This function draws from the conditional posterior of rho, given
#' that the prior on rho is a uniform prior on the standard deviation.
#'
#' @param muI Numeric vector, log of individuals mean values.
#' @param xib Numeric vector, result of X %*% Beta, same length as muI.
#'
#' @return Numeric
#'
postRho <- function(muI, xib){
  #n is the length of muI
  n <- length(muI)

  #Set posterior parameters
  postA <- (n-1)/2
  postB <- sum((muI - xib)^2)/2

  #Draw value for rho and return
  return(rgamma(1, shape = postA, rate = postB))
}

#' Sample a value from the full conditional posterior of mu_i
#'
#' In our model the data are drawn from LogN(mu_i + log(c_ij), tau_i). The prior for mu_i
#' is given as N(x_i^T %*% beta, rho). This function draws from the conditional posterior of mu_i.
#'
#' Additionally, note that in order to vectorize the remainder of the MCMC algorithm
#' this function returns the sampled value repeated for length(yij)
#'
#' @param yij Numeric vector, cycle lengths for a single individual
#' @param cij Positive Integer vector, a sampled vector of length(yij) where the corresponding
#'  values in cij indicate a sampled number of TRUE cycles in each cycle length given by yij
#' @param taui Numeric > 0, A sampled precision for the yijs
#' @param xib Numeric, result of multiplying x_i^T %*% beta (single value, not vector)
#' @param rho Numeric > 0, sampled prior precision of mu_i
#'
#' @return Numeric vector, repeated sampled value of length(yij)
#'
postMui <- function(yij, cij, taui, xib, rho){
  #Ni is the length of yij
  Ni <- length(yij)

  #Set posterior mean and precision
  postPre <- Ni*taui + rho
  postMean <- (rho*xib + taui*sum(log(yij/cij)))/postPre

  #Draw from posterior and return
  dr <- rnorm(1, mean = postMean, sd = sqrt(1/postPre))
  return(rep(dr, Ni))
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
  return(rep(dr, Ni))
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
