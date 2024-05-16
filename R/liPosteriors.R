#' Compute random draw from full conditional posterior for lambda_i in Li algorithm.
#'
#' This function calculates a random draw from the full conditional posterior distribution
#' for lambda_i in the Li algorithm, given the observed values y_ij, the indicators s_ij,
#' and the prior hyperparameters priorK and priorG.
#'
#' @param yij Vector of observed values for individual i.
#' @param sij Vector of cycle skip indicators for individual i.
#' @param priorK Prior hyperparameter kappa.
#' @param priorG Prior hyperparameter gamma.
#'
#' @return A random draw from the posterior distribution of lambda_i.
#'
postLambdai <- function(yij, sij, priorK, priorG){
  #n is the number of sampled values
  n <- length(yij)

  #set posterior K and G
  postK <- priorK + sum(yij)
  postG <- priorG + n + sum(sij)

  #Return gamma draw
  dr <- rgamma(1, shape = postK, rate = postG)
  return(rep(dr, n))
}

#' Compute M-H draw for pi_i in Li algorithm
#'
#' This performs a Metropolis-Hastings draw for pi_i, assuming s_ij follows a truncated geometric distribution with parameters
#' pi_i and S. The proposal distribution for pi_i is Beta(alpha, beta).
#'
#' @param sij Vector of cycle skip indicators for individual i
#' @param currentPii Current value of pi_i
#' @param priorA Hyperparameter alpha.
#' @param priorB Hyperparameter beta.
#' @param S Maximum number of skips allowed in algorithm
#'
#' @return Draw for pi_i, repeated for the number of observations from individual i
#'
postPii <- function(sij, currentPii, priorA, priorB, S){
  #n is the number of cycles for this individual
  n <- length(sij)
  postA <- priorA + sum(sij)
  postB <- priorB + n

  #Sample a possible pii
  possPii <- rbeta(1, postA, postB)

  #Eval M-H prob
  qOld <- dbeta(currentPii, postA, postB)
  qNew <- dbeta(currentPii, postA, postB)
  pOld <- exp(sum(log((currentPii^sij)*(1-currentPii)/(1-currentPii^(S+1)))))
  pNew <- exp(sum(log((possPii^sij)*(1-possPii)/(1-possPii^(S+1)))))
  a <- min(1, (pNew/pOld)*(qOld/qNew))

  #Flip coin to return new or old
  if(as.logical(rbinom(1, 1, a))){
    return(rep(possPii, n))
  }else{
    return(rep(currentPii, n))
  }
}

#' Compute random draw from full conditional posterior for s_ij in Li algorithm.
#'
#' This function calculates a random draw from the full conditional posterior distribution
#' for s_ij in the Li algorithm, given the observed values yijs, the parameter pi,
#' the lambda_i value, and the truncation level S.
#'
#' @param yijs Vector of observed values for s_ij.
#' @param pii Probability parameter pi.
#' @param lambdai Value of lambda_i.
#' @param S Truncation level.
#'
#' @return A random draw from the posterior distribution of s_ij.
#'
postSij <- function(yijs, pii, lambdai, S){
  probs <- sapply(0:S, function(s){
    #Get truncated geometric probability
    p <- (pii^s)*(1-pii)/(1-pii^(S+1))

    #Get likelihood
    lik <- dpois(yijs, lambda = lambdai*(s+1))

    return(p*lik)
  })
  probs <- pmax(probs, rep(0, length(probs)))

  return(glmnet::rmult(probs)-1)
}
