#' Monte Carlo estimate of negative marginal log-likelihood of Li model
#'
#' This function calculates a Monte Carlo estimate of the negative marginal log-likelihood
#' of the given hyperparameters for the generative model from Li et al. (2022).
#' It samples M instances of the parameters from the given distributions
#' and averages the the likelihoods, giving a marginal likelihood for the hyperparameters.
#'
#' @param pars Named numeric vector of hyperparameters containing the
#' elements: kappa, gamma, alpha, beta. NOTE: MUST BE IN CORRECT ORDER.
#'   - \code{kappa}: Numeric value, shape parameter of Gamma distribution for Lambda_i.
#'   - \code{gamma}: Numeric value, rate parameter of Gamma distribution for Lambda_i.
#'   - \code{alpha}: Numeric value, shape1 parameter of Beta distribution for Pi_i.
#'   - \code{beta}:  Numeric value, shape2 parameter of Beta distribution for Pi_i.
#' @param S Integer, maximum number of allowed skips in the model.
#' @param M Integer specifying the number of Monte Carlo iterations.
#' @param cycleDat Data frame containing information about individuals and their tracked cycles.
#' @param verbose Logical with default FALSE. If true, prints extra info while running.
#' @param ... Does nothing.
#'
#' @return Numeric value representing the Monte Carlo estimate of the negative marginal log-likelihood.
#'
#' @references Li, Kathy, et al. "A predictive model for next cycle start date that accounts for adherence in menstrual self-tracking." Journal of the American Medical Informatics Association 29.1 (2022): 3-11.

likVec <- function(pars = c(kappa = 180,
                            gamma = 6,
                            alpha = 2,
                            beta = 20), S = 10, M = 1000, cycleDat, verbose = FALSE, ...){
  if(verbose){
    message(
      paste0('\nParameters: ', pars)
    )
  }
  kappa <- max(pars[1], 0)
  gamma <- max(pars[2], 0)
  alpha <- max(pars[3], 0)
  beta <- max(pars[4], 0)

  #Simulate parameters one per M iteration
  numCyc <- nrow(cycleDat)

  lamIs <- rgamma(M, shape = kappa, rate = gamma)
  piIs <- rbeta(M, shape1 = alpha, shape2 = beta)

  #Get s+1 raised to cycles
  sTOcyc <- rep(1:(S+1), times = rep(numCyc, S+1))^cycleDat$TrackedCycles
  sTOcyc <- rep(sTOcyc, M)[order(rep(rep(1:(S+1), times = rep(numCyc, S+1)), M))]

  #Get pi*exp(-lam)^S
  pilamTOs <- (piIs*exp(-lamIs))^rep(0:S, times = rep(M, S+1))
  pilamTOs <- rep(pilamTOs, times = rep(numCyc, M*(S+1)))

  #Repeat cycle vec M times, and lambda/pi unique pairs, cycle number of times
  cycM <- rep(cycleDat$TrackedCycles, M)
  lamM <- rep(lamIs, times = rep(numCyc, M))
  piM <- rep(piIs, times = rep(numCyc, M))

  #TRY 3, calculate first half, then repeat to calculate second half, then combine
  firstHalf <- ((1-piM)/(1-piM^(S+1)))*dpois(cycM, lamM)

  secondHalf <- rowSums(matrix(pilamTOs*sTOcyc, ncol = S+1))

  lcijmLiks <- log(firstHalf*secondHalf)


  #Now, combine over j
  individualSplits <- as.logical(diff(rep(cycleDat$Individual, M)))
  dimLiks <- exp(diff(c(0,
                        cumsum(lcijmLiks)[individualSplits],
                        sum(lcijmLiks))))

  #Average over Ms
  diLiks <- rowMeans(matrix(dimLiks, ncol = M))

  #Return sum of negative log likelihoods per individual
  negLLik <- -sum(log(diLiks))

  if(verbose){
    message(paste0('\nCurrent Value: ',negLLik))
  }

  return(negLLik)
}

#' Perform hyperparameter inference assuming the model given in Li et al. (2022) on a cycle length dataset.
#'
#' This function performs hyperparameter inference on a given dataset of individuals
#' and their tracked cycles, assuming the model specified in Li et al. (2022).
#' Default starting values for hyperparameters and optimization tuning parameters
#' are those given in Li et al.
#'
#' @param Y A vector of observed cycle lengths.
#' @param cluster A vector indicating the individual cluster/group membership for each observation Y.
#' @param S Maximum number of possible skipped cycles (see Li et al. for details).
#' @param startingParams A vector of starting values for hyperparameters
#'                      (default values from Li et al.).
#'
#' @return A list containing the results of hyperparameter inference.
#'
#' @references Li, Kathy, et al. "A predictive model for next cycle start date that accounts for adherence in menstrual self-tracking." Journal of the American Medical Informatics Association 29.1 (2022): 3-11.
liInference <- function(Y,cluster, S = 10,
                        startingParams = c(kappa = 180,
                                              gamma = 6,
                                              alpha = 2,
                                              beta = 20)){
  #Create cycleDat
  cycleDat <- data.frame('TrackedCycles' = Y, 'Individual' = cluster)

  #Perform inference using ADAM
  ret <- optimg::optimg(par = startingParams, fn = likVec, gr = NULL, S = S, cycleDat = cycleDat,
                        method = 'ADAM', control = list(reltol = .001),
                        maxit = 1000, tol = 1e-3, Interval = .01, verbose = FALSE)
  par <- ret$par
  names(par) <- c('kappa', 'gamma', 'alpha', 'beta')
  return(par)
}
