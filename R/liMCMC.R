
#' Runs MCMC algorithm for performing inference using the model from Li et al. (2022)
#'
#' @description
#' This function performs inference on cycle length data, assuming the model from Li et al. (2022). It is important to note
#' that Li et al. does not actually use this algorithm as they target a particular analytic posterior predictive distribution,
#' and solve directly. However, we are targeting a different posterior and thus use this MCMC to perform inference.
#'
#' @inheritParams skipTrack.MCMC
#'
#' @param S Integer. The maximum number of skips to consider possible.
#' @param hyperparams Named numeric vector of hyperparameters containing the
#' elements: kappa, gamma, alpha, beta. NOTE: MUST BE IN CORRECT ORDER.
#'   - \code{kappa}: Numeric value, shape parameter of Gamma distribution for Lambda_i.
#'   - \code{gamma}: Numeric value, rate parameter of Gamma distribution for Lambda_i.
#'   - \code{alpha}: Numeric value, shape1 parameter of Beta distribution for Pi_i.
#'   - \code{beta}:  Numeric value, shape2 parameter of Beta distribution for Pi_i.
#' @param initialParams A list of initial parameter values for the MCMC algorithm.
#' Default values are provided for pi, lambdais, piis, ss.
#' @param ... For catching unused arguments (like li = TRUE)
#'
#' @return A list containing the MCMC draws for each parameter at each iteration. Each element
#' in the list is itself a list containing:
#' \describe{
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, lambdais, piis, ss.}
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, lambdas, pis.}
#'   \item{kappa}{Fixed value of hyperparameter kappa.}
#'   \item{gamma}{Fixed value of hyperparameter gamma.}
#'   \item{alpha}{Fixed value of hyperparameter alpha.}
#'   \item{beta}{Fixed value of hyperparamter beta.}
#'   \item{S}{Fixed input value S.}
#'   \item{indFirst}{A logical vector indicating the first occurrence of each individual.}
#' }
#'
#' @seealso \code{\link{gibbsStepLi}}
#'
#' @references Li, Kathy, et al. "A predictive model for next cycle start date that accounts for adherence in menstrual self-tracking." Journal of the American Medical Informatics Association 29.1 (2022): 3-11.
#'
liMCMC <- function(Y,
                   cluster,
                   S,
                   hyperparams = c(kappa = 180, gamma = 6, alpha = 2, beta = 20),
                   initialParams = list(pi = c(1/3, 1/3, 1/3),
                                               lambdais = rep(30,
                                                          length(unique(cycleDat$Individual))),
                                               piis = rep(.2,
                                                           length(unique(cycleDat$Individual))),
                                               ss = sample(0:S,
                                                           nrow(cycleDat),
                                                           replace = TRUE)),
                   reps = 1000, ...){
  #Create cycleDat
  cycleDat <- data.frame('TrackedCycles' = Y, 'Individual' = cluster)

  #Deal with initial params
  ip  <-  list(pi = c(1/3, 1/3, 1/3),
               lambdais = rep(30,
                              length(unique(cycleDat$Individual))),
               piis = rep(.2,
                          length(unique(cycleDat$Individual))),
               ss = sample(0:S,
                           nrow(cycleDat),
                           replace = TRUE))
  for(nm in names(ip)){
    if(nm %in% names(initialParams)){
      ip[nm] <- initialParams[nm]
    }
  }

  initialParams <- ip

  #Organize data into initial list
  iDat <- data.frame('Individual' = unique(cycleDat$Individual),
                     'lambdas' = initialParams$lambdais,
                     'pis' = initialParams$piis)
  ijDat <- data.frame('Individual' = cycleDat$Individual,
                      'ys' = cycleDat$TrackedCycles,
                      'ss' = initialParams$ss,
                      'lambdais' = sapply(cycleDat$Individual, function(ind){iDat$lambdas[iDat$Individual == ind]}),
                      'piis' = sapply(cycleDat$Individual, function(ind){iDat$pis[iDat$Individual == ind]}))
  fullDraws <- vector('list', reps + 1)
  fullDraws[[1]] <- list(ijDat = ijDat, iDat = iDat,
                         kappa = as.numeric(hyperparams['kappa']),
                         gamma = as.numeric(hyperparams['gamma']),
                         alpha = as.numeric(hyperparams['alpha']),
                         beta = as.numeric(hyperparams['beta']),
                         S = S,
                         indFirst = !duplicated(ijDat$Individual))
  #Progress bar
  pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)

  #Do gibbs steps
  for(t in 1:reps){
    fullDraws[[t+1]] <- do.call('gibbsStepLi', fullDraws[[t]])
    utils::setTxtProgressBar(pb, t)
  }
  return(fullDraws)
}

#' Gibbs Step Li - One MCMC step for the Li Model
#'
#'
#'
#' @param ijDat A data.frame with parameters at the individual-observation level: Individual, ys, lambdais, piis, ss.
#' @param iDat A data.frame with parameters at the individual level: Individual, lambdas, pis.
#' @param kappa Fixed value of hyperparameter kappa.
#' @param gamma Fixed value of hyperparameter gamma.
#' @param alpha Fixed value of hyperparameter alpha.
#' @param beta Fixed value of hyperparamter beta.
#' @param S Fixed input value S.
#' @param indFirst A logical vector indicating the first occurrence of each individual.
#'
#' @return A list containing one MCMC draws for each parameter. Elements are:
#' \describe{
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, lambdais, piis, ss.}
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, lambdas, pis.}
#'   \item{kappa}{Fixed value of hyperparameter kappa.}
#'   \item{gamma}{Fixed value of hyperparameter gamma.}
#'   \item{alpha}{Fixed value of hyperparameter alpha.}
#'   \item{beta}{Fixed value of hyperparamter beta.}
#'   \item{S}{Fixed input value S.}
#'   \item{indFirst}{A logical vector indicating the first occurrence of each individual.}
#' }
#'
#' @references Li, Kathy, et al. "A predictive model for next cycle start date that accounts for adherence in menstrual self-tracking." Journal of the American Medical Informatics Association 29.1 (2022): 3-11.
#'
gibbsStepLi <- function(ijDat, iDat, kappa, gamma, alpha, beta, S, indFirst){
  #Now i level
  newLambdais <- lapply(iDat$Individual, function(ind){
    postLambdai(yij = ijDat$ys[ijDat$Individual == ind],
                sij = ijDat$ss[ijDat$Individual == ind],
                priorK = kappa, priorG = gamma)
  })
  newLambdais <- do.call('c', newLambdais)
  newPiis <- lapply(iDat$Individual, function(ind){
    postPii(sij = ijDat$ss[ijDat$Individual == ind],
            currentPii = iDat$pis[iDat$Individual == ind],
            priorA = alpha, priorB = beta, S = S)
  })
  newPiis <- do.call('c', newPiis)

  iDatNew <- data.frame(Individual = iDat$Individual,
                        lambdas = newLambdais[indFirst],
                        pis = newPiis[indFirst])

  ijDatNew <- ijDat
  ijDatNew$lambdais <- newLambdais
  ijDatNew$piis <- newPiis

  ijDatNew$ss <- postSij(ijDatNew$ys, pii = ijDatNew$piis,
                         lambdai = ijDatNew$lambdais, S = S)

  return(list(ijDat = ijDatNew, iDat = iDatNew, kappa = kappa,
              gamma = gamma, alpha = alpha, beta = beta, S = S, indFirst = indFirst))
}
