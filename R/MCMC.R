#' Perform one chain of MCMC sampling for the skipTrack model.
#'
#' This function runs a single Markov Chain Monte Carlo (MCMC) chain to update parameters in
#' the skipTrack hierarchical model.
#'
#' @param Y A vector of observed cycle lengths.
#' @param cluster A vector indicating the individual cluster/group membership for each observation Y.
#' @param X A matrix (length(Y) x length(Beta)) of covariates for cycle length mean. Default is a vector of 1's.
#' @param Z A matrix (length(Y) x length(Gamma)) of covariates for cycle length precision. Default is a vector of 1's.
#' @param numSkips The maximum number of skips to allow. Default is 10.
#' @param reps The number of MCMC iterations (steps) to perform. Default is 1000.
#' @param fixedSkips If TRUE cycle skip information (cijs) is not updated in sample steps and the inputs are instead assumed to be true.
#' @param initialParams A list of initial parameter values for the MCMC algorithm.
#' Default values are provided for pi, muis, tauis, rho, cijs, alphas, Beta, Gamma, phi, rhoBeta, rhoGamma, and rhoPhi.
#' @param verbose logical. If true progress bars and additional info are printed to the console.
#'
#' @return A list containing the MCMC draws for each parameter at each iteration. Each element
#' in the list is itself a list containing:
#' \describe{
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, cijs, muis, tauis.}
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, mus, taus, thetas.}
#'   \item{rho}{Updated value of the global parameter rho.}
#'   \item{pi}{Updated value of the global parameter pi.}
#'   \item{Xi}{Matrix of covariates for cycle length mean.}
#'   \item{Z}{Matrix of covariates for cycle length precision.}
#'   \item{Beta}{Updated matrix of coefficients for cycle length mean.}
#'   \item{Gamma}{Updated matrix of coefficients for cycle length precision.}
#'   \item{priorAlphas}{Vector of prior alpha values for updating pi.}
#'   \item{indFirst}{A logical vector indicating the first occurrence of each individual.}
#'   \item{rhoBeta}{Hyperprior parameter rhoBeta, used to update Beta.}
#'   \item{rhoGamma}{Value of the proposal parameter rhoGamma.}
#'   \item{phi}{Updated value of the parameter phi.}
#'   \item{rhoPhi}{Value of the proposal parameter rhoPhi.}
#'   \item{fixedSkips}{Logical. Indicates if skips were fixed.}
#' }
#'
#' @seealso \code{\link{sampleStep}}
#'
skipTrack.MCMC <- function(Y,cluster,
                          X = matrix(1, nrow = length(cluster)),
                          Z = matrix(1, nrow = length(cluster)),
                          numSkips = 10,
                          reps = 1000,
                          fixedSkips = FALSE,
                          initialParams = list(pi = rep(1/(numSkips+1), numSkips+1),
                                               muis = rep(log(30),
                                                          length(unique(cluster))),
                                               tauis = rep(5,
                                                          length(unique(cluster))),
                                               bs = rep(0, length(unique(cluster))),
                                               rho = 1,
                                               cijs = rep(1, length(Y)),
                                               alphas = rep(1, numSkips +1),
                                               Beta = matrix(rep(0, ncol(as.matrix(X))),1),
                                               Gamma = matrix(rep(0, ncol(as.matrix(Z))),1),
                                               rhoBeta = .01,
                                               rhoGamma = 1000,
                                               phi = .01,
                                               rhoPhi = 1000),
                          verbose = FALSE){
  #Set initial params default list
  ip <- list(pi = rep(1/(numSkips+1), numSkips+1),
             muis = rep(log(30),
                        length(unique(cluster))),
             tauis = rep(5,
                         length(unique(cluster))),
             bs = rep(0, length(unique(cluster))),
             rho = 1,
             cijs = rep(1, length(Y)),
             alphas = rep(1, numSkips +1),
             Beta = matrix(rep(0, ncol(as.matrix(X))),1),
             Gamma = matrix(rep(0, ncol(as.matrix(Z))),1),
             rhoBeta = .01,
             rhoGamma = 1000,
             phi = .01, rhoPhi = 1000)

  #Replace anything that needs replacing
  for(i in 1:length(initialParams)){
    nm <- names(initialParams)[i]
    ip[nm] <- initialParams[nm]
  }

  #Replace whole thing
  initialParams <- ip

  #Checks for X and Z
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  #Make Z unique by individual (as currently expected by sampleStep)
  Z <- as.matrix(unique(cbind(cluster, as.data.frame(Z))))[,-1, drop = F]

  if(nrow(Z) != length(unique(cluster))){
    stop('Z cannot contain any time-varying covariates. Each row in Z associated with an individual must be identical to every other row for that individual.')
  }

  #Organize data into initial list
  iDat <- data.frame('Individual' = unique(cluster),
                     'taus' = initialParams$tauis,
                     'bs' = initialParams$bs,
                     'thetas' = exp(Z %*% t(initialParams$Gamma)))
  ijDat <- data.frame('Individual' = cluster,
                      'ys' = Y,
                      'cijs' = initialParams$cijs,
                      'mijs' = log(Y/initialParams$cijs),
                      'muijs' = sapply(cluster, function(ind){initialParams$muis[iDat$Individual == ind]}))
  fullDraws <- vector('list', reps + 1)
  fullDraws[[1]] <- list(ijDat = ijDat,
                         iDat = iDat,
                         rho = initialParams$rho,
                         pi = initialParams$pi,
                         X = X,
                         Z = Z,
                         Beta = initialParams$Beta,
                         Gamma = initialParams$Gamma,
                         priorAlphas = initialParams$alphas,
                         indFirst = !duplicated(ijDat$Individual),
                         rhoBeta = initialParams$rhoBeta,
                         rhoGamma = initialParams$rhoGamma,
                         phi = initialParams$phi,
                         rhoPhi = initialParams$rhoPhi,
                         fixedSkips = fixedSkips)
  #Progress bar
  if(verbose){pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)}

  #Do gibbs steps
  for(t in 1:reps){
    fullDraws[[t+1]] <- do.call('sampleStep', fullDraws[[t]])

    #Remove unecessary info from draws just used (to save space!!)
    fullDraws[[t]]$X <- NULL
    fullDraws[[t]]$Z <- NULL
    fullDraws[[t]]$indFirst <- NULL
    fullDraws[[t]]$ijDat$ys <- NULL

    if(verbose){utils::setTxtProgressBar(pb, t)}
  }
  return(fullDraws)
}

#' Perform a single step of the MCMC sampling process for skipTrack
#'
#' This function performs a single step of the Markov Chain Monte Carlo (MCMC) algorithm to update parameters
#' in a hierarchical model used for identifying skips in menstrual cycle tracking.
#'
#' @param ijDat A data.frame with individual-observation level parameters: Individual, ys, cijs, muis, tauis.
#' @param iDat A data.frame with individual level parameters: Individual, mus, taus, thetas.
#' @param rho Updated value of the global parameter rho.
#' @param pi Updated value of the global parameter pi.
#' @param X A matrix (numIndividuals x length(Beta)) of covariates for cycle length mean. Default is a vector of 1's. NOTE THE DIFFERENCE (from skipTrack.MCMC) IN EXPECTED DIMENSION OF X
#' @param Z A matrix (numIndividuals x length(Gamma)) of covariates for cycle length precision. Default is a vector of 1's. NOTE THE DIFFERENCE (from skipTrack.MCMC) IN EXPECTED DIMENSION OF Z
#' @param Beta Matrix (1 x length(Beta)) of coefficients for cycle length mean.
#' @param Gamma Matrix of (1 x length(Gamma)) coefficients for cycle length precision.
#' @param priorAlphas Vector of prior alpha values for updating pi.
#' @param indFirst A logical vector indicating the first occurrence of each individual.
#' @param rhoBeta Updated value of the global parameter rhoBeta.
#' @param rhoGamma Value of the proposal parameter rhoGamma.
#' @param phi Value of the parameter phi.
#' @param rhoPhi Value of the proposal parameter rhoPhi.
#' @param fixedSkips Logical. If TRUE cycle skip information (cijs) is not updated in sample steps and the inputs are instead assumed to be true.
#'
#' @return A list containing updated parameters after performing a single MCMC step.
#' The list includes:
#' \describe{
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, cijs, muis, tauis.}
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, mus, taus, thetas.}
#'   \item{rho}{Updated value of the global parameter rho.}
#'   \item{pi}{Updated value of the global parameter pi.}
#'   \item{X}{Matrix of covariates for cycle length mean.}
#'   \item{Z}{Matrix of covariates for cycle length precision.}
#'   \item{Beta}{Updated matrix of coefficients for cycle length mean.}
#'   \item{Gamma}{Updated matrix of coefficients for cycle length precision.}
#'   \item{priorAlphas}{Vector of prior alpha values for updating pi.}
#'   \item{indFirst}{A logical vector indicating the first occurrence of each individual.}
#'   \item{rhoBeta}{Hyperprior parameter rhoBeta, used to update Beta.}
#'   \item{rhoGamma}{Value of the proposal parameter rhoGamma.}
#'   \item{phi}{Updated value of the parameter phi.}
#'   \item{rhoPhi}{Value of the proposal parameter rhoPhi.}
#'   \item{fixedSkips}{Logical. Indicates if skips were fixed.}
#' }
#'
sampleStep <- function(ijDat, iDat, rho, pi,
                       X, Z, Beta, Gamma,
                       priorAlphas, indFirst,
                       rhoBeta, rhoGamma, phi, rhoPhi, fixedSkips){
  #Start with high level (without Gamma and thetais as those are connected)
  newBeta <- postBeta(rhoBeta = rhoBeta, X = X, b = iDat$b,
                      m = ijDat$mijs, tau = iDat$taus, indFirst = indFirst)

  #Set XiBeta based on newBeta
  XBeta <- X %*% t(newBeta)

  #Draw new bs (a good value of rho depends on bs being somewhat reasonable)
  newBs <- sapply(iDat$Individual, function(ind){
    postB(taui = iDat$taus[iDat$Individual == ind],
          mi = ijDat$mijs[ijDat$Individual == ind],
          XiBeta = XBeta[ijDat$Individual == ind,],
          rho = rho)
  })

  #Continue with high level information
  newRho <- postRho(b = newBs)
  newPi <- postPi(ci = ijDat$cijs, priorAlphas = priorAlphas)

  #Calculate new muijs
  newMuijs <- sapply(iDat$Individual, function(ind){
    XiBeta = XBeta[ijDat$Individual == ind,]
    return(XiBeta + newBs[iDat$Individual == ind])
  })
  newMuijs <- do.call('c', newMuijs)

  newTauis <- lapply(iDat$Individual, function(ind){
    postTaui(yij = ijDat$ys[ijDat$Individual == ind],
             cij = ijDat$cijs[ijDat$Individual == ind],
             mui = newMuijs[ijDat$Individual == ind],
             thetai = iDat$thetas[iDat$Individual == ind],
             phi = phi)
  })
  newTauis <- do.call('c', newTauis)

  #High level Gamma Things
  newGamList <- postGamma(taui = newTauis, Z = Z, currentGamma = Gamma, phi = phi,
                          rhoGamma = rhoGamma)
  newGamma <- newGamList$Gamma
  newThetas <- newGamList$thetai

  #High Level, Phi
  newPhi <- postPhi(taui = newTauis, thetai = newThetas, currentPhi = phi,
                    rhoPhi = rhoPhi)

  #Create new i level information
  iDatNew <- iDat
  iDatNew$bs <- newBs
  iDatNew$taus <- newTauis
  iDatNew$thetas <- newThetas

  #New ij level information
  ijDatNew <- ijDat
  ijDatNew$muijs <- newMuijs


  if(fixedSkips){
    ijDatNew$cijs <- ijDat$cijs
    ijDatNew$mijs <- ijDat$mijs
  }else{
    ijDatNew$cijs <- postCij(ijDatNew$ys, pi = newPi,
                           muis = ijDatNew$muijs,
                           tauis = sapply(ijDatNew$Individual,
                                          function(ind){iDatNew$taus[iDatNew$Individual == ind]}))
    ijDatNew$mijs <- log(ijDatNew$ys/ijDatNew$cijs)
  }

  return(list(ijDat = ijDatNew, iDat = iDatNew, rho = newRho,
              pi = newPi, X = X, Z = Z, Beta = newBeta,
              Gamma = newGamma, priorAlphas = priorAlphas, indFirst = indFirst,
              rhoBeta = rhoBeta, rhoGamma = rhoGamma, phi = newPhi, rhoPhi = rhoPhi,
              fixedSkips = fixedSkips))
}
