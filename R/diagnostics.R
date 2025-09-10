#' skipTrack MCMC Diagnostics
#'
#' Takes model results from skipTrack.fit and uses genMCMCDiag to get generalized mcmc diagnostics
#'
#' @param stFit A list of MCMC results from the skipTrack.fit function.
#' @param param A character string specifying the parameter for which diagnostics are to be calculated.
#'   Must be one of: 'rho', 'phi', 'Betas', 'Gammas', 'muis', 'tauis', or 'cijs'.
#' @param proximityMap An optional parameter specifying the proximity-map for calculating diagnostics. See package genMCMCDiag for details. Default is NULL.
#' @param ... Additional parameters to be passed to the genDiagnostic function.
#' @inheritDotParams genMCMCDiag::genDiagnostic diagnostics distance verbose
#'
#' @return A mcmcDiag object of MCMC diagnostics for the specified parameter
#'
#' @details If the parameter is 'rho' or 'phi' (the univariate parameters),
#'   the function extracts the specified parameter from the MCMC results and calculates
#'   diagnostics using the genDiagnostic function with the
#'   standard proximityMap. If the parameter is any of the other available options, the
#'   function extracts the corresponding values and calculates diagnostics using the genDiagnostic
#'   function with the specified or default proximityMap ('lanfear') and hammingDist as the distance function.
#'
#'   Details on genDiagnostic can be found in the genMCMCDiag package.
#'
#' @seealso \code{\link[genMCMCDiag:genDiagnostic]{genDiagnostic}}, \code{\link{skipTrack.fit}}
#'
#' @export
#'
#' @examples
#' #Simulated data
#' simDat <- skipTrack.simulate(n = 100, skipProb = c(.7, .2, .1))
#'
#' #Run model fit (should typically run with much more than 50 reps)
#' modFit <- skipTrack.fit(Y = simDat$Y, cluster = simDat$cluster, chains = 2, reps = 50)
#'
#' #Get diagnostics for cijs
#' skipTrack.diagnostics(modFit, 'cijs')
skipTrack.diagnostics <- function(stFit, param = c('rho', 'phi', 'Betas', 'Gammas', 'tauis', 'cijs'),
                                  proximityMap = NULL, ...){
  #If class is skipTrack.model, extract chains, otherwise assume we have chains
  if('skipTrack.model' %in% class(stFit)){
    #Get fit results
    stFit <- stFit$fit
  }

  #select first parameter
  param <- param[1]

  if(param %in% c('rho', 'phi')){#If param is one of the univariate parameters, get diagnostics

    #Extract specified parameter
    mcmcExt <- lapply(stFit, function(chain){
      #Get draws of specified parameter
      draws <- sapply(chain, getElement, name = param)

      #Return in expected format
      return(draws)
    })

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, proximityMap = 'standard', ...))

  }else if(param == 'cijs'){ #Continued alternative proximityMaps specific to skipTrack

    #Extract cijs
    mcmcExt <- lapply(stFit, function(chain){
      #Get list of cij draws
      draws <- lapply(chain, function(d){
        return(d$ijDat$cijs)
      })

      #Return in expected format
      return(draws)
    })

    #set proximityMap if not specified
    if(is.null(proximityMap)){
      proximityMap <- 'lanfear'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, proximityMap = proximityMap,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'Betas'){ #Continued alternative proximityMaps specific to skipTrack

    #Extract cijs
    mcmcExt <- lapply(stFit, function(chain){
      #Get list of Beta draws
      draws <- lapply(chain, function(d){
        return(d$Beta)
      })

      #Return in expected format
      return(draws)
    })

    #set proximityMap if not specified
    if(is.null(proximityMap)){
      proximityMap <- 'lanfear'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, proximityMap = proximityMap,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'Gammas'){ #Continued alternative proximityMaps specific to skipTrack

    #Extract cijs
    mcmcExt <- lapply(stFit, function(chain){
      #Get list of Beta draws
      draws <- lapply(chain, function(d){
        return(d$Gamma)
      })

      #Return in expected format
      return(draws)
    })

    #set proximityMap if not specified
    if(is.null(proximityMap)){
      proximityMap <- 'lanfear'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, proximityMap = proximityMap,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'muijs'){
    #Extract muis
    mcmcExt <- lapply(stFit, function(chain){
      #Get list of mui draws
      draws <- lapply(chain, function(d){
        return(d$ijDat$muijs)
      })

      #Return in expected format
      return(draws)
    })

    #set proximityMap if not specified
    if(is.null(proximityMap)){
      proximityMap <- 'lanfear'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, proximityMap = proximityMap,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'tauis'){
    #Extract tauis
    mcmcExt <- lapply(stFit, function(chain){
      #Get list of taui draws
      draws <- lapply(chain, function(d){
        return(d$iDat$taus)
      })

      #Return in expected format
      return(draws)
    })

    #set proximityMap if not specified
    if(is.null(proximityMap)){
      proximityMap <- 'lanfear'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, proximityMap = proximityMap,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'sijs'){ #proximityMap for LI inference

    #Extract sijs
    mcmcExt <- lapply(stFit, function(chain){
      #Get list of sij draws
      draws <- lapply(chain, function(d){
        return(d$ijDat$ss)
      })

      #Return in expected format
      return(draws)
    })

    #set proximityMap if not specified
    if(is.null(proximityMap)){
      proximityMap <- 'lanfear'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, proximityMap = proximityMap,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else{ #Throw error if param is not recognized
    stop("param must be character string 'rho', 'phi', 'Betas', 'Gammas', 'muis', 'tauis', or 'cijs'")
  }
}
