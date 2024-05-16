#Line to get CRAN to accept package, (ggplot causes problems)
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ys", "cijs", 'taus',
                                                        'count', 'Beta', 'Lower',
                                                        'Upper', 'chain'))

#' Visualize Results from skipTrack.fit
#'
#' This function takes results from skipTrack.fit and produces several helpful visualizations.
#'
#' @param stFit A list containing MCMC results obtained from skipTrack.fit.
#'
#' @return A list of three ggplot2 objects:
#'   \itemize{
#'     \item cijOverLength - Scatter plot of estimated Cij values against reported cycle length.
#'     \item cijOverTaus - Scatter plot of estimated Cij values against estimated individual precisions, colored by cycle length.
#'     \item cijDens - Density plot of Y values overlayed with a density plot of Y values separated by estimated cij value.
#'   }
#'
#' @examples
#' #Simulated data
#' simDat <- skipTrack.simulate(n = 100, skipProb = c(.7, .2, .1))
#'
#' #Run model fit (should typically run with much more than 50 reps)
#' modFit <- skipTrack.fit(Y = simDat$Y, cluster = simDat$cluster, chains = 2, reps = 50)
#'
#' #Visualize results
#' skipTrack.visualize(modFit)
#'
#' @seealso
#' \code{\link{skipTrack.fit}} for generating MCMC results.
#'
#' @export
#'
skipTrack.visualize <- function(stFit){
  #Get fit results
  stFit <- stFit$fit

  #Collect important variables into dataframes

  #Creates a dataframe with chain/draw specific rhos
  rhoDF <- lapply(1:length(stFit), function(chainI){
    chainIrhos <- sapply(stFit[[chainI]], getElement, 'rho')
    return(data.frame('t' = 1:length(chainIrhos), 'rho' = chainIrhos, 'chain' = chainI))
  })
  rhoDF <- do.call('rbind', rhoDF)

  #Creates a dataframe with chain/draw specific betas and gammas
  betaDF <- lapply(1:length(stFit), function(chainI){
    #Extract and arrange
    chainIbetas <- sapply(stFit[[chainI]], getElement, 'Beta')
    if(is.matrix(chainIbetas)){
      chainIbetas <- t(chainIbetas)
    }else{
      chainIbetas <- matrix(chainIbetas, ncol = 1)
    }

    #Add iteration and chain info
    m <- cbind(matrix(1:nrow(chainIbetas), ncol = 1),
               chainIbetas, matrix(chainI, nrow = nrow(chainIbetas)))

    #Turn into df and name correctly
    df <- as.data.frame(m)
    names(df) <- c('t', paste0('Beta', 0:(ncol(chainIbetas)-1)), 'chain')
    return(df)
  })
  betaDF <- do.call('rbind', betaDF)

  gammaDF <- lapply(1:length(stFit), function(chainI){
    chainIgammas <- t(sapply(stFit[[chainI]], getElement, 'Gamma'))
    m <- cbind(matrix(1:nrow(chainIgammas), ncol = 1),
               chainIgammas, matrix(chainI, nrow = nrow(chainIgammas)))
    df <- as.data.frame(m)
    names(df) <- c('t', paste0('Gamma', 0:(ncol(chainIgammas)-1)), 'chain')
    return(df)
  })
  gammaDF <- do.call('rbind', gammaDF)


  #Creates a dataframe with chain/draw specific muis and tauis
  iList <- lapply(1:length(stFit), function(chainI){
    chainImus <- lapply(1:length(stFit[[chainI]]), function(t){
      mus <- stFit[[chainI]][[t]]$iDat$mus
      return(mus)
    })
    chainItaus <- lapply(1:length(stFit[[chainI]]), function(t){
      taus <- stFit[[chainI]][[t]]$iDat$taus
      return(taus)
    })
    chainImus <- rowMeans(do.call('cbind', chainImus))
    chainItaus <- rowMeans(do.call('cbind', chainItaus))
    chainIiDF <- stFit[[chainI]][[1]]$iDat
    chainIiDF$mus <- chainImus
    chainIiDF$taus <- chainItaus
    chainIiDF$chain <- chainI
    return(chainIiDF)
  })
  overallmuis <- lapply(iList, function(chain){
    return(chain$mus)
  })
  overalltauis <- lapply(iList, function(chain){
    return(chain$taus)
  })
  overallmuis <- rowMeans(do.call('cbind', overallmuis))
  overalltauis <- rowMeans(do.call('cbind', overalltauis))
  iOverDF <- iList[[1]]
  iOverDF$mus <- overallmuis
  iOverDF$taus <- overalltauis
  iOverDF$chain <- NULL

  iChainDF <- do.call('rbind', iList)

  #Creates a dataframe with chain specific average cijs,
  #and a dataframe with overall cijs
  ijList <- lapply(1:length(stFit), function(chainI){
    chainIcijs <- lapply(1:length(stFit[[chainI]]), function(t){
      cijs <- stFit[[chainI]][[t]]$ijDat$cijs
      return(cijs)
    })
    chainIcijs <- do.call('cbind', chainIcijs)
    chainIijDF <- stFit[[chainI]][[1]]$ijDat
    chainIijDF$cijs <- rowMeans(chainIcijs)
    chainIijDF$chain <- chainI
    return(chainIijDF)
  })

  overallcijs <- lapply(ijList, function(chain){
    return(chain$cijs)
  })
  overallcijs <- rowMeans(do.call('cbind', overallcijs))
  cijOverDF <- ijList[[1]]
  cijOverDF$cijs <- overallcijs

  #Add muis and tauis to cijOverDF
  cijOverDF$mus <- sapply(cijOverDF$Individual, function(ind){
    return(iOverDF$mus[iOverDF$Individual == ind])
  })
  cijOverDF$taus <- sapply(cijOverDF$Individual, function(ind){
    return(iOverDF$taus[iOverDF$Individual == ind])
  })

  cijChainDF <- do.call('rbind', ijList)

  ################
  #Make some plots
  ################

  #Cijs vs cycle length overall
  cijOverLength <- ggplot2::ggplot(data = cijOverDF,
                                ggplot2::aes(x = ys, y = cijs)) + ggplot2::geom_point()
  cijOverLength <- cijOverLength + ggplot2::theme_minimal() + ggplot2::ggtitle('Estimated Cijs vs. Ys')
  cijOverLength <- cijOverLength + ggplot2::xlab('Reported Cycle Length') + ggplot2::ylab('Estimated Cij Values')
  cijOverLength <- cijOverLength + ggplot2::theme(plot.title = ggtext::element_textbox_simple())

  #Cijs vs tauis overall
  cijOverTaus <- ggplot2::ggplot(data = cijOverDF,
                                ggplot2::aes(x = taus, y = cijs, col = ys)) + ggplot2::geom_point()
  cijOverTaus <- cijOverTaus + ggplot2::theme_minimal() + ggplot2::ggtitle('Estimated Cijs vs. Estimated Tauis')
  cijOverTaus <- cijOverTaus + ggplot2::xlab('Estimated Individual Precisions') + ggplot2::ylab('Estimated Cij Values')
  cijOverTaus <- cijOverTaus + ggplot2::labs(col = 'Cycle Length')
  cijOverTaus <- cijOverTaus + ggplot2::theme(plot.title = ggtext::element_textbox_simple())

  #Density of Overall Mean by Chain
  #muByChain <- ggplot2::ggplot(data = muDF,
  #                             ggplot2::aes(x = exp(mu),
  #                                          col = as.factor(chain), fill = as.factor(chain)))
  #muByChain <- muByChain + ggplot2::geom_density(alpha = .3) + ggplot2::theme_minimal()
  #muByChain <- muByChain + ggplot2::ggtitle('Sampled Density of Overall Mean by Chain') + ggplot2::xlab('Overall Mean')
  #muByChain <- muByChain + ggplot2::theme(legend.position = 'none')

  #Density of cycle lengths categorized by skip vs density of overall
  cijDens <- ggplot2::ggplot(data = cijOverDF,
                             ggplot2::aes(x = ys)) + ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(count)))
  cijDens <- cijDens + ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(count),
                                                          fill = as.factor(round(cijs))),
                                             alpha = .4, color = NA)
  cijDens <- cijDens + ggplot2::theme_minimal() + ggplot2::theme(legend.position = 'none',
                                                                 plot.title = ggtext::element_textbox_simple())
  cijDens <- cijDens + ggplot2::ggtitle('Cycle Density By Skip Categories')

  #Overall 95% Posterior Intervals for betas
  betaQuants <- lapply(0:(ncol(betaDF)-3), function(i){
    ql <- quantile(betaDF[,paste0('Beta', i)], .05)
    qu <- quantile(betaDF[,paste0('Beta', i)], .95)
    return(data.frame('Lower' = ql, 'Upper' = qu, 'Beta' = i))
  })
  betaQuants <- do.call('rbind', betaQuants)

  betaQuantsChain <- lapply(0:(ncol(betaDF)-3), function(i){
    byChain <- lapply(1:max(betaDF$chain), function(ch){
      ql <- quantile(betaDF[betaDF$chain == ch,paste0('Beta', i)], .05)
      qu <- quantile(betaDF[betaDF$chain == ch,paste0('Beta', i)], .95)
      return(data.frame('Lower' = ql, 'Upper' = qu, 'Beta' = i, 'chain' = ch))
    })
    return(do.call('rbind', byChain))
  })
  betaQuantsChain <- do.call('rbind', betaQuantsChain)

  betaPI <- ggplot2::ggplot(data = betaQuants) + ggplot2::geom_errorbar(ggplot2::aes(x = Beta,
                                                                                     ymin = Lower,
                                                                                     ymax = Upper))
  betaPI <- betaPI + ggplot2::theme_minimal() + ggplot2::ggtitle('Estimated Beta Posterior Intervals')

  betaPIch <- ggplot2::ggplot(data = betaQuantsChain) + ggplot2::geom_errorbar(ggplot2::aes(x = Beta,
                                                                                     ymin = Lower,
                                                                                     ymax = Upper,
                                                                                     col = as.factor(chain)),
                                                                               position = 'dodge')
  betaPIch <- betaPIch + ggplot2::theme_minimal() + ggplot2::ggtitle('Estimated Beta Posterior Intervals - By Chain')
  betaPIch <- betaPIch + ggplot2::theme(legend.position = 'none')


  return(list(cijOverLength, cijOverTaus, cijDens))
}
