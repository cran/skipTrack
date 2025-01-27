#' Simulate user-tracked menstrual cycle data for multiple individuals
#'
#' This function generates synthetic data for user-tracked menstrual cycles given a generative model,
#' skip probabilities, maximum cycles and covariates (depending on the model).
#' It supports built-in models ('skipTrack', 'li', 'mixture') and custom models written as functions.
#'
#' @param n Number of individuals to simulate data for.
#' @param model model for data simulation. Can be a character ('skipTrack', 'li', 'mixture') or a custom function.
#' @param skipProb Vector of probabilities for number of true cycles per tracked cycle. For
#' example, (.7, .2, .1) means that 70% of observed cycles will contain one true cycle, 20%
#' will contain 2 true cycles and 10% will contain 3 true cycles. Default is NULL. If
#' model == 'li', skipProb values are set and user input will be ignored.
#' @param maxCycles Maximum number of cycles for generating skip cycles. Default is the length of skipProb.
#' If model == 'li', this must be specified; if model == 'skipTrack' or 'mixture', leave as default.
#' @param trueBetas Optional. True values for the mean regression coefficients (not counting intercept which is automatic based on the model).
#' @param trueGammas Optional. True values for the precision regression coefficients (not counting intercept which is automatic based on the model). Precision covariates not available for model == 'li'.
#' @param overlap Optional. Number of (non-intercept) columns shared between X and Z. Columns are shared from left to right.
#' @param avgCyclesPer Average number of cycles contributed by each individual. Actual number is drawn from Poisson for each person. Default is 7.
#'
#' @return A list containing:
#' \describe{
#'  \item{'Y'}{Tracked cycles from the simulated data.}
#'  \item{'cluster'}{Individual identifiers from the simulated data.}
#'  \item{'X'}{Covariate matrix for Betas (mean cycle length).}
#'  \item{'Z'}{Covariate matrix for Gammas (regularity).}
#'  \item{'Beta'}{True beta coefficients.}
#'  \item{'Gamma'}{True gamma coefficients.}
#'  \item{'NumTrue'}{Number of true cycles in each tracked cycle. Order matches Y.}
#'  \item{'Underlying'}{Subset of the simulated data containing individual level information. For 'skipTrack' - individual mean and precision for log(cycle lengths), for 'li' - individual mean for cycle lengths, for 'mixture' - individual mean for cycle lengths}
#'}
#'
#'
#' @examples
#' # Example simulation from the SkipTrack model
#' resultSt <- skipTrack.simulate(1000, model = 'skipTrack', skipProb = c(.7, .2, .1))
#' hist(resultSt$Y, breaks = 5:200)
#'
#' # Example simulation from the Li model
#' resultLi <- skipTrack.simulate(1000, model = 'li', maxCycles = 3)
#' hist(resultLi$Y, breaks = 5:200)
#'
#' #Example simulation from the mixture model
#' resultMix <- skipTrack.simulate(1000, model = 'mixture', skipProb = c(.7, .2, .1))
#' hist(resultMix$Y, breaks = 5:200)
#'
#' @seealso \code{\link{stSim}}, \code{\link{liSim}}, \code{\link{mixSim}}
#'
#' @references Li, Kathy, et al. "A predictive model for next cycle start date that accounts for adherence in menstrual self-tracking." Journal of the American Medical Informatics Association 29.1 (2022): 3-11.
#'
#' @export
skipTrack.simulate <- function(n,
                         model = c('skipTrack', 'li', 'mixture'),
                         skipProb = NULL,
                         maxCycles = length(skipProb),
                         trueBetas = NULL,
                         trueGammas = NULL,
                         overlap = 0,
                         avgCyclesPer = 7){

  #Checks for built in models, custom models are on their own.
  if(is.character(model)){
    #If model is li, make sure length(skipProb) is 1
    if((model[1] == 'li') & !is.null(skipProb)){
      stop('with model li, skipProb should be NULL and maxCycles should be specified.')
    }

    #If model is skipTrack, make sure maxCycles is length(skipProb)
    if((model[1] == 'skipTrack') & maxCycles != length(skipProb)){
      stop(
        'for model skipTrack, maxCycles should be the length of skipProb'
      )
    }
  }

  #n gives the number of individuals, for each individual simulate tracked cycles
  #given model
  if(is.character(model)){
    if(model[1] == 'skipTrack'){
      simDat <- lapply(1:n, stSim, skipProb = skipProb, maxCycles = maxCycles,
                       trueBetas = trueBetas,
                       trueGammas = trueGammas,
                       overlap = overlap,
                       avgCyclesPer = avgCyclesPer)
      simDat <- do.call('rbind', simDat)
    }else if(model[1] == 'li'){
      simDat <- lapply(1:n, liSim, skipProb = skipProb, maxCycles = maxCycles,
                       trueBetas = trueBetas, trueGammas = trueGammas,
                       avgCyclesPer = avgCyclesPer)
      simDat <- do.call('rbind', simDat)
    }else if(model[1] == 'mixture'){
      simDat <- lapply(1:n, mixSim,
                       skipProb = skipProb,
                       maxCycles = maxCycles,
                       trueBetas = trueBetas,
                       trueGammas = trueGammas,
                       overlap = overlap,
                       avgCyclesPer = avgCyclesPer)
      simDat <- do.call('rbind', simDat)
    }else{
      stop('specified model is unknown.')
    }
  }else if(is.function(model)){
    simDat <- lapply(1:n, model, skipProb = skipProb, maxCycles = maxCycles)
    simDat <- do.call('rbind', simDat)
  }else{
    stop('model must be one of the specified characters, or a function taking skipProb and maxCycles')
  }

  #Get X and Z
  X <- as.matrix(simDat[,grepl('X|Individual', names(simDat))])[,-1, drop = F]
  Z <- as.matrix(simDat[,grepl('Z|Individual', names(simDat))])[,-1, drop = F]
  #X <- unique(X)[,-1]
  #Z <- unique(Z)[,-1]

  #Get trueBetas and trueGammas based on settings
  if(is.null(trueBetas)){
    trueBetas <- simDat$Beta0[1]
  }else{
    trueBetas <- c(simDat$Beta0[1], trueBetas)
  }
  if(is.null(trueGammas)){
    trueGammas <- simDat$Gamma0[1]
  }else{
    trueGammas <- c(simDat$Gamma0[1], trueGammas)
  }

  #Return as list with specific components separated out
  return(list('Y' = simDat$TrackedCycles,
              'cluster' = simDat$Individual,
              'X' = X,
              'Z' = Z,
              'Beta' = trueBetas,
              'Gamma' = trueGammas,
              'NumTrue' = simDat$NumTrue,
              'Underlying' = simDat[, grepl('Mean|Prec', names(simDat)), drop = F]))
}

#' Simulate user tracked menstrual cycle data for an individual, based on the skipTrack model.
#'
#' This function generates synthetic data for user tracked menstrual cycles for a
#' single individual. For this model Beta_0 = log(30), Gamma_0 = 5.5, and phi = .01.
#'
#' @param i Individual identifier. Character, numeric or integer.
#' @param skipProb Vector of probabilities for number of true cycles per tracked cycle. For
#' example, (.7, .2, .1) means that 70% of observed cycles will contain one true cycle, 20%
#' will contain 2 true cycles and 10% will contain 3 true cycles.
#' @param maxCycles Maximum number of true cycles per tracked cycle. Ignored for this model.
#' @param trueBetas Optional. True values for the mean regression coefficients (not counting intercept which is automatic based on the model).
#' @param trueGammas Optional. True values for the precision regression coefficients (not counting intercept which is automatic based on the model).
#' @param overlap Optional. Number of (non-intercept) columns shared between X and Z. Columns are shared from left to right.
#' @param avgCyclesPer Average number of cycles contributed by each individual. Actual number is drawn from Poisson for each person. Default is 7.
#'
#' @return
#' \describe{
#'   \item{'Individual'}{Individual identifiers.}
#'   \item{'TrackedCycles'}{Tracked cycles.}
#'   \item{'NumTrue'}{Number of true values.}
#'   \item{'LogMean'}{Individual's mean of log(Y).}
#'   \item{'LogPrec'}{Individual's precision of log(Y)}
#'   \item{'Beta0'}{Beta0 true value.}
#'   \item{'Gamma0'}{Gamma0 true value.}
#'   \item{'X0',...,'XN'}{Covariate matrix for Mean, where N is the length of trueBetas.}
#'   \item{'Z0',...,'ZM'}{Covariate matrix for precision, where M is the length of trueGammas.}
#' }
#'
#' @seealso \code{\link{skipTrack.simulate}}
stSim <- function(i, skipProb, maxCycles, trueBetas, trueGammas, overlap, avgCyclesPer){
  #For each individual, generate the number of (tracked) cycles from poisson(7)
  #(restricted to > 0)
  numCycles <- max(rpois(1, avgCyclesPer), 1)

  #If trueBetas or trueGammas don't exist, set mean/precision to given average,
  #otherwise, create the number of requested covariates and record effects
  if(is.null(trueBetas)){
    lm <- log(30)
    xi <- NULL
  }else{
    xi <- matrix(rnorm(length(trueBetas), 0), nrow = 1)
    lm <- log(30) + xi %*% trueBetas
  }
  if(is.null(trueGammas)){
    precm <- exp(5.5)
    zi <- NULL
  }else{
    #Which x to overlap?
    if(overlap == 0){
      whichX <- 0
    }else{
      whichX <- 1:overlap
    }
    zi <- matrix(c(xi[1,whichX],
                   rnorm(length(trueGammas)-overlap, 0)), nrow = 1)
    precm <- exp(5.5 + zi %*% trueGammas)
  }

  #For each individual sample a mean (on the log scale) and precision (on the log scale)
  phi0 <- .01 #Constant goes here for phi
  prec <- max(4, rgamma(1, shape = precm*phi0, rate = phi0)) #Don't let precision get absurdly low
  lmean <- rnorm(1, lm, .13)

  #Sample c (true cycles per tracked cycle) values for number of cycles
  cs <- sample(1:maxCycles, numCycles, replace = TRUE, prob = skipProb)

  #Sample tracked cycle lengths
  ys <- round(rlnorm(numCycles, meanlog = lmean + log(cs), sdlog = sqrt(1/prec)))

  xi <- as.data.frame(cbind(1, xi))
  zi <- as.data.frame(cbind(1, zi))
  names(xi) <- paste0('X', 0:(ncol(xi)-1))
  names(zi) <- paste0('Z', 0:(ncol(zi)-1))

  #Return as data.frame
  df <- data.frame('Individual' = i, 'TrackedCycles' = ys, 'NumTrue' = cs,
                   'LogMean' = lmean, 'LogPrec' = prec, 'Beta0' = log(30), 'Gamma0' = 5.5)

  df <- cbind(df, xi, zi)
  return(df)
}

#' Simulate user tracked menstrual cycle data for an individual using the li model.
#'
#' This function generates synthetic data for user tracked menstrual cycles for a
#' single individual using the li model. For this model Beta0 = log(30), and Gamma0 doesn't really make sense.
#'
#' @param i Individual identifier. Character, numeric or integer.
#' @param skipProb Vector, ignored for this model.
#' @param maxCycles Integer, Maximum possible number of true cycles per tracked cycle.
#' @param trueBetas Optional. True values for generated mean regression coefficients.
#' @param trueGammas NULL, left for consistency. Will throw error if specified.
#' @param avgCyclesPer Average number of cycles contributed by each individual. Actual number is drawn from Poisson for each person. Default is 7.
#'
#' @return
#' \describe{
#'   \item{'Individual'}{Individual identifiers.}
#'   \item{'TrackedCycles'}{Tracked cycles.}
#'   \item{'NumTrue'}{Number of true values.}
#'   \item{'SkipProb'}{Individual's probability of skipping tracking a cycle}
#'   \item{'Mean'}{Individual's mean values.}
#'   \item{'Beta0'}{Beta0 true value.}
#'   \item{'Gamma0}{NA}
#'   \item{'Z0'}{1}
#'   \item{'X0',...,'XN'}{Covariate matrix for Mean, where N is the length of trueBetas.}
#' }
#' @seealso \code{\link{skipTrack.simulate}}
liSim <- function(i, skipProb, maxCycles, trueBetas, trueGammas = NULL, avgCyclesPer){
  if(!is.null(trueGammas)){
    warning('Li data generation model does not take covs for precision, trueGamma input will be ignored')
  }

  #For each individual, generate the number of (tracked) cycles from poisson(7)
  #(restricted to > 0)
  numCycles <- max(rpois(1, avgCyclesPer), 1)

  #For each individual sample a mean
  indMean <- rgamma(1, 180, 6)

  #For each individual sample a skip probability
  indSkip <- rbeta(1, 2, 20)

  #For each tracked cycle, sample the number of skipped cycles, cap at maxCycles
  numSkips <- rgeom(numCycles, 1-indSkip)
  numSkips <- pmin(numSkips, rep(maxCycles-1, numCycles))

  #If there are trueBeta values, adjust individual mean for those
  if(is.null(trueBetas)){
    chnge <- 0
    xi <- NULL
  }else{
    xi <- matrix(rnorm(length(trueBetas), 0), nrow = 1)
    chnge <- xi %*% trueBetas
  }

  #Adjust based on change (on log scale)
  indMean <- as.numeric(exp(log(indMean) + chnge))

  #For each tracked cycle, draw a length (dependent on numSkips)
  ys <- rpois(numCycles, indMean*(1+numSkips))

  #Return as data.frame
  #NOTE: Beta0 is log(30) here as skipTrack model assumes a distribution on log(Y) not just Y,
  df <- data.frame('Individual' = i, 'TrackedCycles' = ys, 'NumTrue' = numSkips + 1,
                   'Mean' = indMean, 'SkipProb' = indSkip,
                   'Z0' = 1, Beta0 = log(30), Gamma0 = NA)

  xi <- as.data.frame(cbind(1, xi))
  names(xi) <- paste0('X', 0:(ncol(xi)-1))


  df <- cbind(df, xi)

  return(df)
}

#' Simulate user tracked menstrual cycle data for an individual using the mixture model.
#'
#' This function generates synthetic data for user tracked menstrual cycles for a
#' single individual using the mixture model. For this model Beta_0 is set to log(30) and Gamma_0
#' is set to 15, although for the skipTrack model this lacks interpretation.
#'
#' @param i Individual identifier. Character, numeric, or integer.
#' @param skipProb Vector of probabilities for the number of true cycles per tracked cycle. For
#' example, (.7, .2, .1) means that 70% of observed cycles will contain one true cycle, 20%
#' will contain 2 true cycles, and 10% will contain 3 true cycles.
#' @param maxCycles Maximum number of true cycles per tracked cycle.
#' @param trueBetas Optional. True values for generated mean regression coefficients.
#' @param trueGammas Optional. True values for the generated precision regression coefficients.
#' @param overlap Optional. Number of (non-intercept) columns shared between X and Z. Columns are shared from left to right.
#' @param avgCyclesPer Average number of cycles contributed by each individual. Actual number is drawn from Poisson for each person. Default is 7.
#'
#' @return
#' \describe{
#'   \item{'Individual'}{Individual identifiers.}
#'   \item{'TrackedCycles'}{Tracked cycles.}
#'   \item{'NumTrue'}{Number of true values.}
#'   \item{'Mean'}{Individual's mean values.}
#'   \item{'Beta0'}{Beta0 true value.}
#'   \item{'Gamma0'}{Gamma0 true value.}
#'   \item{'X0',...,'XN'}{Covariate matrix for Mean, where N is the length of trueBetas.}
#'   \item{'Z0',...,'ZM'}{Covariate matrix for precision, where M is the length of trueGammas.}
#' }
#'
#' @seealso \code{\link{skipTrack.simulate}}
mixSim <- function(i, skipProb, maxCycles, trueBetas, trueGammas, overlap, avgCyclesPer){
  #For each individual, generate the number of (tracked) cycles from poisson(7)
  #(restricted to > 0)
  numCycles <- max(rpois(1, avgCyclesPer), 1)

  #Per cycle generate numTrue
  cs <- sample(1:maxCycles, numCycles, replace = TRUE, prob = skipProb)

  #If trueBetas or trueGammas == NULL, set mean/precision to parameters specifically,
  #otherwise, create the number of requested covariates and record effects
  if(is.null(trueBetas)){
    m <- 30
    xi <- NULL
  }else{
    xi <- matrix(rnorm(length(trueBetas), 0), nrow = 1)
    m <- 30 + xi %*% trueBetas
  }

  #Get individual mean based on trueBetas
  indMean <- round(rgamma(1, m*5, 5))

  #Create probabilities affected by Gammas which sort individuals into 3 regularity categories
  if(is.null(trueGammas)){
    p <- .5
    zi <- NULL
  }else{
    #Which x to overlap?
    if(overlap == 0){
      whichX <- 0
    }else{
      whichX <- 1:overlap
    }
    zi <- matrix(c(xi[1,whichX],
                   rnorm(length(trueGammas)-overlap, 0)), nrow = 1)
    linComp <- zi %*% trueGammas
    p <- exp(linComp)/(1+exp(linComp))
  }

  #Draw category given p
  indCat <- rbinom(1, 2, p)

  #Different behavior depending on category
  if(indCat == 0){
    ys <- rpois(numCycles, lambda = indMean*cs)
  }else if(indCat == 1){
    ys <- sapply(cs, function(drC){
      indMean*drC + sum(sample((-2:2), drC, replace = TRUE, prob = c(.05, .25,.4, .25, .05)))
    })
  }else if(indCat == 2){
    ys <- indMean*cs
  }

  #Create X Z dfs
  xi <- as.data.frame(cbind(1, xi))
  zi <- as.data.frame(cbind(1, zi))
  names(xi) <- paste0('X', 0:(ncol(xi)-1))
  names(zi) <- paste0('Z', 0:(ncol(zi)-1))


  df <- data.frame('Individual' = i, 'TrackedCycles' = ys, 'NumTrue' = cs,
                   'Mean' = indMean, Beta0 = log(30), Gamma0 = 15)

  #Attach covariate info
  df <- cbind(df, xi, zi)

  return(df)
}
