SingleBLRM.Weight.WB <- function(){

  #  bivariate normal prior for (logAlpha,logBeta)
  #  Prior = (mean.logAlpha, mean.logBeta, sd.logAlpha, sd.logBeta, corr)

  #  Define covariance matrices and model paramters for each mixture component
  for (j in 1:Nmix) {
    cova[j,1,1] <- Prior[j,3]*Prior[j,3]
    cova[j,2,2] <- Prior[j,4]*Prior[j,4]
    cova[j,1,2] <- Prior[j,3]*Prior[j,4]*Prior[j,5]
    cova[j,2,1] <- cova[j,1,2]
    prec[j,1:2,1:2] <- inverse(cova[j,1:2,1:2])

    #  Prior distribution
    logAlphaBetaAll[j,1:2] ~ dmnorm(Prior[j,1:2],prec[j,1:2,1:2])
  }

  # latent mixture index
  Z ~ dcat(Pmix[1:Nmix,1])
  for (j in 1:Nmix) {
    Zcat[j] <- equals(Z,j)
  }
  logAlphaBeta[1] <-  logAlphaBetaAll[Z,1]
  logAlphaBeta[2] <-  logAlphaBetaAll[Z,2]

  logAlphaEBeta[1] <-  logAlphaBetaAll[Z,1]
  logAlphaEBeta[2] <-  exp(logAlphaBetaAll[Z,2])

  logAlphaBeta.XY <- logAlphaBeta[1]*logAlphaBeta[2]



  #  define Puctoffs1, with 1 as last entry
  for (j in 1:(Ncat-1)) { Pcutoffs1[j] <- Pcutoffs[j] }
  Pcutoffs1[Ncat] <- 1

  #  Probabilities of Tox: pTox
  for (i in 1:Ndoses)  {
    lin[i] <- logAlphaBeta[1] + exp(logAlphaBeta[2])*log(Doses[i,1]/DoseRef)
    logit(pTox[i]) <- lin[i]


    # for each dose, indicators for different toxicity categories
    # (means of these indicators correspond to category probabilities)
    pCat[i,1] <- step(Pcutoffs1[1] - pTox[i])
    for (j in 2:Ncat) {
      pCat[i,j] <- step(Pcutoffs1[j] - pTox[i]) - sum(pCat[i,1:(j-1)])
    }

    #  Prediction likelihood
    Ntox.pred[i] ~ dbin(pTox[i],New.cohortsize)

    # Predictive probability of 0:New.cohortsize DLTs among New.cohortsize patients

    for (j in 1:(New.cohortsize+1)) {
      pred.Tox[i,j] <- equals(j-1,Ntox.pred[i])
    }

  }

  #  binomial likelihoods
  for (i in 1:Ncohorts)  {
    lin1[i] <- logAlphaBeta[1] + exp(logAlphaBeta[2])*log(DosesAdm[i,1]/DoseRef)
    logit(p1Tox[i]) <- lin1[i]

    #Using zero trick for the likelihood for any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)

    #zero[i] <- 0
    minus.loglik[i] <- -Ntox[i,1]*log(p1Tox[i]) - (Npat[i,1] - Ntox[i,1])*log(1 - p1Tox[i])
    zero[i] ~ dpois(minus.loglik[i])
  }
}


