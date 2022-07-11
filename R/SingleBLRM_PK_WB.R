SingleBLRM.PK.WB <- function(){

  #  bivariate normal prior for (logAlpha,logBeta)
  #  Prior = (mean.logAlpha, mean.logBeta, sd.logAlpha, sd.logBeta, corr)

  #  PK-Toxicity model
  #  Define covariance matrices and model paramters for each mixture component
  for (j in 1:Nmix1) {
    cova[j,1,1] <- Prior[j,3]*Prior[j,3]
    cova[j,2,2] <- Prior[j,4]*Prior[j,4]
    cova[j,1,2] <- Prior[j,3]*Prior[j,4]*Prior[j,5]
    cova[j,2,1] <- cova[j,1,2]
    prec[j,1:2,1:2] <- inverse(cova[j,1:2,1:2])

    #  Prior distribution
    logAlphaBetaAll[j,1:2] ~ dmnorm(Prior[j,1:2],prec[j,1:2,1:2])
  }

  # latent mixture index
  Z ~ dcat(Pmix[1:Nmix1,1])
  for (j in 1:Nmix1) {
    Zcat[j] <- equals(Z,j)
  }
  logAlphaBeta[1] <-  logAlphaBetaAll[Z,1]
  logAlphaBeta[2] <-  logAlphaBetaAll[Z,2]

  logAlphaEBeta[1] <-  logAlphaBetaAll[Z,1]
  logAlphaEBeta[2] <-  exp(logAlphaBetaAll[Z,2])

  logAlphaBeta.XY <- logAlphaBeta[1]*logAlphaBeta[2]

  #Dose-PK model
  #  Define covariance matrices and model paramters for each mixture component
  for (j in 1:Nmix2){
    covd[j,1,1] <- Prior.PK[j,3]*Prior.PK[j,3]
    covd[j,2,2] <- Prior.PK[j,4]*Prior.PK[j,4]
    covd[j,1,2] <- Prior.PK[j,3]*Prior.PK[j,4]*Prior.PK[j,5]
    covd[j,2,1] <- covd[j,1,2]
    prec.PK[j,1:2,1:2] <- inverse(covd[j,1:2,1:2])

    #  Prior distribution
    logGammaAll[j,1:2] ~ dmnorm(Prior.PK[j,1:2],prec.PK[j,1:2,1:2])
  }

  # latent mixture index
  Z.PK ~ dcat(Pmix.PK[1:Nmix2,1])
  for (j in 1:Nmix2) {
    Zcat.PK[j] <- equals(Z.PK,j)
  }

  logGamma[1] <-  logGammaAll[Z,1]
  logGamma[2] <-  logGammaAll[Z,2]

  logEGamma[1] <-  logGammaAll[Z,1]
  logEGamma[2] <-  exp(logGammaAll[Z,2])

  logGamma.XY <- logGamma[1]*logGamma[2]


  #  define Puctoffs1, with 1 as last entry
  for (j in 1:(Ncat-1)) { Pcutoffs1[j] <- Pcutoffs[j] }
  Pcutoffs1[Ncat] <- 1

  #  Probabilities of Tox: pTox
  for(l in 1:Nout){
  for (i in 1:Ndoses)  {
    lin[i,l] <- logAlphaBeta[1] + (equals(monotone.DLT,1)*exp(logAlphaBeta[2])+ equals(monotone.DLT,2)*logAlphaBeta[2])*mu[i,l] #/(mean(mu[,])))
    logit(pTox[i,l]) <- lin[i,l]

    mu[i,l] <-  logGamma[1] + (equals(monotone.DE,1)*exp(logGamma[2])+ equals(monotone.DE,2)*logGamma[2])*log(Doses[i,1]/DoseRef) + inprod(Xout[l,1:Ncov],beta[1:Ncov,1])


    # for each dose, indicators for different toxicity categories
    # (means of these indicators correspond to category probabilities)
    pCat[i,1,l] <- step(Pcutoffs1[1] - pTox[i,l])
    for (j in 2:Ncat) {
      pCat[i,j,l] <- step(Pcutoffs1[j] - pTox[i,l]) - sum(pCat[i,1:(j-1),l])
    }

    #  Prediction likelihood
    Ntox.pred[i,l] ~ dbin(pTox[i,l],New.cohortsize)

    # Predictive probability of 0:New.cohortsize DLTs among New.cohortsize patients

    for (j in 1:(New.cohortsize+1)) {
      pred.Tox[i,j,l] <- equals(j-1,Ntox.pred[i,l])
    }

  }
  }

  for(i in 1:Ncov){

    prec.beta[i] <- pow(Prior.beta[i,2], -2)
    beta[i,1] ~ dnorm(Prior.beta[i,1], prec.beta[i])

  }

 #Sampling variability
  
  for(i in 1:Ndoses){
  prec.sd.PK[i] <- pow(Prior.sd.PK[i,2],-2)
  sd.PK[i] ~ dlnorm(Prior.sd.PK[i,1],prec.sd.PK[i])
  tau.PK[i] <- pow(sd.PK[i], -2)
  }
  
 #Random effect for PK model  

  prec.re.sd.PK[1] <- pow(Prior.re.sd.PK[1,2],-2)
  prec.re.sd.PK[2] <- pow(Prior.re.sd.PK[2,2],-2)
  
  sd.re.PK[1] ~ dlnorm( Prior.re.sd.PK[1,1],prec.re.sd.PK[1])
  sd.re.PK[2] ~ dlnorm( Prior.re.sd.PK[2,1],prec.re.sd.PK[2])
  rho.re.pk   ~ dunif(Prior.re.rho[1],Prior.re.rho[2])

  cov.re.PK[1,1]  <- pow(sd.re.PK[1],2)
  cov.re.PK[2,2]  <- pow(sd.re.PK[2],2)
  cov.re.PK[1,2]  <- sd.re.PK[1]*sd.re.PK[2]*rho.re.pk 
  cov.re.PK[2,1]  <- cov.re.PK[1,2]
  prec.re.PK[1:2,1:2] <- inverse(cov.re.PK[1:2,1:2])
  
  
  zero[1] <- 0
  zero[2] <- 0

  #  binomial likelihoods for DLT and Normal likelihood for PK parameter
  for (i in 1:Ncohorts)  {
    lin1[i] <- logAlphaBeta[1] + (equals(monotone.DLT,1)*exp(logAlphaBeta[2])+ equals(monotone.DLT,2)*logAlphaBeta[2])*mu1[i,1] #/mean(mu1[,]))
    logit(p1Tox[i]) <- lin1[i]

    re.PK[i, 1:2] ~ dmnorm(zero[1:2], prec.re.PK[1:2,1:2])

     mu1[i,1] <-  logGamma[1] + (equals(monotone.DE,1)*exp(logGamma[2])+ equals(monotone.DE,2)*logGamma[2])*log(DosesAdm[i,1]/DoseRef) + inprod(X[i,1:Ncov],beta[1:Ncov,1]) + equals(rand.m,1)*re.PK[i, 1] + equals(rand.slope,1)*equals(rand.m,1)*re.PK[i, 2]*log(DosesAdm[i,1]/DoseRef)  


    Ntox[i,1] ~ dbin(p1Tox[i],Npat[i,1])
    PK[i,1] ~ dnorm(mu1[i,1], tau.PK[Dose.Grp[i]])
  }
}


