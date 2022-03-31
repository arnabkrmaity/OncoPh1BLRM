MAP.BLRM.WB <- function(){
  # mean mu.ex of exchangeability distribution
  # Priors: Prior.muA, Prior.muB
  muA.prec <- pow( Prior.muA[2], -2 )
  muB.prec <- pow( Prior.muB[2], -2 )
  
  mu.unconst[1] ~ dnorm( Prior.muA[1],muA.prec)
  mu.unconst[2] ~ dnorm( Prior.muB[1],muB.prec)
  # for numerical stability:
  # constrain to -10 and +10 (muA), and -5 and +5 (mu-logBeta)
  mu.ex[1] <- mu.unconst[1] + step(-10-mu.unconst[1])*(-10-mu.unconst[1])
  +step(mu.unconst[1]-10)*(10-mu.unconst[1])
  
  mu.ex[2] <- mu.unconst[2] + step(-5-mu.unconst[2])*(-5-mu.unconst[2]) +
    step(mu.unconst[2]-5)*(5-mu.unconst[2])
  
  # covariance matrix cov.ex of exchangeability distribution; 
  # parameters tau1,tau2,rho    
  # Priors: Prior.tauA, Prior.tauB, Prior.rho       
  
  # Assuming *Nexch* exchangeability distribution
  # with same mu.ex, same rho, but different tau (differential discounting)
  # tau.index: vector of length (Nstudy+1) representing exch distribution 
  rho    ~ dunif(Prior.rho[1],Prior.rho[2])
  
  for (j in 1:Nexch) {
    tau[j,1] ~ dlnorm( Prior.tauA[j,1],prec.tauA[j])
    tau[j,2] ~ dlnorm( Prior.tauB[j,1],prec.tauB[j])
    prec.tauA[j] <- pow(Prior.tauA[j,2],-2)
    prec.tauB[j] <- pow(Prior.tauB[j,2],-2)
    cov.ex[j,1,1]      <- pow(tau[j,1],2)
    cov.ex[j,2,2]      <- pow(tau[j,2],2)
    cov.ex[j,1,2]      <- tau[j,1]*tau[j,2]*rho
    cov.ex[j,2,1]      <- cov.ex[j,1,2]
    prec.ex[j,1:2,1:2] <- inverse(cov.ex[j,1:2,1:2])
  }
  
  # logAB = (logAlpha, logBeta):  parameters for each trial
  # *Nstudy+1* is the MAP prediction
  # random effects for all strata and prediction (in position *Nstudy+1*)
  zero[1] <- 0
  zero[2] <- 0
  for (j in 1:(Nstudy+1)) {
    re[j,1:2] ~ dmnorm(zero[1:2],prec.ex[tau.index[j,1],1:2,1:2])
    logAB[j,1] <-  mu.ex[1] + re[j,1]
    logAB[j,2] <-  mu.ex[2] + re[j,2]
    
    logAB.XY[j] <- logAB[j,1]*logAB[j,2]
  }
  #Parameters from MAP distribution
  MAP.logAB[1] <- logAB[Nstudy+1,1]
  MAP.logAB[2] <- logAB[Nstudy+1,2]
  
  # get E(logA x logB), to get BVN correlation if MCMC
  
  MAP.logAB.XY <- MAP.logAB[1]*MAP.logAB[2]
  
  
  # likelihood/sampling model 
  for (i in 1:Nobs) {
    lin1[i] <- logAB[Study[i,1], 1] + exp(logAB[Study[i,1], 2])*log(DosesAdm[i,1]/DoseRef)
    logit(p1Tox[i]) <- lin1[i]
    Ntox[i,1] ~ dbin(p1Tox[i],Npat[i,1])                            
  }
  
  #  define Puctoffs1, with 1 as last entry
  for (j in 1:(Ncat-1)) { 
    Pcutoffs1[j] <- Pcutoffs[j] 
  }
  Pcutoffs1[Ncat] <- 1
  
  # for all studies and doses
  #   P: (Nstudy+1) x Ndoses 
  #   Pcat: (Nstudy+1) x Ndoses x  Ncat
  
  for (j in 1:(Nstudy+1)){
    for (i in 1:Ndoses){ 
      lin[j,i] <- logAB[j, 1] + exp(logAB[j, 2])*log(Doses[i,1]/DoseRef)
      logit(pTox[j,i]) <- lin[j,i]
      
      # for each dose, indicators for different toxicity categories
      # (means of these indicators correspond to category probabilities)   
      pCat[j,i,1] <- step(Pcutoffs1[1] - pTox[j,i])
      for (k in 2:Ncat) {
        pCat[j,i,k] <- step(Pcutoffs1[k] - pTox[j,i]) - sum(pCat[j,i,1:(k-1)])
      }
    }
  }
  
  # pTox for MAP
  for (i in 1:Ndoses){ 
    MAP.pTox[i] <- pTox[Nstudy+1,i]
  }
  
  # pCat for MAP
  for (i in 1:Ndoses){
    for (k in 1:Ncat) {
      MAP.pCat[i,k] <- pCat[Nstudy+1,i,k]
    }
  }
  
 for (i in 1:Ndoses)  {
  #  Prediction likelihood
  Ntox.pred[i] ~ dbin(pTox[Nstudy+1,i],New.cohortsize)
  
  # Predictive probability of 0:New.cohortsize DLTs among New.cohortsize patients
  
  for (j in 1:(New.cohortsize+1)) {
    pred.Tox[i,j] <- equals(j-1,Ntox.pred[i])
                     }
  
                      }
}


