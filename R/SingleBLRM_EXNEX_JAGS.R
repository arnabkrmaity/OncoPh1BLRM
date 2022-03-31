## JAGS code for BLRM plus EXNEX function

## single agent model

SingleBLRMEXNEX<-function(){
  
  ## create a vector of two zeros for the dmnorm means under EX
  for (i in 1:2) {zero2[i] <- 0}
  
  ##################################
  # prior covariance matrix for NEX
  cova[1,1] <- Prior.NEX[3]*Prior.NEX[3]
  cova[2,2] <- Prior.NEX[4]*Prior.NEX[4]
  cova[1,2] <- Prior.NEX[3]*Prior.NEX[4]*Prior.NEX[5]
  cova[2,1] <- cova[1,2]
  prec.NEX[1:2,1:2] <- inverse(cova[1:2,1:2]) ## convert to 1/variance for parameterization of dmnorm
  ##################################
  
  #####################################################
  # prior for EX
  prec.mu1[1]<-pow(Prior.EX.mu1[2],-2)
  prec.mu2[1]<-pow(Prior.EX.mu2[2],-2)
  
  mu1[1] ~ dnorm(Prior.EX.mu1[1],prec.mu1[1]) 
  mu2[1] ~ dnorm(Prior.EX.mu2[1],prec.mu2[1])
  
  
  rho  ~ dunif(Prior.rho[1],Prior.rho[2])
  
  prec.tau1[1] <- pow(Prior.EX.tau1[2],-2)
  prec.tau2[1] <- pow(Prior.EX.tau2[2],-2)
  
  
  tau.norm[1] ~ dnorm(Prior.EX.tau1[1],prec.tau1[1])# use a normal distribution and convert to log normal as it is more stable.
  tau.norm[2] ~ dnorm(Prior.EX.tau2[1],prec.tau2[1])
  
  tau[1] <- max(exp(tau.norm[1]),0.000001) ## convert to log normal and set the max to a small value to avoid inversion errors
  tau[2] <- max(exp(tau.norm[2]),0.000001) ## convert to log normal
  
  

  cov.ex[1,1]      <- tau[1]*tau[1]
  cov.ex[2,2]      <- tau[2]*tau[2]
  cov.ex[1,2]      <- tau[1]*tau[2]*rho
  cov.ex[2,1]      <- cov.ex[1,2]
  prec.EX.tau[1:2,1:2] <- inverse(cov.ex[1:2,1:2])  
  ############################################################
  
  
  
  
  
  
  ## log-odds parameters under EX
  for(j in 1:Nstrata){
    
    reTmp[j,1:2] ~ dmnorm(zero2[1:2],prec.EX.tau[1:2,1:2]) ## Determine the random effects for log alpha and beta
    re[1,j]<-reTmp[j,1] # log alpha random effects 
    re[2,j]<-reTmp[j,2] # log beta random effects 
    
    log.alphabeta.EXNEX[1,j] <- mu1[1] + re[1,j] ## update the prior log alpha with the random effects for strata j
    log.alphabeta.EXNEX[2,j] <- mu2[1] + re[2,j] ## update the prior log beta with the random effects for strata j
    
    alpha.EXNEX[1,j] <- exp(log.alphabeta.EXNEX[1,j]) ## create a alpha vector with the EX result
    beta.EXNEX[1,j] <- exp(log.alphabeta.EXNEX[2,j])## create a beta vector with the EX result
    log.alpha.EXNEX[1,j]<-log.alphabeta.EXNEX[1,j]
    log.beta.EXNEX[1,j]<-log.alphabeta.EXNEX[2,j]
  }
  
  
  # log-odds parameters under NEX
  for (j in 1:Nstrata) {
    
    tmp.NEX[j,1:2] ~ dmnorm(Prior.NEX[1:2],prec.NEX[1:2,1:2])
    log.alphabeta.EXNEX[3,j] <- tmp.NEX[j,1]
    log.alphabeta.EXNEX[4,j] <- tmp.NEX[j,2]
    
    alpha.EXNEX[2,j] <- exp(log.alphabeta.EXNEX[3,j]) ## add alpha for NEX to alpha vector
    beta.EXNEX[2,j] <- exp(log.alphabeta.EXNEX[4,j]) ## add beta for NEX to beta vector
    log.alpha.EXNEX[2,j]<-log.alphabeta.EXNEX[3,j] 
    log.beta.EXNEX[2,j]<-log.alphabeta.EXNEX[4,j]
    
  }
  
  
  ## choose an alpha and beta
  
  # latent mixture indicators:
  # exch.index: categorial 1,...,Nmix=Nexch+1
  # exch: Nstrata x Nmix matrix of 0/1 elements
  for (j in 1:Nstrata) {
    exch.index[j] ~ dcat(pMix[1:Nmix])## create a categorial distribution of 1:Nmix indicating which EX distribution or NEX distribution to use with probability pMix for the category. E.g., Nmix = 2 and pMix = c(0.9,0.1) then exch.index = 1 90% of the time and 2 10%. This is then used to select the alpha and beta parameters that corespond to the EX or NEX distributions
    
    ## Indicator to determine which mixture probability was used
    for (jj in 1:Nmix) {
      exch[j,jj] <- equals(exch.index[j],jj) ## logical when the exch = EX or EXNEX
    }
    
  }
  
  # pick alpha and beta
  for (j in 1:Nstrata) {
    
    log.alpha[j]<-log.alpha.EXNEX[exch.index[j],j]
    log.beta[j]<-log.beta.EXNEX[exch.index[j],j]
    alpha[j]<-alpha.EXNEX[exch.index[j],j]
    beta[j]<-beta.EXNEX[exch.index[j],j]
    
  }
  
  
  
  # Likelihood 
  for(j in 1:Nstrata){
    
    for (i in 1:Ncohorts[j]){
      
      Ntox[j,i] ~ dbin(Pr.Tox1[j,i],Npat[j,i]) ## number of toxicities at dose level i are distributed by a binomal 
      logit(Pr.Tox1[j,i])<- log.alpha[j]+beta[j]*log(DosesAdm[j,i]/DoseRef) ## Logistic regression model for the probability of a DLT
      
    }
    
  }
  
  ## will return strata specific results for the provisional dose levels. Assumes the same provisional dose levels are used across all studies 
  for (j in 1:Nstrata) { 
    
    # for each dose: probabilities of toxicity, category probabilities, risks
    for (i in 1:Ndoses) { ## for each provisional dose
      lin[j,i]<-log.alpha[j]+beta[j]*log(Doses[i]/DoseRef) ## log odds of dlt across all provisional dose levels
      logit(Pr.Tox[j,i]) <- lin[j,i] ## probability of a dlt for a given strata and at a given dose level
      Pr.Cat[j,i,1] <- step(Pint[1]-Pr.Tox[j,i]) ## prob of underdosing
      Pr.Cat[j,i,2] <- step(Pint[2]-Pr.Tox[j,i])-step(Pint[1]-Pr.Tox[j,i]) ## prob of target
      Pr.Cat[j,i,3] <- step(1-Pr.Tox[j,i])-step(Pint[2]-Pr.Tox[j,i]) ## prob of unacceptable tox
      ProvDoseLelCat[j,i]<-step(-Pr.Tox[j,i])+Doses[i] ## display what dose was used. Used to see what happens when strata have different number of dose levels.
    }
    
  }
  
  }
  
  


