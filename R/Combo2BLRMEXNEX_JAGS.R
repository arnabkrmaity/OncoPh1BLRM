## JAGS code for BLRM plus EXNEX function
## for doublet combination

combo2BLRMEXNEX<-function(){

  ## create a vector of two zeros for the dmnorm means under EX
    for (i in 1:2) {zero2[i] <- 0}
  
  ## Agent 1

  ##################################
  # prior covariance matrix for NEX
  cova1[1,1] <- Prior.NEX.Agent1[3]*Prior.NEX.Agent1[3]
  cova1[2,2] <- Prior.NEX.Agent1[4]*Prior.NEX.Agent1[4]
  cova1[1,2] <- Prior.NEX.Agent1[3]*Prior.NEX.Agent1[4]*Prior.NEX.Agent1[5]
  cova1[2,1] <- cova1[1,2]
  prec.NEX.Agent1[1:2,1:2] <- inverse(cova1[,]) ## convert to 1/variance for parameterization of dmnorm
  ##################################
  
  #####################################################
  # prior for EX
  prec.mu1.Agent1[1]<-pow(Prior.EX.mu1.Agent1[2],-2)
  prec.mu2.Agent1[1]<-pow(Prior.EX.mu2.Agent1[2],-2)

  mu1.Agent1[1] ~ dnorm(Prior.EX.mu1.Agent1[1],prec.mu1.Agent1[1]) ## prior for alpha
  mu2.Agent1[1] ~ dnorm(Prior.EX.mu2.Agent1[1],prec.mu2.Agent1[1]) ## prior for beta
  
  
    rho.Agent1  ~ dunif(Prior.rho.Agent1[1],Prior.rho.Agent1[2])
  
    prec.tau1.Agent1[1] <- pow(Prior.EX.tau1.Agent1[2],-2)
    prec.tau2.Agent1[1] <- pow(Prior.EX.tau2.Agent1[2],-2)
	
    tau.norm.Agent1[1] ~ dnorm( Prior.EX.tau1.Agent1[1],prec.tau1.Agent1[1])
	tau.norm.Agent1[2] ~ dnorm( Prior.EX.tau2.Agent1[1],prec.tau2.Agent1[1])
	
	tau.Agent1[1] <- max(exp(tau.norm.Agent1[1]),0.000001)## convert to log normal and set the max to a small value to avoid inversion errors
	tau.Agent1[2] <- max(exp(tau.norm.Agent1[2]),0.000001) ## convert to log normal
	

	

    cov.ex.Agent1[1,1]      <- pow(tau.Agent1[1],2)
    cov.ex.Agent1[2,2]      <- pow(tau.Agent1[2],2)
    cov.ex.Agent1[1,2]      <- tau.Agent1[1]*tau.Agent1[2]*rho.Agent1
    cov.ex.Agent1[2,1]      <- cov.ex.Agent1[1,2]
    prec.EX.tau.Agent1[1:2,1:2] <- inverse(cov.ex.Agent1[1:2,1:2])  
  ############################################################
    
  
  ## log-odds parameters under EX
  for(j in 1:Nstrata){
  
  reTmp.Agent1[j,1:2] ~ dmnorm(zero2[1:2],prec.EX.tau.Agent1[1:2,1:2]) ## Determine the random effects for log alpha and beta
re.Agent1[1,j]<-reTmp.Agent1[j,1] # log alpha random effects 
re.Agent1[2,j]<-reTmp.Agent1[j,2] # log beta random effects 

    log.alphabeta.EXNEX.Agent1[1,j] <- mu1.Agent1[1] + re.Agent1[1,j] ## update the prior log alpha with the random effects for strata j
	log.alphabeta.EXNEX.Agent1[2,j] <- mu2.Agent1[1] + re.Agent1[2,j] ## update the prior log beta with the random effects for strata j
	
  alpha.EXNEX.Agent1[1,j] <- exp(log.alphabeta.EXNEX.Agent1[1,j]) ## create a alpha vector with the EX result
  beta.EXNEX.Agent1[1,j] <- exp(log.alphabeta.EXNEX.Agent1[2,j])## create a beta vector with the EX result
  log.alpha.EXNEX.Agent1[1,j]<-log.alphabeta.EXNEX.Agent1[1,j]
  log.beta.EXNEX.Agent1[1,j]<-log.alphabeta.EXNEX.Agent1[2,j]
  }
  
  
    # log-odds parameters under NEX
  for (j in 1:Nstrata) {
  
      tmp.NEX.Agent1[j,1:2] ~ dmnorm(Prior.NEX.Agent1[1:2],prec.NEX.Agent1[1:2,1:2])
	  log.alphabeta.EXNEX.Agent1[3,j] <- tmp.NEX.Agent1[j,1]
	  log.alphabeta.EXNEX.Agent1[4,j] <- tmp.NEX.Agent1[j,2]
	  
  alpha.EXNEX.Agent1[2,j] <- exp(log.alphabeta.EXNEX.Agent1[3,j]) ## add alpha for NEX to alpha vector
  beta.EXNEX.Agent1[2,j] <- exp(log.alphabeta.EXNEX.Agent1[4,j]) ## add beta for NEX to beta vector
    log.alpha.EXNEX.Agent1[2,j]<-log.alphabeta.EXNEX.Agent1[3,j] 
  log.beta.EXNEX.Agent1[2,j]<-log.alphabeta.EXNEX.Agent1[4,j]
  
  }
  
  

  
  
  ###########
  ## Agent 2 
  ###########
  
  ##################################
  # prior covariance matrix for NEX
  cova2[1,1] <- Prior.NEX.Agent2[3]*Prior.NEX.Agent2[3]
  cova2[2,2] <- Prior.NEX.Agent2[4]*Prior.NEX.Agent2[4]
  cova2[1,2] <- Prior.NEX.Agent2[3]*Prior.NEX.Agent2[4]*Prior.NEX.Agent2[5]
  cova2[2,1] <- cova2[1,2]
  prec.NEX.Agent2[1:2,1:2] <- inverse(cova2[1:2,1:2]) ## convert to 1/variance for parameterization of dmnorm
  
  
  ##################################
  
  #####################################################
  # prior for EX
  prec.mu1.Agent2[1]<-pow(Prior.EX.mu1.Agent2[2],-2)
  prec.mu2.Agent2[1]<-pow(Prior.EX.mu2.Agent2[2],-2)

  mu1.Agent2[1] ~ dnorm(Prior.EX.mu1.Agent2[1],prec.mu1.Agent2[1]) ## prior for alpha
  mu2.Agent2[1] ~ dnorm(Prior.EX.mu2.Agent2[1],prec.mu2.Agent2[1]) ## prior for beta
  
  
    rho.Agent2  ~ dunif(Prior.rho.Agent2[1],Prior.rho.Agent2[2])
  
    prec.tau1.Agent2[1] <- pow(Prior.EX.tau1.Agent2[2],-2)
    prec.tau2.Agent2[1] <- pow(Prior.EX.tau2.Agent2[2],-2)
	
    tau.norm.Agent2[1] ~ dnorm( Prior.EX.tau1.Agent2[1],prec.tau1.Agent2[1])
	tau.norm.Agent2[2] ~ dnorm( Prior.EX.tau2.Agent2[1],prec.tau2.Agent2[1])
	
	tau.Agent2[1] <- max(exp(tau.norm.Agent2[1]),0.000001)## convert to log normal and set the max to a small value to avoid inversion errors
	tau.Agent2[2] <- max(exp(tau.norm.Agent2[2]),0.000001) ## convert to log normal
	

    cov.ex.Agent2[1,1]      <- pow(tau.Agent2[1],2)
    cov.ex.Agent2[2,2]      <- pow(tau.Agent2[2],2)
    cov.ex.Agent2[1,2]      <- tau.Agent2[1]*tau.Agent2[2]*rho.Agent2
    cov.ex.Agent2[2,1]      <- cov.ex.Agent2[1,2]
    prec.EX.tau.Agent2[1:2,1:2] <- inverse(cov.ex.Agent2[1:2,1:2])  
  ############################################################
    
  
  ## log-odds parameters under EX
  for(j in 1:Nstrata){
  
  reTmp.Agent2[j,1:2] ~ dmnorm(zero2[1:2],prec.EX.tau.Agent2[1:2,1:2]) ## Determine the random effects for log alpha and beta
re.Agent2[1,j]<-reTmp.Agent2[j,1] # log alpha random effects 
re.Agent2[2,j]<-reTmp.Agent2[j,2] # log beta random effects 

    log.alphabeta.EXNEX.Agent2[1,j] <- mu1.Agent2[1] + re.Agent2[1,j] ## update the prior log alpha with the random effects for strata j
	log.alphabeta.EXNEX.Agent2[2,j] <- mu2.Agent2[1] + re.Agent2[2,j] ## update the prior log beta with the random effects for strata j
	
  alpha.EXNEX.Agent2[1,j] <- exp(log.alphabeta.EXNEX.Agent2[1,j]) ## create a alpha vector with the EX result
  beta.EXNEX.Agent2[1,j] <- exp(log.alphabeta.EXNEX.Agent2[2,j])## create a beta vector with the EX result
  log.alpha.EXNEX.Agent2[1,j]<-log.alphabeta.EXNEX.Agent2[1,j]
  log.beta.EXNEX.Agent2[1,j]<-log.alphabeta.EXNEX.Agent2[2,j]
  }
  
  
  
    # log-odds parameters under NEX
  for (j in 1:Nstrata) {
  
      tmp.NEX.Agent2[j,1:2] ~ dmnorm(Prior.NEX.Agent2[1:2],prec.NEX.Agent2[1:2,1:2])
	  log.alphabeta.EXNEX.Agent2[3,j] <- tmp.NEX.Agent2[j,1]
	  log.alphabeta.EXNEX.Agent2[4,j] <- tmp.NEX.Agent2[j,2]
	  
  alpha.EXNEX.Agent2[2,j] <- exp(log.alphabeta.EXNEX.Agent2[3,j]) ## add alpha for NEX to alpha vector
  beta.EXNEX.Agent2[2,j] <- exp(log.alphabeta.EXNEX.Agent2[4,j]) ## add beta for NEX to beta vector
    log.alpha.EXNEX.Agent2[2,j]<-log.alphabeta.EXNEX.Agent2[3,j] 
  log.beta.EXNEX.Agent2[2,j]<-log.alphabeta.EXNEX.Agent2[4,j]
  
  }
  
  
    
  
  ##################################
  ## Interaction Parameter
  ##################################
  
  
  ## eta Under EX 
    etaPriorPrec.mu[1] <- pow(Prior.EX.eta.mu[2], -2)
    mu.eta[1] ~ dnorm(Prior.EX.eta.mu[1],etaPriorPrec.mu[1]) ## prior for eta
	
	## prior for tau.eta (heterogeneity in the interactions among the studies)
    etaPriorPrec.tau[1] <- pow(Prior.EX.eta.tau[2], -2)	
	
    tau.norm.eta[1] ~ dnorm(Prior.EX.eta.tau[1],etaPriorPrec.tau[1]) 
	tau.eta<-max(exp(tau.norm.eta[1]),0.000001)## convert to log normal
	etaPriorPrec[1] <- pow(tau.eta, -2)

	
  ## log-odds parameters under EX
    for(j in 1:Nstrata){
	
	reTmp.eta[j,1] ~ dnorm(zero2[1],etaPriorPrec[1]) ## Determine the random effects for the interaction
	  
	EXNEX.eta[1,j] <- mu.eta[1] + reTmp.eta[j,1] ## update the interaction prior with the random effects for strata j

  }
  
  
  ## eta Under NEX
    etaPriorPrec.NEX[1] <- pow(Prior.NEX.eta[2], -2)
	
    for (j in 1:Nstrata) {
   ## interaction under NEX

    eta.NEX[1,j] ~ dnorm(Prior.NEX.eta[1],etaPriorPrec.NEX[1]) ## prior for alpha
	
	EXNEX.eta[2,j]<-eta.NEX[1,j]	
		
  }
  


 #################################
 ## Latent Mixture indicators
 #################################
 
 
  ## choose an alpha, beta, and eta 
  
    # latent mixture indicators:
  # exch.index: categorial 1,...,Nmix=Nexch+1
  ## create a categorial distribution of 1:Nmix indicating which EX distribution or NEX distribution to use with probability pMix for the category. E.g., Nmix = 2 and pMix = c(0.9,0.1) then exch.index = 1 90% of the time and 2 10%. This is then used to select the alpha and beta parameters that corespond to the EX or NEX distributions
  for (j in 1:Nstrata) {
    exch.index[j] ~ dcat(pMix[1:Nmix])
      
   
      ## Indicator to determine whether EX or NEX is used
	  ## 1 = EX and 2 = NEX
    for (jj in 1:Nmix) {
      exch[j,jj] <- equals(exch.index[j],jj) ## logical when the exch = EX or EXNEX
    }
	
    }
  
  # pick alpha and beta for Agent 1 and Agent 2 and an eta
  for (j in 1:Nstrata) {
  
  	log.alpha.Agent1[j]<-log.alpha.EXNEX.Agent1[exch.index[j],j]
	log.beta.Agent1[j]<-log.beta.EXNEX.Agent1[exch.index[j],j]
	alpha.Agent1[j]<-alpha.EXNEX.Agent1[exch.index[j],j]
	beta.Agent1[j]<-beta.EXNEX.Agent1[exch.index[j],j]

	log.alpha.Agent2[j]<-log.alpha.EXNEX.Agent2[exch.index[j],j]
	log.beta.Agent2[j]<-log.beta.EXNEX.Agent2[exch.index[j],j]
	alpha.Agent2[j]<-alpha.EXNEX.Agent2[exch.index[j],j]
	beta.Agent2[j]<-beta.EXNEX.Agent2[exch.index[j],j]
	
	eta[j]<-EXNEX.eta[exch.index[j],j]
	OddsFactor[j] <- exp(eta[j])
	
  }
  
  

   
  

  # revised data input
 #binomial likelihoods (data) for combo
  
  for(i in 1:Nstrata){ ## select strata/trial data with combination therapy
  
  for (j in 1:Ncohorts[i]){ ## loop over the dose levels tested for a given strata/trial
  
 
    lin1.2[i,j] <- log.alpha.Agent1[i] + beta.Agent1[i]*log(DosesAdm1[i,j]/DoseRef1) 
	    lin1[i,j] <- lin1.2[i,j] 
		#+ step(-20 - lin1.2[i,j]) * (-20 - lin1.2[i,j]) +  step(lin1.2[i,j] - 20) * (20 - lin1.2[i,j]) # if x >= 20 or <= -20 rescale to 20 or -20
    odds11[i,j] <- exp(lin1[i,j]) 

    
    lin2.2[i,j] <- log.alpha.Agent2[i] + beta.Agent2[i]*log(DosesAdm2[i,j]/DoseRef2) 
    lin2[i,j] <- lin2.2[i,j] 
	#+ step(-20 - lin2.2[i,j]) * (-20 - lin2.2[i,j]) +  step(lin2.2[i,j] - 20) * (20 - lin2.2[i,j]) # if x >= 20 or <= -20 rescale to 20 or -20
    odds22[i,j] <- exp(lin2[i,j])
    
    odds121.0[i,j] <- odds11[i,j] + odds22[i,j] + odds11[i,j]*odds22[i,j]         
    odds121[i,j] <-  odds121.0[i,j]*exp(eta[i]*(DosesAdm1[i,j]/DoseRef1)*(DosesAdm2[i,j]/DoseRef2))
    Pr.Tox1[i,j] <- odds121[i,j]/(1 + odds121[i,j])
	
	Ntox[i,j] ~ dbin(Pr.Tox1[i,j], Npat[i,j]) 
    
    
  } 
  
}



  
## Marginal odds for each Strata/Agent

  for (j in 1:Nstrata) {
  for (i in 1:Ndoses1){
    #  Marginal odds for Agent 1
    lin1.1[j,i] <- log.alpha.Agent1[j] + beta.Agent1[j]*log(Doses1[i]/DoseRef1)
    odds1[j,i] <- exp(lin1.1[j,i]) 
  }
  }
  
   for (j in 1:Nstrata) {
  for(i in 1:Ndoses2){      
    #  Marginal odds for Agent 2 
    lin2.1[j,i] <- log.alpha.Agent2[j] + beta.Agent2[j]*log(Doses2[i]/DoseRef2)
    odds2[j,i] <- exp(lin2.1[j,i])
  } 
  }
  

  ## will return strata specific results for the provisional dose levels. Assumes the same provisional dose levels are used across all studies 
 for (k in 1:Nstrata) { 
  #  combination model and probability of each category for each combination
  for(j in 1:Ndoses1){ ## number of provisional doses for Agent 1
    for(i in 1:Ndoses2){ ## number of provisional doses for Agent 2 
      #  Odds without and with interaction 
      odds12.0[k,j,i] <- odds1[k,j] + odds2[k,i] + odds1[k,j]*odds2[k,i]
      odds12[k,j,i] <- odds12.0[k,j,i]*exp(eta[k]*(Doses1[j]/DoseRef1)*(Doses2[i]/DoseRef2))
      Pr.Tox[k,j,i] <- odds12[k,j,i]/(1 + odds12[k,j,i]) ## probablity of DLT at dose i, j
      	  
	  
	Pr.Cat[k,j,i,1] <- step(Pint[1]-Pr.Tox[k,j,i]) ## indicator for underdosing (mean is the prob of the interval)
    Pr.Cat[k,j,i,2] <- step(Pint[2]-Pr.Tox[k,j,i])-step(Pint[1]-Pr.Tox[k,j,i]) ##  indicator for target tox (mean is the prob of the interval)
    Pr.Cat[k,j,i,3] <- step(1-Pr.Tox[k,j,i])-step(Pint[2]-Pr.Tox[k,j,i]) ##  indicator for overdose (mean is the prob of the interval)

	  
     }
    }
	}
  
    


 
 
  
  
  
}
  