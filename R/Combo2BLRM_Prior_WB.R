Combo2BLRM.Prior.WB <- function(){
  
  #  bivariate normal prior for (logAlpha1,logBeta1) [marginal model for Agent1]
  #  Prior1 = (mean.logAlpha1, mean.logBeta1, sd.logAlpha1, sd.logBeta1, corr1)
  
  #  Covariance and precision matrix for agent 1 for each mixture component
  
  for(j in 1:Nmix1){
    cova1[j,1,1] <- Prior1[j,3]*Prior1[j,3]
    cova1[j,2,2] <- Prior1[j,4]*Prior1[j,4]
    cova1[j,1,2] <- Prior1[j,3]*Prior1[j,4]*Prior1[j,5]
    cova1[j,2,1] <- cova1[j,1,2]
    prec1[j,1:2,1:2] <- inverse(cova1[j,,])
    
    #  Prior distribution for logAlpha1,logBeta1 for each mixture component
    logAlphaBeta1All[j,1:2] ~ dmnorm(Prior1[j,1:2],prec1[j,1:2,1:2])
  }
  
  # latent mixture index
  Z1 ~ dcat(Pmix1[1:Nmix1,1]) 
  
  for (j in 1:Nmix1) {
    Z1cat[j,1] <- equals(Z1,j)
  }
  
  logAlphaBeta1[1] <-  logAlphaBeta1All[Z1,1]
  logAlphaBeta1[2] <-  logAlphaBeta1All[Z1,2]
  
  logAlphaEBeta1[1] <-  logAlphaBeta1All[Z1,1]
  logAlphaEBeta1[2] <-  exp(logAlphaBeta1All[Z1,2])
  
  
  logAlphaBeta1.XY <- logAlphaBeta1[1]*logAlphaBeta1[2]
  
  
  #  bivariate normal prior for (logAlpha2,logBeta2) [marginal model for Agent2]
  #  Prior2 = (mean.logAlpha2, mean.logBeta2, sd.logAlpha2, sd.logBeta2, corr2)
  
  #  Covariance and precision matrix for agent 2 for each mixture component
  
  for(j in 1:Nmix2){
    cova2[j,1,1] <- Prior2[j,3]*Prior2[j,3]
    cova2[j,2,2] <- Prior2[j,4]*Prior2[j,4]
    cova2[j,1,2] <- Prior2[j,3]*Prior2[j,4]*Prior2[j,5]
    cova2[j,2,1] <- cova2[j,1,2]
    prec2[j,1:2,1:2] <- inverse(cova2[j,,])
    
    #  Prior distribution for logAlpha2,logBeta2 for each mixture component
    logAlphaBeta2All[j,1:2] ~ dmnorm(Prior2[j,1:2],prec2[j,1:2,1:2])
  }
  
  # latent mixture index
  Z2 ~ dcat(Pmix2[1:Nmix2,1]) 
  
  for (j in 1:Nmix2){
    Z2cat[j,1] <- equals(Z2,j)
  }
  
  logAlphaBeta2[1] <-  logAlphaBeta2All[Z2,1]
  logAlphaBeta2[2] <-  logAlphaBeta2All[Z2,2]
  
  logAlphaEBeta2[1] <-  logAlphaBeta2All[Z2,1]
  logAlphaEBeta2[2] <-  exp(logAlphaBeta2All[Z2,2])
  
  logAlphaBeta2.XY <- logAlphaBeta2[1]*logAlphaBeta2[2]
  
  
  # Normal prior for interaction parameter (eta) 
  
  for(i in 1:Nmix.eta){
  etaPriorPrec[i] <- pow(etaPrior[i,2], -2)
  eta12[i]  ~ dnorm(etaPrior[i,1], etaPriorPrec[i])
  }
  
  # latent mixture index
  Z.eta ~ dcat(P.eta.mix[1:Nmix.eta,1]) 
  
  for (j in 1:Nmix.eta){
    Z.etacat[j,1] <- equals(Z.eta,j)
  }
  
  eta <-  eta12[Z.eta]
  OddsFactor <- exp(eta) 
  
  #  define Puctoffs1, with 1 as last entry
  for (j in 1:(Ncat-1)) { Pcutoffs1[j] <- Pcutoffs[j] }
  Pcutoffs1[Ncat] <- 1
  
  for (i in 1:Ndoses1){
    #  Marginal odds for Agent 1
    lin1.1[i] <- logAlphaBeta1[1] + exp(logAlphaBeta1[2])*log(Doses1[i,1]/DoseRef1)
    odds1[i] <- exp(lin1.1[i]) 
  }
  for(i in 1:Ndoses2){      
    #  Marginal odds for Agent 2 
    lin2.1[i] <- logAlphaBeta2[1] + exp(logAlphaBeta2[2])*log(Doses2[i,1]/DoseRef2)
    odds2[i] <- exp(lin2.1[i])
  } 
  
  #  combination model and probability of each category for each combination
  for(i in 1:Ndoses1){
    for(j in 1:Ndoses2){
        #  Odds without and with interaction 
        odds12.0[i,j] <- odds1[i] + odds2[j] + odds1[i] * odds2[j]
        odds12[i,j] <- odds12.0[i,j]*exp(eta*(Doses1[i,1]/DoseRef1)*(Doses2[j,1]/DoseRef2))
        P12[i,j] <- odds12[i,j]/(1 + odds12[i,j])
        
        #  Prediction likelihood
        Ntox.pred[i,j] ~ dbin(P12[i,j],New.cohortsize)
        
        
        # for each dose, indicators for different toxicity categories
        # (means of these indicators correspond to category probabilities)   
        pCat[i,j,1] <- step(Pcutoffs1[1] - P12[i,j])
        for (l in 2:Ncat) {
          pCat[i,j,l] <- step(Pcutoffs1[l] - P12[i,j]) - sum(pCat[i,j,1:(l-1)])
        }
        
        #  Predictive probability of 0:New.cohortsize DLTs among New.cohortsize patients
        
        for (ll in 1:(New.cohortsize+1)){
          pred.Tox[i,j,ll] <- equals(ll-1,Ntox.pred[i,j])
        }
        
      }
    }
  
 
}
