Combo3BLRM.WB <- function(){
  
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
  
  
  #  bivariate normal prior for (logAlpha3,logBeta3) [marginal model for Agent3]
  #  Prior3 = (mean.logAlpha3, mean.logBeta3, sd.logAlpha3, sd.logBeta3, corr3)
  
  #  Covariance and precision matrix for agent 3 for each mixture component
  
  for(j in 1:Nmix3){ 
  cova3[j,1,1] <- Prior3[j,3]*Prior3[j,3]
  cova3[j,2,2] <- Prior3[j,4]*Prior3[j,4]
  cova3[j,1,2] <- Prior3[j,3]*Prior3[j,4]*Prior3[j,5]
  cova3[j,2,1] <- cova3[j,1,2]
  prec3[j,1:2,1:2] <- inverse(cova3[j,,])
  
  #  Prior distribution for logAlpha2,logBeta2 for each mixture component
  logAlphaBeta3All[j,1:2] ~ dmnorm(Prior3[j,1:2],prec3[j,1:2,1:2])
  
  }
  
  # latent mixture index
  Z3 ~ dcat(Pmix3[1:Nmix3,1]) 
  
  for (j in 1:Nmix3){
    Z3cat[j,1] <- equals(Z3,j)
  }
  
  logAlphaBeta3[1] <-  logAlphaBeta3All[Z3,1]
  logAlphaBeta3[2] <-  logAlphaBeta3All[Z3,2]
  
  logAlphaEBeta3[1] <-  logAlphaBeta3All[Z3,1]
  logAlphaEBeta3[2] <-  exp(logAlphaBeta3All[Z3,2])
  
  
  logAlphaBeta3.XY <- logAlphaBeta3[1]*logAlphaBeta3[2]
  
  # Normal prior for interaction parameters (eta12, eta13, eta23 and eta123) 
  etaPriorPrec12 <- pow(etaPrior12[2], -2)
  eta12  ~ dnorm(etaPrior12[1], etaPriorPrec12)
  OddsFactor12 <- exp(eta12)
  
  etaPriorPrec13 <- pow(etaPrior13[2], -2)
  eta13  ~ dnorm(etaPrior13[1], etaPriorPrec13)
  OddsFactor13 <- exp(eta13)
  
  etaPriorPrec23 <- pow(etaPrior23[2], -2)
  eta23  ~ dnorm(etaPrior23[1], etaPriorPrec23)
  OddsFactor23 <- exp(eta23)
  
  etaPriorPrec123 <- pow(etaPrior123[2], -2)
  eta123  ~ dnorm(etaPrior123[1], etaPriorPrec123)
  OddsFactor123 <- exp(eta123)
  
  
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
  
  for(i in 1:Ndoses3){      
    #  Marginal odds for Agent 3 
    lin3.1[i] <- logAlphaBeta3[1] + exp(logAlphaBeta3[2])*log(Doses3[i,1]/DoseRef3)
    odds3[i] <- exp(lin3.1[i])
  } 
  
  #  combination model and probability of each category for each combination
  for(i in 1:Ndoses1){
    for(j in 1:Ndoses2){
      for(k in 1:Ndoses3){
      #  Odds without and with interaction 
      odds12.0[i,j,k] <- odds1[i] + odds2[j] + odds3[k] + odds1[i]*odds2[j]+ odds1[i]*odds3[k]+ odds2[j]*odds3[k] + odds1[i]*odds2[j]*odds3[k]
      odds12[i,j,k] <- odds12.0[i,j,k] * 
        exp(eta12*(Doses1[i,1]/DoseRef1)*(Doses2[j,1]/DoseRef2) + eta13*(Doses1[i,1]/DoseRef1)*(Doses3[k,1]/DoseRef3) + eta23*(Doses2[j,1]/DoseRef2)*(Doses3[k,1]/DoseRef3) + eta123*(Doses1[i,1]/DoseRef1)*(Doses2[j,1]/DoseRef2)*(Doses3[k,1]/DoseRef3))
      P12[i,j,k] <- odds12[i,j,k]/(1 + odds12[i,j,k])
      
      #  Prediction likelihood
      Ntox.pred[i,j,k] ~ dbin(P12[i,j,k],New.cohortsize)
      
      
      # for each dose, indicators for different toxicity categories
      # (means of these indicators correspond to category probabilities)   
      pCat[i,j,k,1] <- step(Pcutoffs1[1] - P12[i,j,k])
      for (l in 2:Ncat) {
        pCat[i,j,k,l] <- step(Pcutoffs1[l] - P12[i,j,k]) - sum(pCat[i,j,k,1:(l-1)])
      }
      
      #  Predictive probability of 0:New.cohortsize DLTs among New.cohortsize patients
      
      for (ll in 1:(New.cohortsize+1)){
        pred.Tox[i,j,k,ll] <- equals(ll-1,Ntox.pred[i,j,k])
      }
      
     }
    }
  }
  
  #binomial likelihoods (data)
  for (i in 1:Ncohorts){
    lin1.2[i] <- logAlphaBeta1[1] + exp(logAlphaBeta1[2])*log(DosesAdm1[i,1]/DoseRef1)
    lin1[i] <- lin1.2[i] + step(-20 - lin1.2[i]) * (-20 - lin1.2[i]) +  
      step(lin1.2[i] - 20) * (20 - lin1.2[i])
    odds11[i] <- exp(lin1[i]) 
    
    lin2.2[i] <- logAlphaBeta2[1] + exp(logAlphaBeta2[2])*log(DosesAdm2[i,1]/DoseRef2)
    lin2[i] <- lin2.2[i] + step(-20 - lin2.2[i]) * (-20 - lin2.2[i]) +  
      step(lin2.2[i] - 20) * (20 - lin2.2[i]) 
    odds22[i] <- exp(lin2[i])
    
    lin3.2[i] <- logAlphaBeta3[1] + exp(logAlphaBeta3[2])*log(DosesAdm3[i,1]/DoseRef3)
    lin3[i] <- lin3.2[i] + step(-20 - lin3.2[i]) * (-20 - lin3.2[i]) +  
      step(lin3.2[i] - 20) * (20 - lin3.2[i]) 
    odds33[i] <- exp(lin3[i])
    
    odds121.0[i] <- odds11[i] + odds22[i] + odds33[i] +odds11[i]*odds22[i] + odds11[i]*odds33[i] + odds22[i]*odds33[i] + odds11[i]*odds22[i]*odds33[i]         
    odds121[i] <-  odds121.0[i] *  
      exp(eta12*(DosesAdm1[i,1]/DoseRef1)*(DosesAdm2[i,1]/DoseRef2) + eta13*(DosesAdm1[i,1]/DoseRef1)*(DosesAdm3[i,1]/DoseRef3) + eta23*(DosesAdm2[i,1]/DoseRef2)*(DosesAdm3[i,1]/DoseRef3) + eta123*(DosesAdm1[i,1]/DoseRef1)*(DosesAdm2[i,1]/DoseRef2)*(DosesAdm3[i,1]/DoseRef3))
    P12.1[i] <- odds121[i]/(1 + odds121[i])
    
    Ntox[i,1] ~ dbin(P12.1[i], Npat[i,1]) 
  }  
  
  
}
  