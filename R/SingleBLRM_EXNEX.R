#' @title SingleBLRM.EXNEX Function 
#'
#' @description This function analyzes dose limiting toxicity (DLT) data using a Bayesian logistic regression model (BLRM) combined with the Exchangeability-nonexchangeability (EXNEX) approach for single agent Phase I Oncology trials. The EXNEX approach allows information to be borrowed across concurrent trials.
#'
#' @param trialData dose-DLT data: a list with length equal to the number of trials. Data for each trial is entered as a list object with 4 elements; DoseAdm, Npat, Ntox, trialName. DoseAdm is a vector of the doses administered. Npat are the number of participants treated at a given dose. Ntox are the number of DLTs observed for a given dose. trialName is the unique name used to identfy a specific trial.
#' @param Doses a vector of Doses for which inferential results will be provided, i.e., the provisional dose levels.
#' @param DoseRef refernce dose for the BLRM.
#' @param Prior.NEX is the nonexchangeability prior for the BLRM parameters, \eqn{log(\alpha)} and \eqn{log(\beta)}. This is a vector of length 5 for the NEX bivariate normal prior for, (\eqn{m_{log(\alpha)}}, \eqn{m_{log(\beta)}}, \eqn{s_{log(\alpha)}}, \eqn{s_{log(\beta)}}, \eqn{Corr(log(\alpha),log(\beta))}).
#' @param Prior.EX.mu1 is the EX prior specification for the \eqn{log(\alpha)} mean parameter. It is a vector of length 2, (\eqn{m_{\mu1}, s_{\mu1}}).
#' @param Prior.EX.mu2 is the EX prior specification for the \eqn{log(\beta)} mean parameter. It is a vector of length 2, (\eqn{m_{\mu2}, s_{\mu2}}).
#' @param Prior.EX.tau1 is the EX prior specification for the \eqn{log(\alpha)} standard deviation parameter. It is a vector of length 2, (\eqn{m_{\tau1}, s_{\tau1}}).
#' @param Prior.EX.tau2 is the EX prior specification for the \eqn{log(\beta)} standard deviation parameter. It is a vector of length 2, (\eqn{m_{\tau2}, s_{\tau2}}).
#' @param Prior.rho is the uniform prior specification for the correlation between \eqn{\tau1} and \eqn{\tau2}. It is a vector of length 2, and the default is set to (-1,1).
#' @param pj.EX is the proportion of exchangeability, resulting in a 1- pj.EX nonexchangeability proportion. Currenlty set to be the same for all strata. For example, pj.EX = 1 is full exchangeability, while pj.EX=0 is a stratified analysis.
#' @param Pint targeted interval probability cutoff, which will be used to define three intervals for the true DLT probability. The default is
#' (0.16, 0.33) defining three categories: underdosing ({[}0--0.16{]}), target ((0.16--0.33)), and 
#' overdosing ({[}0.33--1.00{]}).
#' @param MCMC four component vector of MCMC parameters n.burnin, n.iterations, n.chains, n.thin. See \link{R2WinBUGS} for details.
#' @param Rnd.seed Random seed for the R random stream. This is important for reproducibility. 
#' @param RSeed.WinBUGS Random seed for the WinBUGS random stream. This is important for reproducibility.
#' @param DIC logical: if TRUE, compute deviance, pD, and DIC. See \link{R2WinBUGS} for details. 
#' @param progress.bar type of progress bar. Possible values are "text", "gui", and "none". Type "text"
#' is displayed on the R console. See \link{rjags} for details.
#'
#' @details Statistical model:
#'
#'EXNEX for binary data:
#'The EXNEX apporach allows the use of co-data where stratum specific (e.g. study specific) parameters can be either exchangeable with parameters from other strata or nonexchangeable. For binary data this is the log-odds parameter, \eqn{\theta_j = log(\frac{\pi_j}{1-\pi_j})}.
#'
#'For each strata, j, we have weights \eqn{p_j} and \eqn{1-p_j}. Where the parameters are exchangeable (EX) with probability \eqn{p_j} and \eqn{\theta_j} follows a normal distribution with exchangeability parameters \eqn{\mu} and \eqn{\tau},
#'
#'\eqn{\theta_j \vert \mu, \tau \sim N(\mu,\tau^2)}
#'
#'and nonexchangeable (NEX) with probability \eqn{1-p_j} and strata-specific priors are used where,
#'
#'\eqn{\theta_j \sim N(m_j,v_j)}.
#'
#'BLRM and EXNEX:
#'
#'The dose-toxicity relationship will be modeled using a 2-parameter BLRM where for dose \eqn{d^*}, the number of participants with a dose limiting toxicity (DLT), \eqn{r_d}, in cohort of size \eqn{n_d} is binomial,
#'
#'\eqn{r_d \vert \pi_d \sim Binomial(\pi_d, n_d)}
#'
#'and for each stratum, j, the dose-toxicity model is logistic,
#'
#'\eqn{logit(\pi_{jd}) = log(\alpha_j) + \beta_j log(\frac{d}{d^*})}
#'
#'where \eqn{d} is the current dose under evaluation and \eqn{d^*} is the reference dose, which is used to scale the doses. The \eqn{\alpha, \beta > 0} are the parameters of the model such that \eqn{\alpha} is the therapy's odds of a DLT at \eqn{d^*}, and \eqn{\beta} is the increase in the log-odds of a DLT for a unit increase in the log-dose.
#'
#'The model is then confined to one EX distribution and one NEX distribution, with *a priori* \eqn{p_j} and \eqn{1-p_j}. Under \eqn{p_j = 1} parameters are fully exchangeable, while \eqn{p_j = 0} would result in a fully stratified analysis with no exchangeability.
#'
#'With probability \eqn{p_j} for EX, the logistic parameter vector
#'
#'\eqn{\theta_j = [log(\alpha_j),log(\beta_j)]}
#'
#'is exchangeable with the other strata parameters and follows a bivaraite normal distribution with mean vector \eqn{\mu} and covariance matrix \eqn{\Sigma};
#'
#'\eqn{\theta_j \vert \mu, \Sigma \sim BVN(\mu,\Sigma)}
#'
#'where
#'
#'\eqn{\mu = (\mu_1, \mu_2)} 
#'
#'and 
#'
#'\eqn{\Sigma = \left({\begin{array}{cc} \tau^2_1 & \rho\tau_1\tau_2 \\ \rho\tau_1\tau_2 &  \tau^2_2\end{array}}\right)}.
#'
#'With probability \eqn{1-p_j}, \eqn{\theta_j} is nonexchangeable with other strata and the prior distribution follows a bivariate normal distribution with mean vector \eqn{m_W} and covariance matrix \eqn{S_W},
#'
#'\eqn{\theta_j \sim BVN(m_W, S_W)}.
#'
#'Prior Specification:
#'The Bayesian approach requires the specification of the prior distributions for the model parameters \eqn{log(\alpha)} and \eqn{log(\beta)}. Below are the prior detials for both EX and NEX.
#'
#'Details on prior derivation can be found in the references below.
#'
#'@references 
#'Neuenschwander B, Wandel S, Roychoudhury S, Bailey S. Robust exchangeability designs for early phase clinical trials with multiple strata. Pharm Stat. 2016;15(2):123-34.
#' 
#'Neuenschwander B, Roychoudhury S, Schmidli H. On the use of co-data in clinical trials. Biopharmaceutical Research. 2016;8(3):345-54
#'
#'@return 
#'The function returns as a list the following main items:
#' \item{P}{Posterior summary of DLT probabilities}
#' \item{Pcat}{Posterior summary of DLT interval probabilities}
#' \item{PcatP}{P and Pcat combined}
#' \item{logAlphaBeta}{Posterior summary of (\eqn{\log(\alpha)}, \eqn{\log(\beta)})}
#' \item{PmixEXNEX}{Posterior mixture weigths}
#' \item{intervalPlot}{Interval plot for DLT rates}
#' \item{postProbPlot}{Estimates and 95\% interval plots for DLT rates}
#' \item{samplingPlot}{Summary of the sampling results including the interval DLT rates}
#' 
#' @examples 
#' \dontrun{
#' ## Example 1: Analysis of three concurrent phase I single agent trials 
#' 
#' threeSingleAgent<-SingleBLRM.EXNEX(
#'  DoseRef=24,
#'  Prior.NEX=c(log(0.13),0,2,1,0),
#'  Prior.EX.mu1=c(log(0.13), 1.982),
#'  Prior.EX.mu2=c(0, 0.992),
#'  Prior.EX.tau1=c(log(0.25), 0.3542),
#'  Prior.EX.tau2=c(log(0.125), 0.3542),
#'  Prior.rho=c(-1,1),
#'  MCMC=c(2500, 7500, 4, 1),
#'  Rnd.seed=1452671,
#'  RSeed.WinBUGS=1367,
#'  DIC=TRUE,
#'  progress.bar="text",
#'  pj.EX=0.5,
#'  Doses=c(4, 6, 12, 18, 24),
#'  trialData=list(
#'
#'    list(DosesAdm=c(6,12,24),
#'         Npat=c(3,3,3),
#'        Ntox=c(0,0,0),
#'         trialName="Trial A"),
#'
#'    list(DosesAdm=c(6,12,24),
#'         Npat=c(3,3,3),
#'         Ntox=c(0,1,0),
#'         trialName="Trial B"),
#'
#'    list(DosesAdm=c(6,10,11,15),
#'         Npat=c(3,3,2,3),
#'         Ntox=c(0,1,0,0),
#'         trialName="Trial C")
#'
#'  )
#')                          
#'}
#'
#'
#' @export


#######################
## BLRM+EXNEX function 
#######################

SingleBLRM.EXNEX<-function(
  
  ################
  ## Prior Values
  ################
  DoseRef=NULL,
  Prior.NEX=NULL,
  Prior.EX.mu1=NULL,
  Prior.EX.mu2=NULL,
  Prior.EX.tau1=NULL,
  Prior.EX.tau2=NULL,
  Prior.rho=c(-1,1),
  Pint=c(0.16,0.33),
  MCMC=c(2500, 7500, 4, 1),
  Rnd.seed=1452671,
  RSeed.WinBUGS=1367,
  DIC=TRUE,
  progress.bar="text",
  
  pj.EX=NULL, ## Proportion of exchangeability, the amount of nonexchangeability is 1- pj.EX, pj.EX = 1 is full exchangeability, while pj.EX=0 is a stratified analysis
  
  ## provisional dose levels to be evaluated
  Doses=NULL,
  
  ## Enter in data for the different trials. Trial Data is a list with a length that is equal to the number of different trials. Each trial is a list with three items, DosesAdm, Npat, and Ntox
  trialData=NULL
){
  
  ## Required libraries
  require(R2jags) ## interface for R and JAGS
  require(ggplot2) ## Used to create figures
  require(tidyr) ## Data manipulation


strataDosesAdm<-lapply(1:length(trialData),function(x){
  
  DosesAdm<-trialData[[x]]$DosesAdm
  
  return(DosesAdm)
  
})

strataNpat<-lapply(1:length(trialData),function(x){
  
  Npat<-trialData[[x]]$Npat
  
  return(Npat)
  
})



strataTox<-lapply(1:length(trialData),function(x){
  
  Ntox<-trialData[[x]]$Ntox
  
  return(Ntox)
  
})

strataLabels<-lapply(1:length(trialData),function(x){
  
  label<-trialData[[x]]$trialName
  
  return(label)
  
})

strataLabels<-do.call("cbind",strataLabels)

###################################################
## Function to combine data from multiple studies 
## into a format suitable for JAGS
###################################################

CombineStudies<-function(strataInfo){
  
  numStrata<-length(strataInfo)
  
  numCohots<-sapply(1:numStrata,function(x){
    numCohots<-length(strataInfo[[x]])
    numCohots
    
  })
  
  
  maxNumCohots<-max(numCohots)
  
  combined<-lapply(1:numStrata,function(x){
    
    if(length(strataInfo[[x]])< maxNumCohots){
      
      combined<-c(strataInfo[[x]],rep(NA,maxNumCohots-length(strataInfo[[x]])))
      
    }else{
      
      combined<-strataInfo[[x]]
      
    }
    
  })
  
  combined<-do.call("rbind",combined)
  
}
###################################################
###################################################


################################################
## Setup the data for the single BLRM with EXNEX
################################################


## matrix of the doses administered
DosesAdm<-CombineStudies(strataDosesAdm)


## matrix of toxicities observed at a given dose/study
Ntox<-CombineStudies(strataTox)



## matrix of participants treated at a given dose/study
Npat<-CombineStudies(strataNpat)



## determine the number of provisional doses
Ndoses<-length(Doses)

## determine the number of cohorts treated for each study
Ncohorts<-vapply(1:nrow(DosesAdm),function(x){
  
  tmp<-sum(!is.na(DosesAdm[x,]))
  
},FUN.VALUE=1)

## determine the number of strata, i.e. number of studies
Nstrata=nrow(DosesAdm)

## determine the number of categories for the mixture = Nexch+1
Nmix=2

## check that pj.EX is valid

if(pj.EX > 1 | pj.EX < 0){
  
  stop(print("Invalid entry for pj.EX. Value must be between 0 and 1."))
  
}



data<-list( Ncohorts=Ncohorts, # number of cohorts treated for study 1, 2, ...
            DoseRef=DoseRef, ## reference dose
            DosesAdm=DosesAdm,
            Npat=Npat,
            Ntox=Ntox,
            Ndoses=Ndoses, ## number of provisional dose levels for a given study
            Doses=Doses,
            Pint=Pint,
            pMix=c(pj.EX,1-pj.EX),
            Nmix=Nmix,
            Nstrata=Nstrata,
            Prior.NEX=Prior.NEX,
            Prior.EX.mu1=Prior.EX.mu1,
            Prior.EX.mu2=Prior.EX.mu2,
            Prior.EX.tau1=Prior.EX.tau1,
            Prior.EX.tau2=Prior.EX.tau2,
            Prior.rho=Prior.rho
)




##################################
# run BLRM+EXNEX function in JAGS
##################################

## select what parameters to output from JAGS
pars <- c("Pr.Cat","Pr.Tox","ProvDoseLelCat","log.alpha","log.beta","exch")



#Initial values (using same function that is in the OncoPh1BLRM pkg)
inits.fun = function(i){
  
  tmp.NEX <- matrix(NA, nrow=Nstrata, ncol=2)
  reTmp <- matrix(NA, nrow=Nstrata, ncol=2)
  
  
  for(i1 in 1:Nstrata){
    tmp.NEX[i1,1:2] <- mvtnorm::rmvnorm(1, sigma = diag(2))
    reTmp[i1,1:2] <- mvtnorm::rmvnorm(1, sigma = diag(2))
    
  }
  
  
  tau.norm<-c(rnorm(1),rnorm(1))
  
  exch.index<-rep(1,Nstrata)
  
  
  return(list(tmp.NEX = tmp.NEX,reTmp=reTmp,exch.index=exch.index))
  
}


# MCMC parameters
n.burnin = MCMC[1]
n.iter = MCMC[2]
n.chains = MCMC[3]
n.thin = MCMC[4]


set.seed(Rnd.seed)

## initial parameters
inits <- lapply(rep(1,n.chains),inits.fun)


fit<-R2jags::jags(
  data=data,
  inits= inits,
  parameters.to.save = pars,
  model=SingleBLRMEXNEX,
  n.chains=n.chains,n.burnin=n.burnin,n.iter=n.iter,n.thin=n.thin,
  jags.seed=RndSeed.WinBUGS,
  DIC= DIC,
  progress.bar=progress.bar
)


###################################################
## summary for the interval probabilities
###################################################
## loop over each strata 
out<-lapply(1:data$Nstrata,function(a){
  
  
## loop over the doses for a given strata
DLInfo<-lapply(1:length(data$Doses),function(b){
  
  doseLevl<-data$Doses[b]

## loop over each probability interval
PIInfo<-lapply(1:3,function(c){
  
 probInterval<- ifelse(c==1,"UD",
         ifelse(c==2, "TT",
                ifelse(c==3, "OD",NA)))
  
  
  ##Subset to a given probability interval for a given strata and dose level
  tmp<-data.frame(round(colMeans(data.frame(fit$BUGSoutput$sims.list$Pr.Cat[,a,b,c])),digits=4))
  
    
    name<-strataLabels[a]
    
    tmp$Trial<-name
    
    tmp$probInterval<-probInterval
    
    tmp$doseLevel<-doseLevl
    
    return(tmp)
    
})

tmp1<-do.call("rbind",PIInfo)

return(tmp1)
})

tmp2<-do.call("rbind",DLInfo)

return(tmp2)


})


Pcat<-do.call("rbind",out)

row.names(Pcat)<-seq(1,nrow(Pcat),1)

colnames(Pcat)<-c("prob","Trial","probInterval","Dose")

## create a copy for the plot
out<-Pcat

##################################################
##################################################

###################################################
## Plot the interval probabilties
###################################################


## add a reference line for EWOC and MTD (i.e., TT > 0.5) probabilities
refLine<-data.frame("probInterval"=c("TT","OD","UD"),"ref"=c(0.5,0.25,NA))


out$fill <- ifelse(out$prob > 0.25 & out$probInterval=="OD","firebrick1",
                              ifelse(out$prob > 0.5 & out$probInterval=="TT","chartreuse2","chartreuse4"))

## plot results, note you will get a warning due to the way fill is used to create different groups and an NA exists for the fill color if the condition doesn't exist
intervalPlot<-ggplot(data=out,aes(x=factor(Dose),y=prob, fill=fill))+
  geom_bar(stat="identity")+ ## create a bar plot
  geom_hline(data=refLine,aes(yintercept=ref), linetype=2, color=c("black"), size=1)+ ## add reference lines
  labs(y="Interval Probabilities",x=paste0("Dose"))+ ## create lables
  ggtitle("")+ ## title plot
  facet_grid(probInterval~Trial,scales="free_x")+ ## facet the results based on the probability interval
  scale_fill_identity()+ ## will the results with the conditional fill
  scale_y_continuous(breaks =c(0,0.25,0.5,0.75,1),limits=c(0,1))+ ## create axis breaks
  theme_bw()+ ## theme black and white
  theme(axis.text=element_text(size=14), ## axis font size
        axis.title=element_text(size=16,face="bold"), ## axis title font sizes
        strip.text = element_text(size = 14)) ## font size for facet







###################################################
## table of interval probabilities
###################################################
## create a copy of the long format
PcatLong<-Pcat
## covert from long to wide
Pcat<- Pcat %>% 
  pivot_wider(names_from=probInterval, values_from=prob) %>%
  as.data.frame()




###################################################
## summary for the probabilities
###################################################
## loop over each strata 
outProb<-lapply(1:data$Nstrata,function(a){
 
    
    doseLevl<-data$Doses
    
    tmp<-data.frame("mean"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$Pr.Tox[,a,])),digits=4))
    
    tmp$SD<-apply(fit$BUGSoutput$sims.list$Pr.Tox[,a,],2,sd)
    
    tmp$median<-apply(fit$BUGSoutput$sims.list$Pr.Tox[,a,],2,median)
    
    tmp$Q2.5<-apply(fit$BUGSoutput$sims.list$Pr.Tox[,a,],2,function(x){quantile(x,probs=0.025)})
    
    tmp$Q97.5<-apply(fit$BUGSoutput$sims.list$Pr.Tox[,a,],2,function(x){quantile(x,probs=0.975)})
      
    tmp<-round(tmp,digits=4)
    

    ## subset to the doses of interest
    
    tmp<-tmp[c(1:length(doseLevl)),]
    
    tmp<-cbind.data.frame("Dose"=doseLevl, tmp)
      

    
    ## add in strata names 
    
    name<-strataLabels[a]
    
    tmp<-cbind.data.frame("Trial"=name, tmp)
      
      return(tmp)
      

  
})


P<-do.call("rbind",outProb)



PcatP<-merge(Pcat,P,by=c("Trial","Dose"))


PcatP<-PcatP[order(PcatP$Trial,PcatP$Dose),]



###############################################
## Plot of the DLT probabilties
###############################################

## plot the results
postProbPlot<-ggplot(data=P,aes(x=factor(Dose),y=mean,color=Trial))+
  geom_point(size=2.5,position = position_dodge(0.3))+ ## create a points for the mean DLT
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0.2,position = position_dodge(0.3),size=1)+ ## create bars for the 95% intervals
  labs(y="DLT Rate",x=paste0("Dose"),color="")+ ## labels
  ggtitle("Posterior Probability - Mean DLT Rate with 2.5 and 97.5 Percentiles")+ ## title
  scale_color_brewer(palette="Set1")+ ## color points
  scale_y_continuous(breaks =c(0,0.25,0.5,0.75,1),limits=c(0,1))+ ## set y axis breaks
  theme_bw()+ ## black and white theme
  theme(legend.position="bottom", ## location of legend
        axis.text=element_text(size=14), ## axis font size
        axis.title=element_text(size=16,face="bold"), ## axis title font size 
        legend.text=element_text(size=14)) ## legend font size



###################################################
## posterior summary for log alpha
###################################################
## loop over each strata 
outLogAlpha<-lapply(1:data$Nstrata,function(a){
  
  tmp<-data.frame("mean"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$log.alpha[,a])),digits=4))
  
  tmp$SD<-sd(fit$BUGSoutput$sims.list$log.alpha[,a])
  
  tmp$median<-median(fit$BUGSoutput$sims.list$log.alpha[,a])
  
  tmp$Q2.5<-quantile(fit$BUGSoutput$sims.list$log.alpha[,a],probs=0.025)
  
  tmp$Q97.5<-quantile(fit$BUGSoutput$sims.list$log.alpha[,a],probs=0.975)
  
  tmp<-round(tmp,digits=4)
  
  name<-strataLabels[a]
  
  tmp<-cbind.data.frame("Trial"=name, tmp)
  
  row.names(tmp)<-"logAlpha"
  
  
  return(tmp)
  
  
  
})


outLogAlpha<-do.call("rbind",outLogAlpha)



###################################################
## posterior summary for log alpha
###################################################
## loop over each strata 
outLogBeta<-lapply(1:data$Nstrata,function(a){
  
  tmp<-data.frame("mean"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$log.beta[,a])),digits=4))
  
  tmp$SD<-sd(fit$BUGSoutput$sims.list$log.beta[,a])
  
  tmp$median<-median(fit$BUGSoutput$sims.list$log.beta[,a])
  
  tmp$Q2.5<-quantile(fit$BUGSoutput$sims.list$log.beta[,a],probs=0.025)
  
  tmp$Q97.5<-quantile(fit$BUGSoutput$sims.list$log.beta[,a],probs=0.975)
  
  tmp<-round(tmp,digits=4)
  
  name<-strataLabels[a]
  
  tmp<-cbind.data.frame("Trial"=name, tmp)
  
  row.names(tmp)<-"logBeta"
  
  
  return(tmp)
  
  
  
})


outLogBeta<-do.call("rbind",outLogBeta)

logAlphaBeta<-rbind.data.frame(outLogAlpha,outLogBeta)




###################################################
## summary for the posterior mixture probabilities
## 1 = EX, 2 = NEX
###################################################
## loop over each strata 


out<-lapply(1:data$Nstrata,function(a){
  
  

      
      ##Determine how often EX or NEX was used
      tmp<-data.frame("PmixEX"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$exch[,a,1])),digits=4))
      tmp$PmixNEX<-round(colMeans(data.frame(fit$BUGSoutput$sims.list$exch[,a,2])),digits=4)
      
      
      name<-strataLabels[a]
      
      tmp$Trial<-name
      
     
      
      return(tmp)
      
   
  
})


PmixEXNEX<-do.call("rbind",out)

row.names(PmixEXNEX)<-seq(1,nrow(PmixEXNEX),1)




##################################################
##################################################

#################################################
## Summarize the interval proibabilites
#################################################

doseLevels<-data$Doses


out<-lapply(1:Nstrata,function(z){
  
  tmp<-lapply(1:Ndoses,function(y){
    
    udResults<-which(fit$BUGSoutput$sims.list$Pr.Cat[,z,y,1]==1)
    udID<-rep("UD",length(udResults))
    
    ttResults<-which(fit$BUGSoutput$sims.list$Pr.Cat[,z,y,2]==1)
    ttID<-rep("TT",length(ttResults))
    
    odResults<-which(fit$BUGSoutput$sims.list$Pr.Cat[,z,y,3]==1)
    odID<-rep("OD",length(odResults))
    
    
    tmp<-data.frame("Dose"=doseLevels[y],"Interval"=c(udID,ttID,odID),"rowLoc"=c(udResults,ttResults,odResults))
    
    return(tmp)
  })
  
  tmp<-do.call("rbind",tmp)
  
  tmp<-data.frame("Strata"=strataLabels[z],tmp)
  
  return(tmp)
  
})


out<-do.call("rbind",out)




out1<-lapply(1:Nstrata,function(z){
  
  tmp<-lapply(1:Ndoses,function(y){
    
    intervalInfo<-out[out$Strata==strataLabels[z] & out$Dose==doseLevels[y],c("Interval", "rowLoc")]
    
    intervalInfo<-intervalInfo[order(intervalInfo$rowLoc,decreasing=FALSE),]
    
    tmp<-data.frame("Dose"=doseLevels[y],"Interval"=intervalInfo$Interval,"prob"=fit$BUGSoutput$sims.list$Pr.Tox[,z,y])
    
  })
  tmp<-do.call("rbind",tmp)
  
  tmp<-data.frame("Strata"=strataLabels[z],tmp)
  
  return(tmp)
  
})

out1<-do.call("rbind",out1)

numSims<-nrow(out1)/Nstrata/Ndoses




######################################################
## Determine the text locations for each of the facets
######################################################

textLoc<-lapply(1:Nstrata,function(z){
  
  tmp<-lapply(1:Ndoses,function(y){
    

    props<-table(round(out1$prob[out1$Strata==strataLabels[z] & out1$Dose==doseLevels[y]],digits=2))/numSims
    
    propMax<-max(props)
    
    
    tmp<-data.frame("Dose"=doseLevels[y],"propMax"=propMax)
    
  })
  tmp<-do.call("rbind",tmp)
  
  tmp<-data.frame("Strata"=strataLabels[z],tmp)
  
  return(tmp)
  
  
})

textLoc<-do.call("rbind",textLoc)


## Add in the facet label information
textInfo<-PcatLong

textInfo$label<-paste("P(",textInfo$probInterval,") = ",round(textInfo$prob,2),sep="")

textInfo<-lapply(1:Nstrata,function(z){
  
  tmp<-lapply(1:Ndoses,function(y){
  
  tmpLab<-textInfo$label[textInfo$Trial==strataLabels[z] & textInfo$Dose==doseLevels[y]]
  
  tmpLab<-paste(tmpLab[1],tmpLab[2],tmpLab[3],sep = "\n")
  
  tmp<-data.frame("Dose"=doseLevels[y],"tmpLab"=tmpLab)
  
})
  
  tmp<-do.call("rbind",tmp)
  
  tmp<-data.frame("Strata"=strataLabels[z],tmp)
  
})

textInfo<-do.call("rbind",textInfo)


dat_text<-merge(textLoc,textInfo,by=c("Strata","Dose"))
dat_text$Interval<-"UD"

# ## text for UD
# dat_text_ud <- Pcat[Pcat$probInterval=="UD",]
# dat_text_ud$prob<-round(dat_text_ud$prob,digits=2)
# colnames(dat_text_ud)<-c("prob", "Strata", "Interval", "Dose")
# dat_text_ud<-merge(dat_text_ud,textLoc,by=c("Strata","Dose"))
# 
# ## text for TT
# dat_text_tt <- Pcat[Pcat$probInterval=="TT",]
# dat_text_tt$prob<-round(dat_text_tt$prob,digits=2)
# colnames(dat_text_tt)<-c("prob", "Strata", "Interval", "Dose")
# dat_text_tt<-merge(dat_text_tt,textLoc,by=c("Strata","Dose"))
# 
# ## text for OD
# dat_text_od <- Pcat[Pcat$probInterval=="OD",]
# dat_text_od$prob<-round(dat_text_od$prob,digits=2)
# colnames(dat_text_od)<-c("prob", "Strata", "Interval", "Dose")
# dat_text_od<-merge(dat_text_od,textLoc,by=c("Strata","Dose"))


samplingPlot<-ggplot(out1,aes(x=prob,fill=Interval))+
  geom_histogram(aes(y = ..count../numSims),binwidth=0.01,alpha=0.7)+
  scale_fill_manual(values=c("firebrick1","chartreuse2","chartreuse4"))+
  labs(x="Probability of Toxicity",y="Percent of Samples")+
  facet_wrap(Strata~Dose,scales="free_y",nrow=Nstrata)+
  geom_vline(xintercept=Pint[1],lty=2)+
  geom_vline(xintercept=Pint[2],lty=2)+
  geom_text(data    = dat_text,mapping = aes(x = Pint[2]+((1-Pint[2])/2), y = propMax-(propMax*0.5), label = gsub("\\\\n", "\n", tmpLab)),size=3.5)+ ## add probability interval info and rescale for each facet
  scale_x_continuous(breaks=seq(0,1,0.2))+
  ggtitle("Posterior Sampling Results")+
  guides(fill=guide_legend(title="Interval"))+
  theme_bw()+
  theme(legend.position = "bottom")
  
 
  
return(list(P=P,Pcat=Pcat,PcatP=PcatP,logAlphaBeta=logAlphaBeta,PmixEXNEX=PmixEXNEX,intervalPlot=intervalPlot,postProbPlot=postProbPlot,samplingPlot=samplingPlot))

}











