#' @title Combo2BLRM.EXNEX Function 
#'
#' @description This function analyzes dose limiting toxicity (DLT) data using a Bayesian logistic regression model (BLRM) combined with the Exchangeability-nonexchangeability (EXNEX) approach for 2 agent Phase I Oncology trials. The EXNEX approach allows information to be borrowed across concurrent trials.
#'
#' @param trialData dose-DLT data: a list with length equal to the number of trials. Data for each trial is entered as a list object with 5 elements; DoseAdm1,DoseAdm2, Npat, Ntox, trialName. DoseAdm1 is a vector of the doses administered for agent 1. DoseAdm2 is a vector of the doses administered for agent 2. Npat are the number of participants treated at a given dose. Ntox are the number of DLTs observed for a given dose. trialName is the unique name used to identfy a specific trial.
#' @param agent1 is a character string with the name/abbreviation for agent 1.
#' @param Doses1 a vector of doses for agent 1 for which inferential results will be provided, i.e., the provisional dose levels.
#' @param DoseRef1 refernce dose for agent 1 used in the BLRM.
#' @param Prior.NEX.Agent1 is the nonexchangeability prior for agent 1 used in the combination BLRM, \eqn{log(\alpha)} and \eqn{log(\beta)}. This is a vector of length 5 for the NEX bivariate normal prior for, (\eqn{m_{log(\alpha)}}, \eqn{m_{log(\beta)}}, \eqn{s_{log(\alpha)}}, \eqn{s_{log(\beta)}}, \eqn{Corr(log(\alpha),log(\beta))}).
#' @param Prior.EX.mu1.Agent1 is the EX prior specification for Agent 1 for the \eqn{log(\alpha)} mean parameter. It is a vector of length 2, (\eqn{m_{\mu1}, s_{\mu1}}).
#' @param Prior.EX.mu2.Agent1 is the EX prior specification for Agent 1 for the \eqn{log(\beta)} mean parameter. It is a vector of length 2, (\eqn{m_{\mu2}, s_{\mu2}}).
#' @param Prior.EX.tau1.Agent1 is the EX prior specification for Agent 1 for the \eqn{log(\alpha)} standard deviation parameter. It is a vector of length 2, (\eqn{m_{\tau1}, s_{\tau1}}).
#' @param Prior.EX.tau2.Agent1 is the EX prior specification for Agent 1 for the \eqn{log(\beta)} standard deviation parameter. It is a vector of length 2, (\eqn{m_{\tau2}, s_{\tau2}}).
#' @param Prior.rho.Agent1 is the uniform prior specification for Agent 1 for the correlation between \eqn{\tau1} and \eqn{\tau2}. It is a vector of length 2, and the default is set to (-1,1).
#' @param agent2 is a character string with the name/abbreviation for agent 2.
#' @param Doses2 a vector of doses for agent 2 for which inferential results will be provided, i.e., the provisional dose levels.
#' @param DoseRef2 refernce dose for agent 2 used in the BLRM.
#' @param Prior.NEX.Agent2 is the nonexchangeability prior for agent 2 used in the combination BLRM, \eqn{log(\alpha)} and \eqn{log(\beta)}. This is a vector of length 5 for the NEX bivariate normal prior for, (\eqn{m_{log(\alpha)}}, \eqn{m_{log(\beta)}}, \eqn{s_{log(\alpha)}}, \eqn{s_{log(\beta)}}, \eqn{Corr(log(\alpha),log(\beta))}).
#' @param Prior.EX.mu1.Agent2 is the EX prior specification for Agent 2 for the \eqn{log(\alpha)} mean parameter. It is a vector of length 2, (\eqn{m_{\mu1}, s_{\mu1}}).
#' @param Prior.EX.mu2.Agent2 is the EX prior specification for Agent 2 for the \eqn{log(\beta)} mean parameter. It is a vector of length 2, (\eqn{m_{\mu2}, s_{\mu2}}).
#' @param Prior.EX.tau1.Agent2 is the EX prior specification for Agent 2 for the \eqn{log(\alpha)} standard deviation parameter. It is a vector of length 2, (\eqn{m_{\tau1}, s_{\tau1}}).
#' @param Prior.EX.tau2.Agent2 is the EX prior specification for Agent 2 for the \eqn{log(\beta)} standard deviation parameter. It is a vector of length 2, (\eqn{m_{\tau2}, s_{\tau2}}).
#' @param Prior.rho.Agent2 is the uniform prior specification for Agent 2 for the correlation between \eqn{\tau1} and \eqn{\tau2}. It is a vector of length 2, and the default is set to (-1,1).
#' @param Prior.EX.eta.mu is the mean EX prior specification for the interaction parameter, i.e., log odds multiplier \eqn{\eta}. It is a vector of length 2, (\eqn{m_{\eta}, s_{\eta}}).
#' @param Prior.EX.eta.tau is the standard deviation EX prior specification for the interaction parameter, i.e., log odds multiplier \eqn{\eta}. It is a vector of length 2, (\eqn{m_{\tau}, s_{\tau}}).
#' @param Prior.NEX.eta is the nonexchangeability prior specification for the interaction parameter, i.e., log odds multiplier \eqn{\eta}. It is a vector of length 2, (\eqn{m_{\eta}, s_{\eta}}).
#' @param pj.EX is the proportion of exchangeability, resulting in a 1- pj.EX nonexchangeability proportion. For example, pj.EX = 1 is full exchangeability, while pj.EX=0 is a stratified analysis.
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
#'Details on the statistical model are provided below and additional details are provided in the references.
#'
#'Prior Specification:
#'The Bayesian approach requires the specification of the prior distributions for the model parameters for each agent and their interaction, e.g., \eqn{log(\alpha)}, \eqn{log(\beta)}, and \eqn{log(\eta)}. These priors need to be specified for both EX and NEX.
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
#'Combination BLRM with 2 agents and EXNEX:
#'
#'The dose-toxicity relationship will be modeled using a 5-parameter BLRM where for a combination dose of agent 1 \eqn{d_1^*} and agent 2 \eqn{d_2^*}, the number of participants with a dose limiting toxicity (DLT), \eqn{r_d}, in cohort of size \eqn{n_d} is binomial,
#'
#'\eqn{r_d \vert \pi_d \sim Binomial(\pi_d, n_d)}
#'
#'and for each stratum, j, and therapy the dose-toxicity model is logistic,
#'
#'\eqn{logit(\pi_{1jd}) = log(\alpha_{1j}) + \beta_{1j} log(\frac{d_1}{d_1^*})}
#'
#'\eqn{logit(\pi_{2jd}) = log(\alpha_{2j}) + \beta_{2j} log(\frac{d_2}{d_2^*})}
#'
#'where \eqn{d_1} and \eqn{d_2} are the current combination doses under evaluation and \eqn{d_1^*} and \eqn{d_2^*} are the reference doses, which are used to scale the doses. The \eqn{\alpha, \beta > 0} are the parameters of the model such that \eqn{\alpha} is the therapy's odds of a DLT at \eqn{d^*}, and \eqn{\beta} is the increase in the log-odds of a DLT for a unit increase in the log-dose.
#'
#'The dose-DLT relationship of the combination is defined as 
#'
#'\eqn{Odds(\pi_{12}(d_1,d_2)) = \frac{\pi_{12}(d_1,d_2)}{1-\pi_{12}(d_1,d_2)}}
#'\eqn{ = exp(\eta \frac{d_1}{d_1^*}\frac{d_2}{d_2^*})[\frac{\pi_1(d_1)+\pi_2(d_2)-\pi_1(d_1)\pi_2(d_2)}{(1-\pi_1(d_1))(1-\pi_2(d_2))}]}
#'
#'where \eqn{\eta} is the interaction paramter between agent 1 and agent 2. The Bayesian approach requires the specification of prior of prior distributions for the 5 model parameters \eqn{log(\alpha_{1})}, \eqn{log(\alpha_{2})}, \eqn{log(\beta_{1})}, \eqn{log(\beta_{2})}, and the interaction parameter \eqn{\eta}.
#'
#'The model is then confined to one EX distribution and one NEX distribution, with *a priori* \eqn{p_j} and \eqn{1-p_j}. Under \eqn{p_j = 1} parameters are fully exchangeable, while \eqn{p_j = 0} would result in a fully stratified analysis with no exchangeability. Note currenty only one weight can be specified and it is applied acorss all strata, i.e., \eqn{j=1}.
#'
#'With probability \eqn{p_j} for EX, the logistic parameter vector for agent 1
#'
#'\eqn{\theta_{1j} = [log(\alpha_{1j}),log(\beta_{1j})]}
#'
#'is exchangeable with the other strata parameters and follows a bivaraite normal distribution with mean vector \eqn{\mu_1} and covariance matrix \eqn{\Sigma_1};
#'
#'\eqn{\theta_{1j} \sim BVN(\mu_{1},\Sigma_{1})}
#'
#'where
#'
#'\eqn{\mu_{1} = (\mu_{11}, \mu_{12})} 
#'
#'and 
#'
#'\eqn{\Sigma_1 = \left({\begin{array}{cc} \tau^2_{11} & \rho\tau_{11}\tau_{12} \\ \rho\tau_{11}\tau_{12} &  \tau^2_{12}\end{array}}\right)}.
#'
#'the logistic parameter vector for agent 2,
#'
#'#'\eqn{\theta_{2j} = [log(\alpha_2j),log(\beta_2j)]}
#'
#'is exchangeable with the other strata parameters and follows a bivaraite normal distribution with mean vector \eqn{\mu_2} and covariance matrix \eqn{\Sigma_2};
#'
#'\eqn{\theta_{2j} \sim BVN(\mu_{2},\Sigma_2)}
#'
#'where
#'
#'\eqn{\mu_2 = (\mu_{21}, \mu_{22})} 
#'
#'and 
#'
#'\eqn{\Sigma_2 = \left({\begin{array}{cc} \tau^2_{21} & \rho\tau_{21}\tau_{22} \\ \rho\tau_{21}\tau_{22} &  \tau^2_{22}\end{array}}\right)}.
#'
#' The interaction paramter \eqn{\eta_{j}} is exchangeable with the other strata and follows a normal distribution,
#'
#'\eqn{\eta_{j} = \sim N(\mu_{\eta}, \tau^2_{\eta})}.
#'
#'With probability \eqn{1-p_j}, \eqn{\theta_{1j}}, \eqn{\theta_{2j}} and, \eqn{\eta_{j}} are nonexchangeable with other strata. The prior distributions for \eqn{\theta_{1j}} and \eqn{\theta_{2j}} follow a bivariate normal distribution with mean vector \eqn{m_{1W}} and covariance matrix \eqn{S_{1W}} for agent 1 and mean vector \eqn{m_{2W}} and covariance matrix \eqn{S_{2W}} for agent 2,
#'
#'\eqn{\theta_{1j} \sim BVN(m_{1W}, S_{1W})},
#'
#'\eqn{\theta_{2j} \sim BVN(m_{2W}, S_{2W})}.
#'
#'The prior distribution for \eqn{\eta_{j}} follows a normal distribution with
#'
#'\eqn{\eta_{j} \sim N(\mu_{\eta_{0}}, \tau^2_{\eta_{0}})}.
#'
#'
#'Details on deriving the priors can be found in the references below.
#'
#'@references 
#'Neuenschwander B, Wandel S, Roychoudhury S, Bailey S. Robust exchangeability designs for early phase clinical trials with multiple strata. Pharm Stat. 2016;15(2):123-34.
#' 
#'Neuenschwander B, Roychoudhury S, Schmidli H. On the use of co-data in clinical trials. Biopharmaceutical Research. 2016;8(3):345-54
#'
#'Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. A Bayesian Industry Approach to Phase I Combination Trials in Oncology. Statistical Methods for Drug Combination Studies. Taylor & Francis, 2015.
#'
#'@return 
#'The function returns as a list the following main items:
#' \item{P}{Posterior summary of DLT probabilities for the different dose combinations}
#' \item{Pcat}{Posterior summary of DLT interval probabilities for the different dose combinations}
#' \item{PcatP}{P and Pcat combined}
#' \item{logAlphaBeta}{Posterior summary of (\eqn{\log(\alpha)}, \eqn{\log(\beta)}) for agent 1 and 2}
#' \item{PmixEXNEX}{Posterior mixture weigths}
#' \item{intervalPlot}{Interval plots for DLT rates. This is output as a list equal to the length of the number of trials/strata.}
#' \item{postProbPlot}{Estimates and 95\% interval plots for DLT rates. This is output as a list equal to the length of the number of trials/strata.}
#' 
#' @examples 
#' \dontrun{
#' ## Example 1: Analysis of three concurrent combination phase I trials with 2 agents
#' 
#' threeTrial2Agent<-Combo2BLRM.EXNEX(
#'  DoseRef1=16,
#'  Prior.NEX.Agent1=c(log(0.15/0.85),0,2,1,0),
#'  Prior.EX.mu1.Agent1=c(log(0.15/0.85),2),
#'  Prior.EX.mu2.Agent1=c(0,1),
#'  Prior.EX.tau1.Agent1=c( log(0.25), log(2)/1.96),
#'  Prior.EX.tau2.Agent1=c( log(0.125), log(2)/1.96),
#'  Prior.rho.Agent1=c(-1,1),
#'  DoseRef2=2,
#'  Prior.NEX.Agent2=c(log(0.1/0.9),0,2,1,0),
#'  Prior.EX.mu1.Agent2=c(log(0.1/0.9),2),
#'  Prior.EX.mu2.Agent2=c(0,1),
#'  Prior.EX.tau1.Agent2=c( log(0.25), log(2)/1.96),
#'  Prior.EX.tau2.Agent2=c( log(0.125), log(2)/1.96),
#'  Prior.rho.Agent2=c(-1,1),
#'  Prior.EX.eta.mu=c(0, 1.121),
#'  Prior.NEX.eta=c(0, 1.121),
#'  Prior.EX.eta.tau=c( log(0.25), log(2)/1.96),
#'  Pint=c(0.16,0.33),
#'  MCMC=c(2500, 7500, 4, 1),
#'  Rnd.seed=1452671,
#'  RSeed.WinBUGS=1367,
#'  DIC=TRUE,
#'  progress.bar="text",
#'  agent1="PF1",
#'  agent2="PF2",
#'  pj.EX=0.5,
#'  Doses1=c(2,4,8,16,36),
#'  Doses2=c(0,1,2),
#'  trialData=list(
#'    
#'   list(DosesAdm1=c(1,2,5,10,25,50),
#'         DosesAdm2=c(2,2,2,2,2,2),
#'         Npat=c(3,3,2,3,6,3),
#'         Ntox=c(0,0,0,0,0,2),
#'         trialName="trialAB1"),
#'    
#'    list(DosesAdm1=c(5,10,10,10,20),
#'         DosesAdm2=c(2,2,2.5,3,2),
#'         Npat=c(3,3,3,3,3),
#'         Ntox=c(0,0,0,1,0),
#'         trialName="TrialAB2"),
#'    
#'    
#'    list(DosesAdm1=c(1,2,4,8),
#'         DosesAdm2=c(0,0,0,0),
#'         Npat=c(3,3,6,3),
#'         Ntox=c(0,0,1,0),
#'         trialName="trialA")
#'  )
#'  
#')                          
#'}
#'
#'
#' @export





#######################
## BLRM+EXNEX function 
#######################

Combo2BLRM.EXNEX<-function(
  
  ## labels
  agent1=NULL,
  agent2=NULL,
  ################
  ## Prior Values
  ################
  DoseRef1=NULL,
  Prior.NEX.Agent1=NULL,
  Prior.EX.mu1.Agent1=NULL,
  Prior.EX.mu2.Agent1=NULL,
  Prior.EX.tau1.Agent1=NULL,
  Prior.EX.tau2.Agent1=NULL,
  Prior.rho.Agent1=c(-1,1),
  DoseRef2=NULL,
  Prior.NEX.Agent2=NULL,
  Prior.EX.mu1.Agent2=NULL,
  Prior.EX.mu2.Agent2=NULL,
  Prior.EX.tau1.Agent2=NULL,
  Prior.EX.tau2.Agent2=NULL,
  Prior.rho.Agent2=c(-1,1),
  Prior.EX.eta.mu=NULL,
  Prior.EX.eta.tau=NULL, # heterogeneity of the interactions among the studies
  Prior.NEX.eta=NULL,
  Pint=c(0.16,0.33),
  MCMC=c(2500, 7500, 4, 1),
  Rnd.seed=1452671,
  RSeed.WinBUGS=1367,
  DIC=TRUE,
  progress.bar="text",
  
  pj.EX=NULL, ## Proportion of exchangeability, the amount of nonexchangeability is 1- pj.EX, pj.EX = 1 is full exchangeability, while pj.EX=0 is a stratified analysis
  
  ## provisional dose levels to be evaluated for each agent
  Doses1=NULL,
  Doses2=NULL,
  
  ## Enter in data for the different trials. Trial Data is a list with a length that is equal to the number of different trials. Each trial is a list with four items DosesAdm1, DosesAdm2, Npat, Ntox, and trialName. If only single agent data is available the doses for the other agent should be entered in as 0.
  trialData=NULL
){
  
  
  ## Required libraries
  require(R2jags) ## interface for R and JAGS
  require(ggplot2) ## Used to create figures
  require(tidyr) ## Data manipulation


strataDosesAdm1<-lapply(1:length(trialData),function(x){
  
  DosesAdm<-trialData[[x]]$DosesAdm1
  
  return(DosesAdm)
  
})

strataDosesAdm2<-lapply(1:length(trialData),function(x){
  
  DosesAdm<-trialData[[x]]$DosesAdm2
  
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


## matrix of the doses administered for Agent 1
DosesAdm1<-CombineStudies(strataDosesAdm1)

## matrix of the doses administered for Agent 2
DosesAdm2<-CombineStudies(strataDosesAdm2)


## matrix of toxicities observed at a given dose/study
Ntox<-CombineStudies(strataTox)



## matrix of participants treated at a given dose/study
Npat<-CombineStudies(strataNpat)



## determine the number of provisional doses for agent 1
Ndoses1<-length(Doses1)

## determine the number of provisional doses for agent 2
Ndoses2<-length(Doses2)

## determine the number of cohorts treated for each study
Ncohorts<-vapply(1:nrow(DosesAdm1),function(x){
  
  tmp<-sum(!is.na(DosesAdm1[x,]))
  
},FUN.VALUE=1)


## determine the number of strata, i.e. number of studies
Nstrata=nrow(DosesAdm1)

## determine the number of categories for the mixture = Nexch+1
Nmix=2


data<-list( Ncohorts=Ncohorts, # number of cohorts treated for study 1, 2, ...
            DoseRef1=DoseRef1, ## reference dose for agent 1
            DoseRef2=DoseRef2, ## reference dose for agent 2
            DosesAdm1=DosesAdm1, # doses administered for agent 1
            DosesAdm2=DosesAdm2, # doses administered for agent 2
            Npat=Npat, # num pts treated, matrix with separate rows for each strata
            Ntox=Ntox, # num pts with tox, matrix wiht separate rows for each strata
            Ndoses1=Ndoses1, ## number of provisional dose levels for agent 1
            Ndoses2=Ndoses2,## number of provisional dose levels for agent 2
            Doses1=Doses1, # provisional doses for agent 1
            Doses2=Doses2, # provisional doses for agent 2
            Pint=Pint, # probability intervals
            pMix=c(pj.EX,1-pj.EX), # mixture probability
            Nmix=Nmix, # number of distributions Nexch+1
            Nstrata=Nstrata, # number of strata
            Prior.NEX.Agent1=Prior.NEX.Agent1,
            Prior.EX.mu1.Agent1=Prior.EX.mu1.Agent1,
            Prior.EX.mu2.Agent1=Prior.EX.mu2.Agent1,
            Prior.EX.tau1.Agent1=Prior.EX.tau1.Agent1,
            Prior.EX.tau2.Agent1=Prior.EX.tau2.Agent1,
            Prior.rho.Agent1=Prior.rho.Agent1,
            Prior.NEX.Agent2=Prior.NEX.Agent2,
            Prior.EX.mu1.Agent2=Prior.EX.mu1.Agent2,
            Prior.EX.mu2.Agent2=Prior.EX.mu2.Agent2,
            Prior.EX.tau1.Agent2=Prior.EX.tau1.Agent2,
            Prior.EX.tau2.Agent2=Prior.EX.tau2.Agent2,
            Prior.rho.Agent2=Prior.rho.Agent2,
            Prior.EX.eta.mu=Prior.EX.eta.mu,
            Prior.NEX.eta=Prior.NEX.eta,
            Prior.EX.eta.tau=Prior.EX.eta.tau # heterogeneity of the interactions among the studies

)




##################################
# run BLRM+EXNEX function in JAGS
##################################

## select what parameters to output from JAGS
pars <- c("Pr.Cat","Pr.Tox","log.alpha.Agent1","log.alpha.Agent2","log.beta.Agent1","log.beta.Agent2","eta","OddsFactor","exch")



#Initial values (using same function that is in the OncoPh1BLRM pkg)
inits.fun = function(i){
  
  tmp.NEX.Agent1 <- matrix(NA, nrow=Nstrata, ncol=2)
  tmp.NEX.Agent2 <- matrix(NA, nrow=Nstrata, ncol=2)
  reTmp.Agent1 <- matrix(NA, nrow=Nstrata, ncol=2)
  reTmp.Agent2 <- matrix(NA, nrow=Nstrata, ncol=2)
  reTmp.eta<-matrix(NA, nrow=Nstrata, ncol=1)
  
  
  for(i1 in 1:Nstrata){
    tmp.NEX.Agent1[i1,1:2] <- mvtnorm::rmvnorm(1, sigma = diag(2))
    tmp.NEX.Agent2[i1,1:2] <- mvtnorm::rmvnorm(1, sigma = diag(2))
    reTmp.Agent1[i1,1:2] <- mvtnorm::rmvnorm(1, sigma = diag(2))
    reTmp.Agent2[i1,1:2] <- mvtnorm::rmvnorm(1, sigma = diag(2))
    reTmp.eta[i1]<-c(rnorm(1))
    
  }
  
  
  tau.norm.Agent1<-c(rnorm(1),rnorm(1))
  tau.norm.Agent2<-c(rnorm(1),rnorm(1))
  exch.index<-rep(1,Nstrata)
  
  
  
  return(list(tmp.NEX.Agent1 = tmp.NEX.Agent1,
              tmp.NEX.Agent2 = tmp.NEX.Agent2,
              reTmp.Agent1=reTmp.Agent1,
              reTmp.Agent2=reTmp.Agent2,
              tau.norm.Agent1=tau.norm.Agent1,
              tau.norm.Agent2=tau.norm.Agent2,
              exch.index=exch.index,
              reTmp.eta=reTmp.eta))
  
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
  model=combo2BLRMEXNEX,
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
DLInfo1<-lapply(1:length(data$Doses1),function(b){
  

  
  DLInfo2<-lapply(1:length(data$Doses2),function(c){
    


## loop over each probability interval
PIInfo<-lapply(1:3,function(d){
  
 probInterval<- ifelse(d==1,"UD",
         ifelse(d==2, "TT",
                ifelse(d==3, "OD",NA)))
  
  
  ##Subset to a given probability interval for a given strata and dose level
  tmp<-data.frame(round(colMeans(data.frame(fit$BUGSoutput$sims.list$Pr.Cat[,a,b,c,d])),digits=4))
  
    
    name<-strataLabels[a]
    
    tmp$Trial<-name
    
    tmp$probInterval<-probInterval
    
    tmp$dose1Level<-data$Doses1[b]
    
    tmp$dose2Level<-data$Doses2[c]
    
    return(tmp)
    
})

tmp1<-do.call("rbind",PIInfo)

return(tmp1)
})

tmp2<-do.call("rbind",DLInfo2)

return(tmp2)


})

tmp3<-do.call("rbind",DLInfo1)

return(tmp3)


})


Pcat<-do.call("rbind",out)

row.names(Pcat)<-seq(1,nrow(Pcat),1)

colnames(Pcat)<-c("prob","Trial","probInterval","Dose1","Dose2")

## create a copy for the plot
out<-Pcat

##################################################
##################################################

###################################################
## Plot the interval probabilties
###################################################






## add a reference line for EWOC and MTD (i.e., TT > 0.5) probabilities
refLine<-data.frame("probInterval"=c("TT","OD"),"ref"=c(0.5,0.25))


out$fill <- ifelse(out$prob > 0.25 & out$probInterval=="OD","firebrick1",
                              ifelse(out$prob > 0.5 & out$probInterval=="TT","chartreuse2","chartreuse4"))



## plot results, note you will get a warning due to the way fill is used to create different groups and an NA exists for the fill color if the condition doesn't exist



interval_plot_fun <- function (data, trial) {
  
  ggplot(data=data[data$Trial==trial,],aes(x=factor(probInterval),y=prob, fill=fill))+
    geom_bar(stat="identity")+ ## create a bar plot
    geom_hline(data=refLine,aes(yintercept=ref), linetype=2, color=c("black"), size=1)+ ## add reference lines
    labs(y="Interval Probabilities",x=paste0("Dose"))+ ## create lables
    ggtitle(trial)+ ## title plot
    facet_grid(Dose1~Dose2,labeller = labeller(Dose1 = setNames(paste(agent1,unique(out$Dose1),sep=":"),unique(out$Dose1)),Dose2=setNames(paste(agent2,unique(out$Dose2),sep=":"),unique(out$Dose2))))+ ## facet the results based on the probability interval
    scale_fill_identity()+ ## will the results with the conditional fill
    scale_y_continuous(breaks =c(0,0.25,0.5,0.75,1),limits=c(0,1))+ ## create axis breaks
    theme_bw()+ ## theme black and white
    theme(axis.text=element_text(size=10), ## axis font size
          axis.title=element_text(size=16,face="bold"), ## axis title font sizes
          strip.text = element_text(size = 10)) ## font size for facet
}

trial<-unique(out$Trial)



intervalPlots <- lapply(trial, interval_plot_fun, data = out)
intervalPlots[1]
 
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
 
  ## loop over the doses for a given strata
  DLInfo1<-lapply(1:length(data$Doses1),function(b){
    
    DLInfo2<-lapply(1:length(data$Doses2),function(c){

    
    tmp<-data.frame("mean"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$Pr.Tox[,a,b,c])),digits=4))
    
    tmp$SD<-sd(fit$BUGSoutput$sims.list$Pr.Tox[,a,b,c])
    
    tmp$median<-median(fit$BUGSoutput$sims.list$Pr.Tox[,a,b,c])
    
    tmp$Q2.5<-quantile(fit$BUGSoutput$sims.list$Pr.Tox[,a,b,c],probs=0.025)
    
    tmp$Q97.5<-quantile(fit$BUGSoutput$sims.list$Pr.Tox[,a,b,c],probs=0.975)
      
    tmp<-round(tmp,digits=4)
    
    

    ## Dose information
    
    tmp$Dose1<-data$Doses1[b]
    
    tmp$Dose2<-data$Doses2[c]
      

    
    ## add in strata names 
    
    name<-strataLabels[a]
    
    tmp<-cbind.data.frame("Trial"=name, tmp)
      
      return(tmp)
      

  
})


tmp1<-do.call("rbind",DLInfo2)

})
  
  tmp2<-do.call("rbind",DLInfo1)

})

P<-do.call("rbind",outProb)

row.names(P)<-seq(1,nrow(P),1)


PcatP<-merge(Pcat,P,by=c("Trial","Dose1","Dose2"))


PcatP<-PcatP[order(PcatP$Trial,PcatP$Dose1,PcatP$Dose2),]



###############################################
## Plot of the DLT probabilties
###############################################

## plot the results
postProbPlot_fun<-function(data,trial){

ggplot(data=data[data$Trial==trial,],aes(x=factor(Dose1),y=mean,color=factor(Dose2)))+
  geom_point(size=2.5,position = position_dodge(0.3))+ ## create a points for the mean DLT
  geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5),width=0.2,position = position_dodge(0.3),size=1)+ ## create bars for the 95% intervals
  labs(y="DLT Rate",x=agent1,color=agent2)+ ## labels
  ggtitle(paste0(trial,"\n Posterior Probability - Mean DLT Rate with 2.5 and 97.5 Percentiles"))+ ## title
  scale_color_brewer(palette="Set1")+ ## color points
  scale_y_continuous(breaks =c(0,0.25,0.5,0.75,1),limits=c(0,1))+ ## set y axis breaks
  theme_bw()+ ## black and white theme
  theme(legend.position="bottom", ## location of legend
        axis.text=element_text(size=14), ## axis font size
        axis.title=element_text(size=16,face="bold"), ## axis title font size 
        legend.text=element_text(size=14)) ## legend font size

}




postProbPlot <- lapply(trial, postProbPlot_fun, data = P)






###################################################
## posterior summary for log alpha
###################################################
## loop over each strata 

## agent 1
outLogAlpha.Agent1<-lapply(1:data$Nstrata,function(a){

  tmp<-data.frame("mean"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$log.alpha.Agent1[,a])),digits=4))
  
  tmp$SD<-sd(fit$BUGSoutput$sims.list$log.alpha.Agent1[,a])
  
  tmp$median<-median(fit$BUGSoutput$sims.list$log.alpha.Agent1[,a])
  
  tmp$Q2.5<-quantile(fit$BUGSoutput$sims.list$log.alpha.Agent1[,a],probs=0.025)
  
  tmp$Q97.5<-quantile(fit$BUGSoutput$sims.list$log.alpha.Agent1[,a],probs=0.975)
  
  tmp<-round(tmp,digits=4)
  

  ## add in strata names 
  
  name<-strataLabels[a]
  
  tmp<-cbind.data.frame("Trial"=name, tmp)
  
  tmp$parameter<-"logAlpha"
  
  tmp$agent<-agent1
  
  return(tmp)
  
  
  
    })
    
   


outLogAlpha.Agent1<-do.call("rbind",outLogAlpha.Agent1)


outLogAlpha.Agent2<-lapply(1:data$Nstrata,function(a){
  
  tmp<-data.frame("mean"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$log.alpha.Agent2[,a])),digits=4))
  
  tmp$SD<-sd(fit$BUGSoutput$sims.list$log.alpha.Agent2[,a])
  
  tmp$median<-median(fit$BUGSoutput$sims.list$log.alpha.Agent2[,a])
  
  tmp$Q2.5<-quantile(fit$BUGSoutput$sims.list$log.alpha.Agent2[,a],probs=0.025)
  
  tmp$Q97.5<-quantile(fit$BUGSoutput$sims.list$log.alpha.Agent2[,a],probs=0.975)
  
  tmp<-round(tmp,digits=4)
  
  
  ## add in strata names 
  
  name<-strataLabels[a]
  
  tmp<-cbind.data.frame("Trial"=name, tmp)
  
  tmp$parameter<-"logAlpha"
  
  tmp$agent<-agent2
  
  return(tmp)
  
  
  
})




outLogAlpha.Agent2<-do.call("rbind",outLogAlpha.Agent2)


###################################################
## posterior summary for log beta
###################################################
## loop over each strata 
outLogBeta.Agent1<-lapply(1:data$Nstrata,function(a){
  
  tmp<-data.frame("mean"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$log.beta.Agent1[,a])),digits=4))
  
  tmp$SD<-sd(fit$BUGSoutput$sims.list$log.beta.Agent1[,a])
  
  tmp$median<-median(fit$BUGSoutput$sims.list$log.beta.Agent1[,a])
  
  tmp$Q2.5<-quantile(fit$BUGSoutput$sims.list$log.beta.Agent1[,a],probs=0.025)
  
  tmp$Q97.5<-quantile(fit$BUGSoutput$sims.list$log.beta.Agent1[,a],probs=0.975)
  
  tmp<-round(tmp,digits=4)
  
  name<-strataLabels[a]
  
  tmp<-cbind.data.frame("Trial"=name, tmp)
  
  tmp$parameter<-"logBeta"
  
  tmp$agent<-agent1
  
  
  return(tmp)
  
  
  
})


outLogBeta.Agent1<-do.call("rbind",outLogBeta.Agent1)


## agent 2
outLogBeta.Agent2<-lapply(1:data$Nstrata,function(a){
  
  tmp<-data.frame("mean"=round(colMeans(data.frame(fit$BUGSoutput$sims.list$log.beta.Agent2[,a])),digits=4))
  
  tmp$SD<-sd(fit$BUGSoutput$sims.list$log.beta.Agent2[,a])
  
  tmp$median<-median(fit$BUGSoutput$sims.list$log.beta.Agent2[,a])
  
  tmp$Q2.5<-quantile(fit$BUGSoutput$sims.list$log.beta.Agent2[,a],probs=0.025)
  
  tmp$Q97.5<-quantile(fit$BUGSoutput$sims.list$log.beta.Agent2[,a],probs=0.975)
  
  tmp<-round(tmp,digits=4)
  
  name<-strataLabels[a]
  
  tmp<-cbind.data.frame("Trial"=name, tmp)
  
  tmp$parameter<-"logBeta"
  
  tmp$agent<-agent2
  
  
  return(tmp)
  
  
  
})


outLogBeta.Agent2<-do.call("rbind",outLogBeta.Agent2)

logAlphaBeta<-rbind.data.frame(outLogAlpha.Agent1,outLogAlpha.Agent2,outLogBeta.Agent1,outLogBeta.Agent2)

row.names(logAlphaBeta)<-seq(1,nrow(logAlphaBeta),by=1)



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


 
  
return(list(P=P,Pcat=Pcat,PcatP=PcatP,logAlphaBeta=logAlphaBeta,PmixEXNEX=PmixEXNEX,intervalPlots=intervalPlots,postProbPlot=postProbPlot))

}











