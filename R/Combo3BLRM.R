#' @title Combo3BLRM Function
#'
#' @description This function performs analysis with triple combination model
#'
#' @param Data Dose DLT data -- a list objest contains DoseAdm1, DoseAdm2, Npat, Ntox, and Weight.
#'             DoseAdm1 is a vector of first agent doses which are administered at patient cohorts. 
#'             DoseAdm2 is a vector of second agent doses which are administered at patient cohorts. 
#'             DoseAdm3 is a vector of third agent doses which are administered at patient cohorts. 
#'             Npat is a vector of same length as DoseAdm, the elements are the cohort sizes at each DoseAdm. 
#'             Ntox is a vector of same length as DoseAdm, the elements are the the number of toxicitites 
#'             observed for each cohort. Weight is a vector of same length as DoseAdm, the elements are 1 
#'             by default, but can be changed depending on how much information one wants borrow form this 
#'             data set; 1 is maximum i.e. borrowing the full information, 0 is minimum i.e. borrowing 
#'             no information, similary, setting 0.5 means borrowing half of the information
#' @param Agent1 Compound 1
#' @param Doses1 Doses for compound 1
#' @param DoseRef1 Reference dose for compound 1
#' @param Prior.logAB.1 Prior for single agent effect of compound 1 (a list object)
#' @param pmix.logAB1 mix probability compound 1
#' @param Agent2 Compound 2
#' @param Doses2 Doses for Compound 2
#' @param DoseRef2 Dose reference for compound 2
#' @param Prior.logAB.2 Prior for single agent effect of compound 2 (a list object)
#' @param pmix.logAB2 mix probability compound 2
#' @param Agent3 Compound 3
#' @param Doses3 Doses for compound 3
#' @param DoseRef3 Reference dose for compound 3
#' @param Prior.logAB.3 Prior for single agent effect of compound 3 (a list object)
#' @param pmix.logAB3 mix probability compound 3
#' @param Prior.eta Prior for interaction parameters (a list object with \eqn{\eta_{12}}, \eqn{\eta_{13}},
#'        \eqn{\eta_{23}}, and \eqn{\eta_{123}})
#' @param Prior.Only PriorOnly
#' @param New.cohortsize new for prediction
#' @param Outfile output file
#' @param Plot Indicator if Plot will be generated
#' @param int.crit Critical value for interval to change color
#' @param int.col Color of each bar
#' @param pars MCMC parameters to save
#' @param MCMC MCMC specifications (nburn, ntot, nchain, nthin)
#' @param Warn Prints convergence warning
#' @param DIC DIC to save
#' @param progress.bar Jags parameter for type of progress bar. Possible values are
#' “text”, “gui”, and “none”. Type “text” is displayed on the R console. Type “gui”
#' is a graphical progress bar in a new window. The progress bar is suppressed if
#' progress.bar is “none”
#' @param Rnd.seed seed for R
#' @param RSeed.WinBUGS seed for WinBugs
#'
#' @details See vignettes. browseVignettes("OncoPh1BLRM") 
#' 
#' @references Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. A Bayesian Industry
#' Approach to Phase I Combination Trials in Oncology. Statistical Methods for Drug Combination
#' Studies. Taylor & Francis, 2015.
#' 
#' 
#' @return \item{Data}{The dataset}
#' \item{Priors}{Priors used for the analysis}
#' \item{logAlphaBeta1}{Posterior summary of log(\eqn{\alpha}) and log(\eqn{\beta}) for agent 1}
#' \item{logAlphaBeta.cor1}{Posterior summary of the correlation parameter for agent 1}
#' \item{logAlphaEBeta1}{Posterior summary of log(\eqn{\alpha}) and \eqn{\beta} for agent 1}
#' \item{logAlphaBeta2}{Posterior summary of log(\eqn{\alpha}) and log(\eqn{\beta}) for agent 2}
#' \item{logAlphaBeta.cor2}{Posterior summary of the correlation parameter for agent 2}
#' \item{logAlphaEBeta2}{Posterior summary of log(\eqn{\alpha}) and \eqn{\beta} for agent 2}
#' \item{logAlphaBeta3}{Posterior summary of log(\eqn{\alpha}) and log(\eqn{\beta}) for agent 3}
#' \item{logAlphaBeta.cor3}{Posterior summary of the correlation parameter for agent 3}
#' \item{logAlphaEBeta3}{Posterior summary of log(\eqn{\alpha}) and \eqn{\beta} for agent 3}
#' \item{eta12}{Posterior summary of intercation paramter \eqn{\eta_{12}})}
#' \item{eta23}{Posterior summary of intercation paramter \eqn{\eta_{23}})}
#' \item{eta13}{Posterior summary of intercation paramter \eqn{\eta_{13}})}
#' \item{eta123}{Posterior summary of intercation paramter \eqn{\eta_{123}})}
#' \item{OddsFactor12}{Posterior summary of the odds factor between first and second agent}
#' \item{OddsFactor23}{Posterior summary of the odds factor between second and third agent}
#' \item{OddsFactor13}{Posterior summary of the odds factor between first and third agent}
#' \item{OddsFactor123}{Posterior summary of the odds factor between first, second, and third agent}
#' \item{P}{Predicted summaries of the dose levels}
#' \item{Pcat}{Posterior predictive probabilities of doses}
#' \item{PcatP}{Posterior predictive probabilities and summaries of doses}
#' \item{pred.Tox}{Predictive toxicity probilities of doses}
#' \item{Post.pmix1}{Posterior mixture weigths for agent 1}
#' \item{Post.pmix2}{Posterior mixture weigths for agent 2}
#' \item{Post.pmix3}{Posterior mixture weigths for agent 3}
#' \item{Rhat}{The MCMC convergence dignostoc Rhat for all the parameters}
#' \item{MCMC}{MCMC informations}
#' \item{R2WB}{}
#' In addition, a text output containing all of the above and two plots -- posterior DLT rates and 
#' interval probabilities of DLT's, are produced. 
#' 
#' 
#' @examples  
#' \dontrun{
#' 
#' 
#' # The following example is a dose escalation study 
#' # where the escalation is done with respect to the first agent doses
#' # The study start at the maximum tolerated dose (MTD) found by the dose escalation study of the first agent
#' # along with the prespecified fixed doses of the other two agents
#'
#' # Form Interaction prior
#' prior.interaction <- function(Odds.median = NULL, Odds.upper = NULL){
#'  mean.eta <- log(Odds.median)
#'  sd.eta <- (log(Odds.upper) - log(Odds.median))/1.96
#'  return(list(mean.eta=mean.eta, sd.eta=sd.eta))
#' }
#'
#' # Prior for eta: Median -- No increase in Odds of DLT, Upper -- 9 fold increase in Odds of DLT
#' peta <- prior.interaction(Odds.median = 1, Odds.upper = 9)
#' etaPrior <- c(peta$mean.eta, peta$sd.eta)
#' etaPrior
#'
#' eta.prior = list(etaPrior12  = etaPrior,
#'                  etaPrior13  = etaPrior,
#'                  etaPrior23  = etaPrior,
#'                  etaPrior123 = etaPrior
#')
#'
#'
#'
#' Doses.Agent1 <- c(8, 16, 32, 64, 128, 256, 512)  # 8 is starting dose, 50 is an intermediate dose to be investigated
#' # 500 is not tolerable
#' DoseRef.Agent1 <- 275  # Expected MTD 50 + (50 + 500)/2
#' MTD.Agent1     <- 128   # However, this is delcared MTD after the dose escalation experiment
#'
#'
#' require(boot)
#' weakly.informative.prior <- c(logit(0.33), 0, 2, 1, 0)  # Default
#' 
#' 
#' # Doses.Agent2   <- c(0.02, 1, 2.5, 5, 10, 30)  # Used for monotherpy
#' DoseRef.Agent2 <- 2.5
#' 
#' # Doses.Agent3   <- c(25, 50, 75, 100, 125, 150)  # Used for monotherpy
#' DoseRef.Agent3 <- 125
#' 
#' Doses.Agent1 <- c(Doses.Agent1, MTD.KaT6)  
#' # For this agent the data are dose escalation data and the MTD
#' Doses.Agent2 <- c(rep(0, length(Doses.Agent1)), DoseRef.Letro) 
#' # For this agent the data are 0's and the prespecified fixed dose
#' Doses.Agent3 <- c(rep(0, length(Doses.Agent1)), DoseRef.Palbo) 
#' # For this agent the data are 0's and the prespecified fixed dose
#'
#'
#' # Triple agent data
#' Data.combo3 <- list(DosesAdm1 = MTD.Agent1,
#'                     DosesAdm2 = DoseRef.Agent2,
#'                     DosesAdm3 = DoseRef.Agent3,
#'                     Npat      = c(6),
#'                     Ntox      = c(0),
#'                     Weight    = c(1)
#' )
#'
#'
#'
#' fit.triple <- Combo3BLRM(Data           = Data.combo3, 
#'                          Agent1         = "Agent1", 
#'                          Doses1         = Doses.Agent1,
#'                          DoseRef1       = DoseRef.Agent1, 
#'                          Prior.logAB.1  = list(weakly.informative.prior), 
#'                          pmix.logAB1    = c(1),
#'                          Agent2         = "Agent2", 
#'                          Doses2         = Doses.Agent2, 
#'                          DoseRef2       = DoseRef.Agent2,
#'                          Prior.logAB.2  = list(weakly.informative.prior), 
#'                          pmix.logAB2    = c(1), 
#'                          Agent3         = "Agent3",
#'                          Doses3         = Doses.Agent3, 
#'                          DoseRef3       = DoseRef.Agent3, 
#'                          Prior.logAB.3  = list(weakly.informative.prior),
#'                          pmix.logAB3    = c(1), 
#'                          Prior.eta      = eta.prior, 
#'                          Prior.Only     = FALSE,
#'                          Pcutoffs       = c(0.16, 0.33), 
#'                          New.cohortsize = 3, 
#'                          Outfile        = "Triple_0_6_128",
#'                          Plot           = TRUE, 
#'                          int.crit       = c(1, 1, 0.25), 
#'                          int.col        = list("lightgreen", "green", c("yellow", "red")), 
#'                          pars           = c("logAlphaBeta1", "logAlphaBeta2", "logAlphaBeta3", 
#'                                             "eta12", "eta13", "eta23", "eta123", "P12", "pCat", "pred.Tox"), 
#'                          MCMC           = c(20000, 30000, 4, 1), 
#'                          Warn           = TRUE,
#'                          DIC            = TRUE, 
#'                          Rnd.seed       = 1452671,
#'                          RSeed.WinBUGS  = 1367)
#' 
#' fit.triple$PcatP[,, 2, 2, ]
#' 
#'}
#' 
#' 
#' @export


Combo3BLRM <- function(Data               = NULL,
                       Agent1             = NULL,
                       Doses1             = NULL,
                       DoseRef1           = NULL,
                       Prior.logAB.1      = NULL,
                       pmix.logAB1        = c(1),
                       Agent2             = NULL,
                       Doses2             = NULL,
                       DoseRef2           = NULL,
                       Prior.logAB.2      = NULL,
                       pmix.logAB2        = c(1),
                       Agent3             = NULL,
                       Doses3             = NULL,
                       DoseRef3           = NULL,
                       Prior.logAB.3      = NULL,
                       pmix.logAB3        = c(1),
                       Prior.eta          = NULL,
                       Prior.Only         = FALSE,
                       Pcutoffs           = c(0.16, 0.33),
                       New.cohortsize     = 3,
                       Outfile            = "Out3",
                       Plot               = TRUE,
                       int.crit           = c(1,1,0.25),
                       int.col            = list( "lightgreen", "green", c("yellow","red")),
                       pars               = c("logAlphaBeta1", "logAlphaBeta2", "logAlphaBeta3", "eta12", "eta13","eta23", "eta123", "P12","pCat","pred.Tox"),
                       MCMC               = c(2500, 7500, 4, 1),
                       Warn               = TRUE,
                       DIC                = TRUE,
                       progress.bar       = progress.bar,
                       Rnd.seed           = 1452671,
                       RSeed.WinBUGS      = 1367
                      )
{

  version = "30-November-2017"
  
  stopifnot("Mixture weights must add to 1" = sum(pmix.logAB1) == 1)
  stopifnot("Mixture weights must add to 1" = sum(pmix.logAB2) == 1)
  stopifnot("Mixture weights must add to 1" = sum(pmix.logAB3) == 1)

  set.seed(Rnd.seed)

  if(Prior.Only){DIC=FALSE}

  # MCMC parameters
  n.burnin = MCMC[1]
  n.iter = MCMC[2]
  n.chains = MCMC[3]
  n.thin = MCMC[4]

  #Reading Data

if(!Prior.Only){

  DosesAdm1  = Data$DosesAdm1
  DosesAdm2  = Data$DosesAdm2
  DosesAdm3  = Data$DosesAdm3
  Ntox       = Data$Ntox*Data$Weight
  Npat       = Data$Npat*Data$Weight

  DosesAdm1[DosesAdm1 ==0] = 0.00001
  DosesAdm2[DosesAdm2 ==0] = 0.00001
  DosesAdm3[DosesAdm3 ==0] = 0.00001

  Ncohorts <- length(DosesAdm1)

  DosesAdm1 <- cbind(NULL, DosesAdm1)
  DosesAdm2 <- cbind(NULL, DosesAdm2)
  DosesAdm3 <- cbind(NULL, DosesAdm3)
  Ntox <- cbind(NULL, Ntox)
  Npat <- cbind(NULL, Npat)
  zero <- rep(0,Ncohorts)

                }

if(Prior.Only){Npat=0; Ntox=0}

  Ndoses1 <- length(Doses1)
  Ndoses2 <- length(Doses2)
  Ndoses3 <- length(Doses3)


  Ncat <- length(Pcutoffs) + 1
  Pcutoffs1 = c(0, Pcutoffs, 1)

  intLabels = paste(Pcutoffs1[1:Ncat], "-", Pcutoffs1[2:(Ncat + 1)], sep = "")
  Doses1.Label = paste(Agent1, "=", Doses1, sep = '')
  Doses2.Label = paste(Agent2, "=", Doses2, sep = '')
  Doses3.Label = paste(Agent3, "=", Doses3, sep = '')
  pred.tox.Label = paste("DLT=",0:New.cohortsize, sep='')

  Doses1.1 =Doses1
  Doses2.1 = Doses2
  Doses3.1 = Doses3

  Doses1[Doses1 ==0] = 0.00001
  Doses2[Doses2 ==0] = 0.00001
  Doses3[Doses3 ==0] = 0.00001

  Doses1 <- cbind(NULL, Doses1)
  Doses2 <- cbind(NULL, Doses2)
  Doses3 <- cbind(NULL, Doses3)

  #Priors

  #Bivariate normal prior for (logAlpha1,logBeta1) [marginal model for Agent1]

  Prior.logAB1 <- c(Prior.logAB.1, list(c(0,0,0.1,0.1,0)))
  Pmix1 <- c(pmix.logAB1, 0)


  for(j in 1:length(Prior.logAB1)){

         Prior1 <- do.call(rbind,Prior.logAB1)

  }

  Pmix1 <- cbind(NULL, Pmix1)
  Nmix1 <- length(Prior.logAB1)


  #Bivariate normal prior for (logAlpha2,logBeta2) [marginal model for Agent2]

  Prior.logAB2 <- c(Prior.logAB.2, list(c(0,0,0.1,0.1,0)))
  Pmix2 <- c(pmix.logAB2, 0)

  for(j in 1:length(Prior.logAB2)){

        Prior2 <- do.call(rbind,Prior.logAB2)

  }

  Pmix2 <- cbind(NULL,Pmix2)
  Nmix2 <- length(Prior.logAB2)

  #Bivariate normal prior for (logAlpha3,logBeta3) [marginal model for Agent3]

  Prior.logAB3 <- c(Prior.logAB.3, list(c(0,0,0.1,0.1,0)))
  Pmix3 <- c(pmix.logAB3, 0)


  for(j in 1:length(Prior.logAB3)){

       Prior3 <- do.call(rbind,Prior.logAB3)

                                  }
  Pmix3 <- cbind(NULL,Pmix3)
  Nmix3 <- length(Prior.logAB3)


  #Eta priors

  etaPrior12  = Prior.eta$etaPrior12
  etaPrior13  = Prior.eta$etaPrior13
  etaPrior23  = Prior.eta$etaPrior23
  etaPrior123 = Prior.eta$etaPrior123



    #Choice of model

  if(Prior.Only){model= Combo3BLRM.Prior.WB; model.index=1}
  if(!Prior.Only & (any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)))){model= Combo3BLRM.Weight.WB; model.index=2}
  if(!Prior.Only & !(any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)))){model= Combo3BLRM.WB; model.index=3}

  #Data for jags
  if(Prior.Only){

    data= list("Ndoses1", "Ndoses2", "Ndoses3", "Ncat", "Nmix1", "Nmix2", "Nmix3",
               "Prior1", "Prior2", "Prior3", "Pmix1", "Pmix2", "Pmix3",
               "etaPrior12", "etaPrior13", "etaPrior23", "etaPrior123",
               "Doses1", "Doses2", "Doses3", "DoseRef1", "DoseRef2", "DoseRef3",
               "Pcutoffs", "New.cohortsize")
  }
  if(!Prior.Only & model.index==3){

    data= list("Ncohorts", "Ndoses1", "Ndoses2", "Ndoses3", "Ncat", "Nmix1", "Nmix2", "Nmix3",
               "Prior1", "Prior2", "Prior3", "Pmix1", "Pmix2", "Pmix3",
               "etaPrior12", "etaPrior13", "etaPrior23", "etaPrior123",
               "Doses1", "Doses2", "Doses3", "DoseRef1", "DoseRef2", "DoseRef3",
               "DosesAdm1", "DosesAdm2", "DosesAdm3", "Ntox", "Npat",
               "Pcutoffs", "New.cohortsize")
  }

  if(!Prior.Only & model.index==2){

    data= list("Ncohorts", "Ndoses1", "Ndoses2", "Ndoses3", "Ncat", "Nmix1", "Nmix2", "Nmix3",
               "Prior1", "Prior2", "Prior3", "Pmix1", "Pmix2", "Pmix3",
               "etaPrior12", "etaPrior13", "etaPrior23", "etaPrior123",
               "Doses1", "Doses2", "Doses3", "DoseRef1", "DoseRef2", "DoseRef3",
               "DosesAdm1", "DosesAdm2", "DosesAdm3", "Ntox", "Npat", "zero",
               "Pcutoffs", "New.cohortsize")
  }



  #Initial values
  inits.fun = function(i){

    logAlphaBeta1All <- matrix(NA, nrow=Nmix1, ncol=2)
    logAlphaBeta2All <- matrix(NA, nrow=Nmix2, ncol=2)
    logAlphaBeta3All <- matrix(NA, nrow=Nmix3, ncol=2)

    for(i1 in 1:Nmix1){logAlphaBeta1All[i1,1:2] <- rmvnorm(1, sigma = diag(2))}
    for(i2 in 1:Nmix2){logAlphaBeta2All[i2,1:2] <- rmvnorm(1, sigma = diag(2))}
    for(i3 in 1:Nmix3){logAlphaBeta3All[i3,1:2] <- rmvnorm(1, sigma = diag(2))}

    eta12 <- rnorm(1, 0, 0.1)
    eta13 <- rnorm(1, 0, 0.1)
    eta23 <- rnorm(1, 0, 0.1)
    eta123 <- rnorm(1, 0, 0.1)

    Z1 <- 1
    Z2 <- 1
    Z3 <- 1

    return(list(logAlphaBeta1All = logAlphaBeta1All,
                logAlphaBeta2All = logAlphaBeta2All,
                logAlphaBeta3All = logAlphaBeta3All,
                eta12 = eta12,
                eta13 = eta13,
                eta23 = eta23,
                eta123 = eta123,
                Z1=Z1,
                Z2=Z2,
                Z3=Z3)
           )
  }

  inits <- lapply(rep(1,n.chains),inits.fun)

  RndSeed.WinBUGS <- RSeed.WinBUGS

  pars = c(pars,"logAlphaEBeta1", "logAlphaEBeta2", "logAlphaEBeta3","logAlphaBeta1.XY", "logAlphaBeta2.XY", "logAlphaBeta3.XY", "OddsFactor12", "OddsFactor13", "OddsFactor23", "OddsFactor123")


  fit = jags(
    data=data,
    inits= inits,
    parameters.to.save = pars,
    model=model,
    n.chains=n.chains,n.burnin=n.burnin,n.iter=n.iter,n.thin=n.thin,
    jags.seed=RndSeed.WinBUGS,
    DIC= DIC,
    progress.bar =progress.bar
    #bugs.directory="C:/Users/roychs04/Documents/Folder/Software/WinBUGS14"
  )


  #Processing output

  vnames = c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat")
  summary = fit$BUGSoutput$summary
  #fit$BUGSoutput$sims.matrix = NULL
  fit$BUGSoutput$sims.array = NULL

  if(Warn){

  logAlphaBeta1.neff <- fit$BUGSoutput$summary[c("logAlphaBeta1[1]", "logAlphaBeta1[2]"), "n.eff"]
  if (any(logAlphaBeta1.neff  < 1000)) {
    cat("\nWARNING: Neff < 1000 for logAlphaBeta1  Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for logAlphaBeta1 Consider increasing the sample size (number of iterations).")
  }


  logAlphaBeta2.neff <- fit$BUGSoutput$summary[c("logAlphaBeta2[1]", "logAlphaBeta2[2]"), "n.eff"]
  if (any(logAlphaBeta2.neff  < 1000)) {
    cat("\nWARNING: Neff < 1000 for logAlphaBeta2  Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for logAlphaBeta2 Consider increasing the sample size (number of iterations).")
  }

  logAlphaBeta3.neff <- fit$BUGSoutput$summary[c("logAlphaBeta3[1]", "logAlphaBeta3[2]"), "n.eff"]
  if (any(logAlphaBeta3.neff  < 1000)) {
    cat("\nWARNING: Neff < 1000 for logAlphaBeta3  Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for logAlphaBeta3 Consider increasing the sample size (number of iterations).")
  }


  eta12.neff <- fit$BUGSoutput$summary[c("eta12"), "n.eff"]
  if (any(eta12.neff  < 1000)) {
    cat("\nWARNING: Neff < 1000 for eta12  Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for eta12 Consider increasing the sample size (number of iterations).")
  }

  eta13.neff <- fit$BUGSoutput$summary[c("eta13"), "n.eff"]
  if (any(eta13.neff  < 1000)) {
    cat("\nWARNING: Neff < 1000 for eta13  Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for eta13 Consider increasing the sample size (number of iterations).")
  }

  eta23.neff <- fit$BUGSoutput$summary[c("eta23"), "n.eff"]
  if (any(eta23.neff  < 1000)) {
    cat("\nWARNING: Neff < 1000 for eta23  Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for eta23 Consider increasing the sample size (number of iterations).")
  }

  eta123.neff <- fit$BUGSoutput$summary[c("eta123"), "n.eff"]
  if (any(eta123.neff  < 1000)) {
    cat("\nWARNING: Neff < 1000 for eta123  Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for eta123 Consider increasing the sample size (number of iterations).")
  }
  }

  logAlphaBeta1 = NULL
  Rhat.logAlphaBeta1 = NULL
  logAlphaBeta1.cor=NULL
  Rhat.logAlphaBeta1.cor=NULL
  logAlphaEBeta1 =NULL
  Rhat.logAlphaEBeta1 =NULL


  logAlphaBeta2 = NULL
  Rhat.logAlphaBeta2 = NULL
  logAlphaBeta2.cor=NULL
  Rhat.logAlphaBeta2.cor=NULL
  logAlphaEBeta2 =NULL
  Rhat.logAlphaEBeta2 =NULL



  logAlphaBeta3 = NULL
  Rhat.logAlphaBeta3 = NULL
  logAlphaBeta3.cor=NULL
  Rhat.logAlphaBeta3.cor=NULL
  logAlphaEBeta3 =NULL
  Rhat.logAlphaEBeta3 =NULL

  eta12 = NULL
  Rhat.eta12 = NULL
  OddsFactor12 =NULL
  Rhat.OddsFactor12=NULL
  eta13 = NULL
  Rhat.eta13 = NULL
  OddsFactor13 =NULL
  Rhat.OddsFactor13=NULL
  eta23 = NULL
  Rhat.eta23 = NULL
  OddsFactor23 =NULL
  Rhat.OddsFactor23=NULL
  eta123 = NULL
  Rhat.eta123 =NULL
  OddsFactor123 =NULL
  Rhat.OddsFactor123=NULL

  P = NULL
  Rhat.P = NULL

  Pcat = NULL
  PcatP = NULL
  pred.Tox =NULL

  Post.pmix1 =NULL
  Post.pmix2 =NULL
  Post.pmix3 = NULL


  if(is.element("logAlphaBeta1", pars)){
    logAlphaBeta1 = BUGS.Select.Output("logAlphaBeta1", summary, cols = vnames)
    Rhat.logAlphaBeta1 = logAlphaBeta1[, "Rhat"]

    logAlphaBeta1 = BUGS.Table2Array(logAlphaBeta1, Labels = list(c("logAlpha1", "logBeta1")), Dim = 2)

    logAlphaEBeta1 = BUGS.Select.Output("logAlphaEBeta1", summary, cols = vnames)
    Rhat.logAlphaEBeta1 = logAlphaEBeta1[, "Rhat"]

    logAlphaEBeta1 = BUGS.Table2Array(logAlphaEBeta1, Labels = list(c("logAlpha1", "beta1")), Dim = 2)

  }

  if(is.element("logAlphaBeta1.XY", pars)){

    logAlphaBeta1.XY = BUGS.Select.Output("logAlphaBeta1.XY", summary, cols = vnames)
    Rhat.logAlphaBeta1.XY.cor = logAlphaBeta1.XY["Rhat"]

    logAlphaBeta1.cor <-  (logAlphaBeta1.XY[1] - logAlphaBeta1[1, 1]*logAlphaBeta1[2,1])/logAlphaBeta1[1, 2]/logAlphaBeta1[2, 2]

  }


  if(is.element("logAlphaBeta2", pars)){
    logAlphaBeta2 = BUGS.Select.Output("logAlphaBeta2", summary, cols = vnames)
    Rhat.logAlphaBeta2 = logAlphaBeta1[, "Rhat"]

    logAlphaBeta2 = BUGS.Table2Array(logAlphaBeta2, Labels = list(c("logAlpha2", "logBeta2")), Dim = 2)

    logAlphaEBeta2 = BUGS.Select.Output("logAlphaEBeta2", summary, cols = vnames)
    Rhat.logAlphaEBeta2 = logAlphaEBeta2[, "Rhat"]

    logAlphaEBeta2 = BUGS.Table2Array(logAlphaEBeta2, Labels = list(c("logAlpha2", "beta2")), Dim = 2)

  }

  if(is.element("logAlphaBeta2.XY", pars)){

    logAlphaBeta2.XY = BUGS.Select.Output("logAlphaBeta2.XY", summary, cols = vnames)
    Rhat.logAlphaBeta2.XY.cor = logAlphaBeta2.XY["Rhat"]

    logAlphaBeta2.cor <-  (logAlphaBeta2.XY[1] - logAlphaBeta2[1, 1]*logAlphaBeta2[2,1])/logAlphaBeta2[1, 2]/logAlphaBeta2[2, 2]

  }


  if(is.element("logAlphaBeta3", pars)){
    logAlphaBeta3 = BUGS.Select.Output("logAlphaBeta3", summary, cols = vnames)
    Rhat.logAlphaBeta3 = logAlphaBeta3[, "Rhat"]

    logAlphaBeta3 = BUGS.Table2Array(logAlphaBeta3, Labels = list(c("logAlpha3", "logBeta3")), Dim = 2)

    logAlphaEBeta3 = BUGS.Select.Output("logAlphaEBeta3", summary, cols = vnames)
    Rhat.logAlphaEBeta3 = logAlphaEBeta2[, "Rhat"]

    logAlphaEBeta3 = BUGS.Table2Array(logAlphaEBeta3, Labels = list(c("logAlpha3", "beta3")), Dim = 2)

  }


  if(is.element("logAlphaBeta3.XY", pars)){

    logAlphaBeta3.XY = BUGS.Select.Output("logAlphaBeta3.XY", summary, cols = vnames)
    Rhat.logAlphaBeta3.XY.cor = logAlphaBeta3.XY["Rhat"]

    logAlphaBeta3.cor <-  (logAlphaBeta3.XY[1] - logAlphaBeta3[1, 1]*logAlphaBeta3[2,1])/logAlphaBeta3[1, 2]/logAlphaBeta3[2, 2]

  }


  if(is.element("eta12", pars)){
    eta12 = BUGS.Select.Output("eta12", summary, cols = vnames)
    Rhat.eta12 = eta12["Rhat"]
    eta12 <- rbind(NULL, eta12)

    rownames(eta12) = "eta12"

  }

  if(is.element("OddsFactor12", pars)){
    OddsFactor12 = BUGS.Select.Output("OddsFactor12", summary, cols = vnames)
    Rhat.OddsFactor12 = OddsFactor12["Rhat"]
    OddsFactor12 <- rbind(NULL, OddsFactor12)

    rownames(OddsFactor12) = "OddsFactor12"

  }

  if(is.element("eta13", pars)){
    eta13 = BUGS.Select.Output("eta13", summary, cols = vnames)
    Rhat.eta13 = eta13["Rhat"]
    eta13 <- rbind(NULL, eta13)

    rownames(eta13) = "eta13"

  }

  if(is.element("OddsFactor13", pars)){
    OddsFactor13 = BUGS.Select.Output("OddsFactor13", summary, cols = vnames)
    Rhat.OddsFactor13 = OddsFactor13["Rhat"]
    OddsFactor13 <- rbind(NULL, OddsFactor13)

    rownames(OddsFactor13) = "OddsFactor13"

  }

  if(is.element("eta23", pars)){
    eta23 = BUGS.Select.Output("eta23", summary, cols = vnames)
    Rhat.eta23 = eta23["Rhat"]
    eta23 <- rbind(NULL, eta23)

    rownames(eta23) = "eta23"

  }

  if(is.element("OddsFactor23", pars)){
    OddsFactor23 = BUGS.Select.Output("OddsFactor23", summary, cols = vnames)
    Rhat.OddsFactor23 = OddsFactor23["Rhat"]
    OddsFactor23 <- rbind(NULL, OddsFactor23)

    rownames(OddsFactor23) = "OddsFactor23"

  }


  if(is.element("eta123", pars)){
    eta123 = BUGS.Select.Output("eta123", summary, cols = vnames)
    Rhat.eta123 = eta123["Rhat"]
    eta123 <- rbind(NULL, eta123)

    rownames(eta123) = "eta123"

  }

  if(is.element("OddsFactor123", pars)){
    OddsFactor123 = BUGS.Select.Output("OddsFactor123", summary, cols = vnames)
    Rhat.OddsFactor123 = OddsFactor123["Rhat"]
    OddsFactor123 <- rbind(NULL, OddsFactor123)

    rownames(OddsFactor123) = "OddsFactor123"

  }




  if (is.element("P12", pars)) {
    P = BUGS.Select.Output("P12", summary, cols = vnames)
    Rhat.P = P[, "Rhat"]
    aiw = round(Beta.ms2ab(P[, "mean"], P[, "sd"])$n, 1)
    P = cbind(P, aiw)
    P = BUGS.Table2Array(P, Labels = list(Doses1.Label, Doses2.Label, Doses3.Label), Dim = c(Ndoses1, Ndoses2, Ndoses3))
    P = aperm(P, c(1,4,2,3))

  }


  if (is.element("pCat", pars)) {
    Pcat = BUGS.Select.Output("pCat", summary, cols = vnames)
    Pcat = BUGS.Table2Array(Pcat, Labels = list(Doses1.Label, Doses2.Label, Doses3.Label, intLabels), Dim = c(Ndoses1, Ndoses2, Ndoses3, Ncat))
    Pcat = Pcat[, , , ,"mean", drop = FALSE]
    Pcat = aperm(Pcat, c(1, 4, 2, 3,5))

  }

 if (is.element("pCat", pars) & is.element("P12", pars)) {

     PcatP <- array(NA, dim=c(Ndoses1, 10, Ndoses2, Ndoses3, 1))

     for(k1 in 1:Ndoses2){
       for(k2 in 1:Ndoses3){

         PcatP[ ,,k1,k2,1] <- cbind(rbind(NULL,Pcat[,,k1,k2,1]), rbind(NULL,P[,,k1,k2]))

                          }
     }

     vnames1 = c(intLabels,vnames,"aiw")

     dimnames(PcatP) = list(Doses1.Label, vnames1, Doses2.Label, Doses3.Label, "")
 }


  if (is.element("pred.Tox", pars)){

    pred.Tox = BUGS.Select.Output("pred.Tox", summary, cols = vnames)
    Rhat.pred.Tox= pred.Tox[, "Rhat"]



    pred.Tox = BUGS.Table2Array(pred.Tox, Labels = list(Doses1.Label, Doses2.Label, Doses3.Label,pred.tox.Label), Dim = c(Ndoses1, Ndoses2, Ndoses3, New.cohortsize+1))
    pred.Tox = pred.Tox[, , , ,"mean", drop = FALSE]
    pred.Tox = aperm(pred.Tox, c(1,4, 2, 3,5))

  }

  if (is.element("Zcat1", pars)){

    Post.pmix1.1 = BUGS.Select.Output("Zcat1", summary, cols = "mean")
    Post.pmix1 <- cbind(NULL,Post.pmix1.1[-Nmix1])
    Post.pmix1.nms <- paste("mix1",1:(Nmix1-1),sep='')
    rownames(Post.pmix1) <- Post.pmix1.nms

  }


  if (is.element("Zcat2", pars)){

    Post.pmix2.1 = BUGS.Select.Output("Zcat2", summary, cols = "mean")
    Post.pmix2 <- cbind(NULL,Post.pmix2.1[-Nmix2])
    Post.pmix2.nms <- paste("mix2",1:(Nmix2-1),sep='')
    rownames(Post.pmix2) <- Post.pmix2.nms

  }

  if (is.element("Zcat3", pars)){

    Post.pmix3.1 = BUGS.Select.Output("Zcat3", summary, cols = "mean")
    Post.pmix3 <- cbind(NULL,Post.pmix3.1[-Nmix3])
    Post.pmix3.nms <- paste("mix3",1:(Nmix3-1),sep='')
    rownames(Post.pmix3) <- Post.pmix3.nms

  }

   Rhat = c(Rhat.logAlphaBeta1, Rhat.logAlphaBeta1.cor, Rhat.logAlphaEBeta1,
            Rhat.logAlphaBeta2, Rhat.logAlphaBeta2.cor, Rhat.logAlphaEBeta2,
            Rhat.logAlphaBeta3, Rhat.logAlphaBeta3.cor, Rhat.logAlphaEBeta3,
            Rhat.eta12, Rhat.eta13, Rhat.eta23, Rhat.eta123,
            Rhat.P)

   Rhat.max = max(Rhat)

 if (Rhat.max > 1.1) {
    warning("\n Warning: at least one of the parameters has a Gelman-Rubin diagnostic Rhat > 1.1")
 }

   Priors = list(Prior.logAB1 = Prior.logAB.1,
                 pmix.logAB1 = pmix.logAB1,
                 Prior.logAB2 = Prior.logAB.2,
                 pmix.logAB2 = pmix.logAB2,
                 Prior.logAB3 = Prior.logAB.3,
                 pmix.logAB3 = pmix.logAB3,
                 etaPrior12 = etaPrior12,
                 etaPrior13 = etaPrior13,
                 etaPrior23 = etaPrior23,
                 etaPrior123= etaPrior123)

   labels = list(intLabels = intLabels)
   Agents = c(Agent1, Agent2, Agent3)

   if(!Prior.Only){
   Data.1 = data.frame(t(do.call(rbind, Data)))
   colnames(Data.1) = c("DoseAdm1", "DoseAdm2","DoseAdm3","Npat","Ntox", "Weight")
   N = c(Ncohorts = Ncohorts, Ndoses1 = Ndoses1, Ndoses2=Ndoses2, Ndoses3=Ndoses3, Nmix1=Nmix1-1, Nmix2=Nmix2-1, Nmix3=Nmix3-1)

   Data11 <- aggregate(Data.1, by=list(Data.1$DoseAdm1, Data.1$DoseAdm2, Data.1$DoseAdm3, Data.1$Weight),sum)[,c(1,2,3,4,8,9)]
   Data1 <- Data11[,c(1,2,3,5,6,4)]
   colnames(Data1) = c("DoseAdm1", "DoseAdm2", "DoseAdm3","Npat","Ntox", "Weight")


   }

   if(Prior.Only){
     Data1 = NULL
     N = c(Ndoses1 = Ndoses1, Ndoses2=Ndoses2, Ndoses3=Ndoses3, Nmix1=Nmix1-1, Nmix2=Nmix2-1, Nmix3=Nmix3-1)
   }

   #Collecting outputs

   outlist = list(Data = list(data = Data1, N = N, labels = labels, Agents = Agents),
                  Priors = Priors,
                  logAlphaBeta1 = logAlphaBeta1,
                  logAlphaBeta1.cor=logAlphaBeta1.cor,
                  logAlphaEBeta1 = logAlphaEBeta1,
                  logAlphaBeta2 = logAlphaBeta2,
                  logAlphaBeta2.cor= logAlphaBeta2.cor,
                  logAlphaEBeta2 = logAlphaEBeta2,
                  logAlphaBeta3 = logAlphaBeta3,
                  logAlphaBeta3.cor= logAlphaBeta3.cor,
                  logAlphaEBeta3 = logAlphaEBeta3,
                  eta12 = eta12,
                  eta13 = eta13,
                  eta23 = eta23,
                  eta123 = eta123,
                  OddsFactor12 = OddsFactor12,
                  OddsFactor13 = OddsFactor13,
                  OddsFactor23 = OddsFactor23,
                  OddsFactor123 = OddsFactor123,
                  P = P,
                  Pcat = Pcat,
                  PcatP = PcatP,
                  pred.Tox = pred.Tox,
                  Post.pmix1 =Post.pmix1,
                  Post.pmix2 =Post.pmix2,
                  Post.pmix3 = Post.pmix3,
                  Rhat = list(Rhat = Rhat, Rhat.max = Rhat.max),
                  MCMC = list(pars.monitored = pars, MCMC = MCMC, RndSeed.R = Rnd.seed, RndSeed.WinBUGS = RndSeed.WinBUGS),
                  R2WB = fit)


   if(Plot){

    for(ii in 1:Ndoses3){

      flilenm <- paste(Outfile,"_", Doses3.Label[ii],".pdf", sep='')

      P1 <- outlist$P
      Pcat1 <- outlist$Pcat

      nms.P1 <- dimnames(P1)
      nms.Pcat1 <- dimnames(Pcat1)

      dim.P1 <- dim(P1)
      dim.Pcat1 <- dim(Pcat1)

      if(Ndoses1 == 1 || Ndoses2==1){

        P.1 <- array(P1[,,,ii], dim= dim.P1[1:3], dimnames = nms.P1[1:3])
        Pcat.1 <- array(Pcat1[,,,ii,"mean"], dim= dim.Pcat1[1:3], dimnames=nms.Pcat1[1:3])

      }

      if(Ndoses1 > 1 & Ndoses2 > 1){

        P.1 <- P1[,,,ii]
        Pcat.1 <- Pcat1[,,,ii,"mean"]

      }

      pdf(file = flilenm, onefile = F)


      Graph2Combo(
        P = P.1,
        Pcat = Pcat.1,
        Pcutoffs=Pcutoffs,
        int.crit= int.crit,
        int.col = int.col,
        drug1.labels= Doses1.1,
        drug2.labels= Doses2.1,
        lab.drug2 = Agent2,
        lab.drug1 = Agent1
      )

      dev.off()


    }

   }

   #Creating Output file with necessary

   sink(paste(Outfile, ".txt", sep = ""), append = F)
   cat("\n --------------------------------------------------------------------------------\n")
   cat("Dose Escalation: Bayesian Inference")
   cat("\n", date())
   cat("\n Working directory: ", getwd())
   cat("\n --------------------------------------------------------------------------------\n")
   cat("Data \n\n")
   cat("\n --------------------------------------------------------------------------------\n")
   cat("\n")
   print(outlist$Data)
   cat("\n")

   cat("\n --------------------------------------------------------------------------------\n")
   cat(" Posterior Results \n\n")
   cat("\n --------------------------------------------------------------------------------\n")
   cat("\n")
   cat("Parameters 1 \n\n")
   cat("\n")
   print(outlist$logAlphaBeta1)
   cat("\n")
   print(outlist$logAlphaEBeta1)
   cat("\n")
   cat("Parameters 2 \n\n")
   cat("\n")
   print(outlist$logAlphaBeta2)
   cat("\n")
   print(outlist$logAlphaEBeta2)
   cat("\n")
   cat("Parameters 3 \n\n")
   cat("\n")
   print(outlist$logAlphaBeta3)
   cat("\n")
   print(outlist$logAlphaEBeta3)
   cat("\n")

   cat("\n")
   cat("Interaction parameter (eta12) \n\n")
   cat("\n")
   print(outlist$eta12)
   cat("\n")
   print(outlist$OddsFactor12)
   cat("\n")
   cat("Interaction parameter (eta13) \n\n")
   cat("\n")
   print(outlist$eta13)
   cat("\n")
   print(outlist$OddsFactor13)
   cat("\n")
   cat("Interaction parameter (eta23) \n\n")
   cat("\n")
   print(outlist$eta23)
   cat("\n")
   print(outlist$OddsFactor23)
   cat("\n")

   cat("Interaction parameter (eta123) \n\n")
   cat("\n")
   print(outlist$eta123)
   cat("\n")
   print(outlist$OddsFactor123)
   cat("\n")
   cat("\n")
   cat("Interval probabilities and DLT rates\n\n")
   cat("\n")
   print(outlist$PcatP)
   cat("\n")

   cat("\n")
   cat("Gelman-Rubin diagnostics \n\n")
   cat("\n")
   print(outlist$Rhat)
   cat("\n")

   cat("\n")
   cat("Mixture Probabilities \n\n")
   cat("\n")
   cat("Agent 1")
   cat("\n")
   print(outlist$Post.pmix1)
   cat("\n")
   cat("Agent 2")
   cat("\n")
   print(outlist$Post.pmix2)
   cat("\n")
   cat("Agent 3")
   cat("\n")
   print(outlist$Post.pmix3)
   cat("\n")

   cat("\n")
   cat("\n --------------------------------------------------------------------------------\n")
   cat(" Model Specification \n\n")
   cat("\n --------------------------------------------------------------------------------\n")
   cat("\n")
   cat("Prior distributions \n\n")
   cat("\n")
   print(outlist$Priors)
   cat("\n")

   cat("\n")
   cat("MCMC parameters \n\n")
   cat("\n")
   print(outlist$MCMC)
   cat("\n")
   sink()

  return(outlist)

}
