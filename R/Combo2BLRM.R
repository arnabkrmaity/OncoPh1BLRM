#' @title Combo2BLRM Function
#'
#' @description This function performs analysis with dual combination model
#'
#' @param Data Dose DLT data -- a list objest contains DoseAdm1, DoseAdm2, Npat, Ntox, and Weight.
#' DoseAdm1 is a vector of first agent doses which are administered at patient cohorts. 
#' DoseAdm2 is a vector of second agent doses which are administered at patient cohorts. 
#' Npat is a vector of same length as DoseAdm, the elements are the cohort sizes at each DoseAdm. 
#' Ntox is a vector of same length as DoseAdm, the elements are the the number of toxicitites 
#' observed for each cohort. Weight is a vector of same length as DoseAdm, the elements are 1 
#' by default, but can be changed depending on how much information one wants borrow form this 
#' data set; 1 is maximum i.e. borrowing the full information, 0 is minimum i.e. borrowing 
#' no information, similary, setting 0.5 means borrowing half of the information
#' @param Agent1 Compound 1 name
#' @param Doses1 Doses for compound 1
#' @param DoseRef1 Dose reference for compound 1
#' @param Prior.logAB.1 Prior for single agent effect of compound 1 (a list object)
#' @param pmix.logAB1 mixing probability of the prior distributions for compound 1. 
#' Usually is set to c(1) however this is a vector of mixing probabilities when 
#' mixture of priors is used. See the package vignettes for details
#' @param Agent2 Compound 2 name
#' @param Doses2 Doses for compound 2
#' @param DoseRef2 Dose reference for compound 2
#' @param Prior.logAB.2 Prior for single agent effect of compound 2 (a list object)
#' @param pmix.logAB2 mixing probability of the prior distributions for compound 2. 
#' Usually is set to c(1) however this is a vector of mixing probabilities when 
#' mixture of priors is used. See the package vignettes for details
#' @param Prior.eta prior for interaction parameter \eqn{\eta}, a vector of mean and standard deviation
#' @param pmix.eta mixing probability of the prior distributions of the interaction parameter. 
#' Usually is set to c(1) however this is a vector of mixing probabilities when 
#' mixture of priors is used.
#' @param Prior.Only either TRUE or FALSE. If TRUE, Data must be equal to NULL meaning that there is 
#' no data and posterior summaries are basically the summaries from the prior distribution. If FALSE,
#' then Data must be provided and posterior will be computed.
#' @param New.cohortsize new for prediction
#' @param Outfile output file
#' @param Plot Indicator if Plot will be generated
#' @param int.crit Critical value for interval to change color
#' @param int.col Color of each bar
#' @param pars MCMC parameters to save
#' @param MCMC MCMC specifications -- nburn: burnin length, ntot: Markov chain length after burnin, 
#' nchain: number of Markov chains, nthin: thinning interval, default is 1 (no thinning)
#' @param Warn Prints convergence warning
#' @param DIC DIC to save
#' @param progress.bar Jags parameter for type of progress bar. Possible values are
#' "text", "gui", and "none". Type "text" is displayed on the R console. Type "gui"
#' is a graphical progress bar in a new window. The progress bar is suppressed if
#' progress.bar is "none"
#' @param Rnd.seed Seed for R
#' @param RSeed.WinBUGS Seed for WinBugs
#'
#' @details See vignettes. browseVignettes("OncoPh1BLRM") 
#' 
#' @references Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. 
#' A Bayesian Industry Approach to Phase I Combination Trials in Oncology. 
#' Statistical Methods for Drug Combination Studies. Taylor & Francis, 2015.
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
#' \item{eta}{Posterior summary of intercation paramter \eqn{\eta})}
#' \item{OddsFactor}{Posterior summary of the odds factor}
#' \item{P}{Predicted summaries of the dose levels}
#' \item{Pcat}{Posterior predictive probabilities of doses}
#' \item{PcatP}{Posterior predictive probabilities and summaries of doses}
#' \item{pred.Tox}{Predictive toxicity probilities of doses}
#' \item{Post.pmix1}{Posterior mixture weigths for agent 1}
#' \item{Post.pmix2}{Posterior mixture weigths for agent 2}
#' \item{Post.pmix2}{Posterior mixture weigths for \eqn{\eta}}
#' \item{Rhat}{The MCMC convergence dignostoc Rhat for all the parameters}
#' \item{MCMC}{MCMC informations}
#' \item{R2WB}{}
#' In addition, a text output containing all of the above and two plots -- posterior DLT rates and 
#' interval probabilities of DLT's, are produced. 
#' 
#' 
#' @examples 
#' \dontrun{
#' # Doses for agent 1
#' Doses1 <- c(0, 3, 4.5, 6)
#' DoseRef1 <- 3
#' 
#' #Doses for agent 2
#' Doses2 <- c(0, 400, 600, 800)
#' DoseRef2 <- 800
#' 
#' # Priors derived using MAP analysis
#' MAP.prior.Agent1 <- c(-2.53, 0.03, 1.02, 0.86, -0.54)
#' MAP.prior.Agent2 <- c(-2.12, 0.69, 0.61, 0.64, -0.20)
#' 
#' # Interaction prior
#' prior.interaction <- function(Odds.median= NULL, Odds.upper= NULL){
#' mean.eta <- log(Odds.median)
#' sd.eta <- (log(Odds.upper) - log(Odds.median))/1.96
#' return(list(mean.eta=mean.eta, sd.eta=sd.eta))
#' }
#' 
#' # Prior for eta: Median: No increase in Odds of DLT Upper: 9 fold increase in Odds of DLT
#' peta <- prior.interaction(Odds.median = 1.1, Odds.upper = 9)
#' etaPrior <- c(peta$mean.eta, peta$sd.eta)
#' 
#' # Duel agent data
#' data.combo2 <- list(DosesAdm1 = c(   3,   3,   6,   3,   6),
#'                     DosesAdm2 = c( 400, 800, 400, 800, 400),
#'                     Npat      = c(   3,   3,   3,   3,   3),
#'                     Ntox      = c(   0,   1,   1,   1,   0),
#'                     Weight    = c(   1,   1,   1,   1,   1)
#' )
#' 
#' # Analysis using dual agent BLRM
#' fit.BLRM.dual <-  Combo2BLRM(Data           = data.combo2 ,
#'                              Agent1         = "DRUG 1",
#'                              Doses1         = Doses1,
#'                              DoseRef1       = DoseRef1,
#'                              Prior.logAB.1  = list(MAP.prior.Agent1),
#'                              pmix.logAB1    = c(1),
#'                              Agent2         = "DRUG 2",
#'                              Doses2         = Doses2,
#'                              DoseRef2       = 800,
#'                              Prior.logAB.2  = list(MAP.prior.Agent2),
#'                              pmix.logAB2    = c(1),
#'                              Prior.eta      = list(etaPrior),
#'                              Prior.Only     = FALSE,
#'                              Pcutoffs       = c(0.16, 0.33),
#'                              New.cohortsize = 3,
#'                              Outfile        = "Output_dual",
#'                              MCMC           = c(2500, 15000, 4, 1),
#'                              DIC            = TRUE,
#'                              Rnd.seed       = 1452671,
#'                              RSeed.WinBUGS  = 1367
#' )
#' }
#' 
#' 
#' 
#' @export




Combo2BLRM <- function(Data               = NULL,
                       Agent1             = NULL,
                       Doses1             = NULL,
                       DoseRef1           = NULL,
                       Prior.logAB.1       = NULL,
                       pmix.logAB1        = c(1),
                       Agent2             = NULL,
                       Doses2             = NULL,
                       DoseRef2           = NULL,
                       Prior.logAB.2      = NULL,
                       pmix.logAB2        = c(1),
                       Prior.eta          = NULL,
                       pmix.eta           = c(1),
                       Prior.Only         = FALSE,
                       Pcutoffs           = c(0.16, 0.33),
                       New.cohortsize     = 3,
                       Outfile            = "Out3",
                       Plot               = TRUE,
                       int.crit           = c(1,1,0.25),
                       int.col            = list( "lightgreen", "green", c("yellow","red")),
                       pars               = c("logAlphaBeta1", "logAlphaBeta2", "eta","P12","pCat","pred.Tox"),
                       MCMC               = c(2500, 7500, 4, 1),
                       Warn               = TRUE,
                       DIC                = TRUE,
                       progress.bar       = "text",
                       Rnd.seed           = 1452671,
                       RSeed.WinBUGS      = 1367
                      )
{

  version = "30-November-2017"
  
  stopifnot("Mixture weights must add to 1" = sum(pmix.logAB1) == 1)
  stopifnot("Mixture weights must add to 1" = sum(pmix.logAB2) == 1)

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
  Ntox       = Data$Ntox*Data$Weight
  Npat       = Data$Npat*Data$Weight

  DosesAdm1[DosesAdm1 ==0] = 0.00001
  DosesAdm2[DosesAdm2 ==0] = 0.00001

  Ncohorts <- length(DosesAdm1)

  DosesAdm1 <- cbind(NULL, DosesAdm1)
  DosesAdm2 <- cbind(NULL, DosesAdm2)

  Ntox <- cbind(NULL, Ntox)
  Npat <- cbind(NULL, Npat)
  zero <- rep(0,Ncohorts)

                }

if(Prior.Only){Npat=0; Ntox=0}

  Ndoses1 <- length(Doses1)
  Ndoses2 <- length(Doses2)


  Ncat <- length(Pcutoffs) + 1
  Pcutoffs1 = c(0, Pcutoffs, 1)

  intLabels = paste(Pcutoffs1[1:Ncat], "-", Pcutoffs1[2:(Ncat + 1)], sep = "")
  Doses1.Label = paste(Agent1, "=", Doses1, sep = '')
  Doses2.Label = paste(Agent2, "=", Doses2, sep = '')
  pred.tox.Label = paste("DLT=",0:New.cohortsize, sep='')

  Doses1.1 =Doses1
  Doses2.1 = Doses2

  Doses1[Doses1 ==0] = 0.00001
  Doses2[Doses2 ==0] = 0.00001

  Doses1 <- cbind(NULL, Doses1)
  Doses2 <- cbind(NULL, Doses2)

  if(min(pmix.logAB1) < 1){pars <- c(pars, "Z1cat")}
  if(min(pmix.logAB2) < 1){pars <- c(pars, "Z2cat")}
  if(min(pmix.eta) < 1){pars <- c(pars, "Z.etacat")}

  #Priors

  #Bivariate normal prior for (logAlpha1,logBeta1) [marginal model for Agent1]

  Prior.logAB1 <- c(Prior.logAB.1, list(c(0,0,0.1,0.1,0)))
  Pmix1 <- c(pmix.logAB1, 0)
  Prior1 <- do.call(rbind,Prior.logAB1)

  Pmix1 <- cbind(NULL, Pmix1)
  Nmix1 <- length(Prior.logAB1)


  #Bivariate normal prior for (logAlpha2,logBeta2) [marginal model for Agent2]

  Prior.logAB2 <- c(Prior.logAB.2, list(c(0,0,0.1,0.1,0)))
  Pmix2 <- c(pmix.logAB2, 0)
  Prior2 <- do.call(rbind,Prior.logAB2)

  Pmix2 <- cbind(NULL,Pmix2)
  Nmix2 <- length(Prior.logAB2)


  #Eta priors

  Prior.eta <- c( Prior.eta, list(c(0,1)))
  P.eta.mix  = c(pmix.eta,0)

  etaPrior <- do.call(rbind,Prior.eta)
  P.eta.mix <- cbind(NULL,P.eta.mix)
  Nmix.eta <- length(Prior.eta)



  #Choice of model

  if(Prior.Only){model= Combo2BLRM.Prior.WB; model.index=1}
  if(!Prior.Only & (any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)))){model= Combo2BLRM.Weight.WB; model.index=2}
  if(!Prior.Only & !(any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)))){model= Combo2BLRM.WB; model.index=3}

  #Data for jags
  if(Prior.Only){

    data= list("Ndoses1", "Ndoses2", "Ncat", "Nmix1", "Nmix2","Nmix.eta",
               "Prior1", "Prior2","Pmix1", "Pmix2","etaPrior","P.eta.mix",
               "Doses1", "Doses2", "DoseRef1", "DoseRef2",
               "Pcutoffs", "New.cohortsize")
  }
  if(!Prior.Only & model.index== 3){

    data= list("Ncohorts", "Ndoses1", "Ndoses2", "Ncat", "Nmix1", "Nmix2", "Nmix.eta",
               "Prior1", "Prior2", "Pmix1", "Pmix2","etaPrior", "P.eta.mix",
               "Doses1", "Doses2", "DoseRef1", "DoseRef2",
               "DosesAdm1", "DosesAdm2", "Ntox", "Npat",
               "Pcutoffs", "New.cohortsize")
  }

  if(!Prior.Only & model.index== 2){

    data= list("Ncohorts", "Ndoses1", "Ndoses2", "Ncat", "Nmix1", "Nmix2", "Nmix.eta",
               "Prior1", "Prior2", "Pmix1", "Pmix2","etaPrior", "P.eta.mix",
               "Doses1", "Doses2", "DoseRef1", "DoseRef2",
               "DosesAdm1", "DosesAdm2", "Ntox", "Npat", "zero",
               "Pcutoffs", "New.cohortsize")
  }


  #Initial values
  inits.fun = function(i){

    logAlphaBeta1All <- matrix(NA, nrow=Nmix1, ncol=2)
    logAlphaBeta2All <- matrix(NA, nrow=Nmix2, ncol=2)

    for(i1 in 1:Nmix1){logAlphaBeta1All[i1,1:2] <- rmvnorm(1, sigma = diag(2))}
    for(i2 in 1:Nmix2){logAlphaBeta2All[i2,1:2] <- rmvnorm(1, sigma = diag(2))}

    eta12 <- rnorm(Nmix.eta, 0, 0.1)

    Z1 <- 1
    Z2 <- 1

    return(list(logAlphaBeta1All = logAlphaBeta1All,
                logAlphaBeta2All = logAlphaBeta2All,
                eta12 = eta12,
                Z1=Z1,
                Z2=Z2
                )
           )
  }

  inits <- lapply(rep(1,n.chains),inits.fun)

  RndSeed.WinBUGS <- RSeed.WinBUGS


  pars = c(pars,"logAlphaEBeta1", "logAlphaEBeta2","logAlphaBeta1.XY", "logAlphaBeta2.XY", "OddsFactor")


  fit = jags(
    data=data,
    inits= inits,
    parameters.to.save = pars,
    model=model,
    n.chains=n.chains,n.burnin=n.burnin,n.iter=n.iter,n.thin=n.thin,
    jags.seed=RndSeed.WinBUGS,
    DIC= DIC,
    progress.bar = progress.bar
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


  eta.neff <- fit$BUGSoutput$summary[c("eta"), "n.eff"]
  if (any(eta.neff  < 1000)) {
    cat("\nWARNING: Neff < 1000 for eta  Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for eta Consider increasing the sample size (number of iterations).")
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


  eta = NULL
  Rhat.eta = NULL
  OddsFactor =NULL
  Rhat.OddsFactor=NULL

  P = NULL
  Rhat.P = NULL

  Pcat = NULL
  PcatP = NULL
  pred.Tox =NULL

  Post.pmix1 =NULL
  Post.pmix2 =NULL
  Post.pmix.eta = NULL


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



  if(is.element("eta", pars)){
    eta = BUGS.Select.Output("eta", summary, cols = vnames)
    Rhat.eta = eta["Rhat"]
    eta <- rbind(NULL, eta)

    rownames(eta) = "eta"

  }

  if(is.element("OddsFactor", pars)){
    OddsFactor = BUGS.Select.Output("OddsFactor", summary, cols = vnames)
    Rhat.OddsFactor = OddsFactor["Rhat"]
    OddsFactor <- rbind(NULL, OddsFactor)

    rownames(OddsFactor) = "OddsFactor"

  }

  if (is.element("P12", pars)) {
    P = BUGS.Select.Output("P12", summary, cols = vnames)
    Rhat.P = P[, "Rhat"]
    aiw = round(Beta.ms2ab(P[, "mean"], P[, "sd"])$n, 1)
    P = cbind(P, aiw)
    P = BUGS.Table2Array(P, Labels = list(Doses1.Label, Doses2.Label), Dim = c(Ndoses1, Ndoses2))
    P = aperm(P, c(1,3,2))

  }


  if (is.element("pCat", pars)) {
    Pcat = BUGS.Select.Output("pCat", summary, cols = vnames)
    Pcat = BUGS.Table2Array(Pcat, Labels = list(Doses1.Label, Doses2.Label, intLabels), Dim = c(Ndoses1, Ndoses2,  Ncat))
    Pcat = Pcat[, , , "mean", drop = FALSE]
    Pcat = aperm(Pcat, c(1, 3, 2, 4))

  }

 if (is.element("pCat", pars) & is.element("P12", pars)) {

     PcatP <- array(NA, dim=c(Ndoses1, 10, Ndoses2, 1))

     for(k1 in 1:Ndoses2){

         PcatP[ ,,k1,1] <- cbind(rbind(NULL,Pcat[,,k1,1]), rbind(NULL,P[,,k1]))

                          }


     vnames1 = c(intLabels,vnames,"aiw")

     dimnames(PcatP) = list(Doses1.Label, vnames1, Doses2.Label, "")
 }


  if (is.element("pred.Tox", pars)){

    pred.Tox = BUGS.Select.Output("pred.Tox", summary, cols = vnames)
    Rhat.pred.Tox= pred.Tox[, "Rhat"]



    pred.Tox = BUGS.Table2Array(pred.Tox, Labels = list(Doses1.Label, Doses2.Label,pred.tox.Label), Dim = c(Ndoses1, Ndoses2, New.cohortsize+1))
    pred.Tox = pred.Tox[, , , "mean", drop = FALSE]
    pred.Tox = aperm(pred.Tox, c(1,3, 2, 4))

  }

  if (is.element("Z1cat", pars)){

    Post.pmix1.1 = BUGS.Select.Output("Z1cat", summary, cols = "mean")
    Post.pmix1 <- cbind(NULL,Post.pmix1.1[-Nmix1])
    Post.pmix1.nms <- paste("mix1",1:(Nmix1-1),sep='')
    rownames(Post.pmix1) <- Post.pmix1.nms

  }


  if (is.element("Z2cat", pars)){

    Post.pmix2.1 = BUGS.Select.Output("Z2cat", summary, cols = "mean")
    Post.pmix2 <- cbind(NULL,Post.pmix2.1[-Nmix2])
    Post.pmix2.nms <- paste("mix2",1:(Nmix2-1),sep='')
    rownames(Post.pmix2) <- Post.pmix2.nms

  }


  if (is.element("Z.etacat", pars)){

    Post.pmix.eta.1 = BUGS.Select.Output("Z.etacat", summary, cols = "mean")
    Post.pmix.eta <- cbind(NULL,Post.pmix.eta.1[-Nmix.eta])
    Post.pmix.eta.nms <- paste("mix.eta",1:(Nmix.eta-1),sep='')
    rownames(Post.pmix.eta) <- Post.pmix.eta.nms

  }


   Rhat = c(Rhat.logAlphaBeta1, Rhat.logAlphaBeta1.cor, Rhat.logAlphaEBeta1,
            Rhat.logAlphaBeta2, Rhat.logAlphaBeta2.cor, Rhat.logAlphaEBeta2,
            Rhat.eta,
            Rhat.P)

   Rhat.max = max(Rhat)

 if (Rhat.max > 1.1) {
    warning("\n Warning: at least one of the parameters has a Gelman-Rubin diagnostic Rhat > 1.1")
 }

   Priors = list(Prior.logAB1 = Prior.logAB.1,
                 pmix.logAB1 = pmix.logAB1,
                 Prior.logAB2 = Prior.logAB.2,
                 pmix.logAB2 = pmix.logAB2,
                 etaPrior = etaPrior
)

   labels = list(intLabels = intLabels)
   Agents = c(Agent1, Agent2)

   if(!Prior.Only){
   Data.1 = data.frame(t(do.call(rbind, Data)))
   colnames(Data.1) = c("DoseAdm1", "DoseAdm2","Npat","Ntox", "Weight")
   N = c(Ncohorts = Ncohorts, Ndoses1 = Ndoses1, Ndoses2=Ndoses2, Nmix1=Nmix1-1, Nmix2=Nmix2-1)

   Data11 <- aggregate(Data.1, by=list(Data.1$DoseAdm1, Data.1$DoseAdm2, Data.1$Weight),sum)[,c(1,2,3,6,7)]
   Data1 <- Data11[,c(1,2,4,5,3)]
   colnames(Data1) = c("DoseAdm1", "DoseAdm2","Npat","Ntox", "Weight")

   }

   if(Prior.Only){
     Data1 = NULL
     N = c(Ndoses1 = Ndoses1, Ndoses2=Ndoses2,Nmix1=Nmix1-1, Nmix2=Nmix2-1, Nmix.eta=Nmix.eta-1)
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
                  eta = eta,
                  OddsFactor = OddsFactor,
                  P = P,
                  Pcat = Pcat,
                  PcatP = PcatP,
                  pred.Tox = pred.Tox,
                  Post.pmix1 = Post.pmix1,
                  Post.pmix2 = Post.pmix2,
                  Post.pmix.eta = Post.pmix.eta,
                  Rhat = list(Rhat = Rhat, Rhat.max = Rhat.max),
                  MCMC = list(pars.monitored = pars, MCMC = MCMC, RndSeed.R = Rnd.seed, RndSeed.WinBUGS = RndSeed.WinBUGS),
                  R2WB = fit)



   if(Plot){

     P1 <- outlist$P
     Pcat1 <- outlist$Pcat

     nms.P1 <- dimnames(P1)
     nms.Pcat1 <- dimnames(Pcat1)

     dim.P1 <- dim(P1)
     dim.Pcat1 <- dim(Pcat1)

     if(Ndoses1 == 1 || Ndoses2==1){

       P.1 <- array(P1, dim= dim.P1[1:3], dimnames = nms.P1[1:3])
       Pcat.1 <- array(Pcat1[,,,"mean"], dim= dim.Pcat1[1:3], dimnames=nms.Pcat1[1:3])

     }

     if(Ndoses1 > 1 & Ndoses2 > 1){

       P.1 <- P1
       Pcat.1 <- Pcat1[,,,"mean"]

     }

     pdf(file = paste(Outfile,'.pdf',sep=''), onefile = F)

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

   cat("\n")
   cat("Interaction parameter (eta) \n\n")
   cat("\n")
   print(outlist$eta)
   cat("\n")
   print(outlist$OddsFactor)
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
   cat("Interaction")
   cat("\n")
   print(outlist$Post.pmix.eta)
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
