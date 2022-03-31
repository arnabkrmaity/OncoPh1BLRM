#' @title SingleBLRM.PK Function
#'
#' @description This function performs analysis with Bayesian logistic regression (BLRM) with PK data for 
#'              single agent dose escalation
#'
#' @param Data  Dose DLT data (a list objest contains DoseAdm,  Npat, Ntox, Weight and X (covariate for exposure))
#' @param Agent Compound name
#' @param Doses Doses for compound
#' @param Doses.Label Label for dose groups
#' @param DoseRef Dose reference for compound
#' @param Prior.logAB Prior for Bayesian Logistic Regession Model for Exposure-DLT (a list object)
#' @param pmix.logAB mix probability  Bayesian Logistic Regession Model for Exposure-DLT
#' @param monotone.logAB Monotonicity of Exposure-DLT model(values TRUE (default)/FALSE)
#' @param Prior.logAB.PK Prior for Dose-exposure model (a list object)
#' @param pmix.logAB.PK mix probability for Dose-exposure model
#' @param monotone.PK Monotonicity of Dose-exposure model(values TRUE (default)/FALSE)
#' @param Prior.sd.PK Prior for standard deviation for PK parameter
#' @param random.model logical: if TRUE, a random effect model will be used for dose-exposure analysis
#' @param random.slope logical:if TRUE, model with random intercept and slope will be used for dose-exposure model 
#' @param Prior.re.sd.PK Prior for subject random effect standard deviation
#' @param Prior.re.rho prior for correlation between slope and interect of random effect
#' @param Xout Output pattern for covariates
#' @param Xout.label Output pattern labels
#' @param Prior.beta Prior for regression coefficient for covariates in Dose-exposure model
#' @param Prior.Only Prior only
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
#' type of progress bar. Possible values are "text", "gui", and "none". Type "text" is
#' displayed on the R console. Type "gui" is a graphical progress bar in a new window.
#' The progress bar is suppressed if progress.bar is "none"
#' @param Rnd.seed Seed for R
#' @param RSeed.WinBUGS Seed for WinBugs
#'
#'
#' @details To be added
#' @references
#' 
#' 
#' @return 
#' \item{Data}{Two list objects of doses and PK (exposure) inputed respectively}
#' \item{Priors}{Input priors}
#' \item{logAlphaBeta}{Posterior summary of (\eqn{\log(\alpha)}, \eqn{\log(\beta)})}
#' \item{logAlphaBeta.cor}{Posterior summary of the correlation between (\eqn{\log(\alpha)}, \eqn{\log(\beta)})}
#' \item{logAlphaEBeta}{Posterior summary of (\eqn{\log(\alpha)}, \eqn{\beta})}
#' \item{logGamma}{Posterior summary of (\eqn{\log(\gamma_0)}, \eqn{\log(\gamma_1)})}
#' \item{logGamma.cor}{Posterior summary of the correlation between (\eqn{\log(\gamma_0)}, \eqn{\log(\gamma_1)})}
#' \item{logEGamma}{Posterior summary of (\eqn{\log(\gamma_0)}, \eqn{\gamma_1})}
#' \item{sd.PK}{Posterior summary of subject random effect standard deviation}
#' \item{sd.re.PK}{Posterior summary of subject random effect standard deviation}
#' \item{beta}{Posterior summary of regression coefficient for covariates in Dose-exposure model}
#' \item{P}{Posterior summary of DLT probabilities}
#' \item{Pcat}{Posterior summary of DLT interval probabilities}
#' \item{PcatP}{P and Pcat combined}
#' \item{pred.Tox}{Predicted probabilities for the number of DLT's for a cohort of size New.cohortSize}
#' \item{Post.pmix}{Posterior mixture weigths (for mixture prior) for exposure-DLT model}
#' \item{Post.pmix.PK}{Posterior mixture weigths (for mixture prior) for Dose-exposure model}
#' \item{Rhat}{Gelman-Rubin MCMC convergence diagnostics}
#' \item{MCMC}{MCMC specifications}
#' \item{R2WB}{WinBUGS object with summaries of MCMC analysis (see \link{R2WinBUGS})}
#' 
#' 
#' @seealso \code{\link{SingleBLRM}}
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#'
#' Data.Dose.PK <- list(DosesAdm = c(0.1, 0.1, 0.3, 0.3, 0.3, 1.0, 1.0, 3.0, 3.0, 3.0,
#'                                   10.0, 10.0, 30.0, 30.0, 30.0, 30.0, 30.0, 50.0, 50.0, 50.0),
#'                      Dose.Grp = c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 6, 6,
#'                                   7, 7, 7),
#'                      Npat     = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#'                      Ntox     = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1),
#'                      Weight   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#'                      PK       = c(2.82, 1.68, 2.87, 8.14, 4.39, 13.30, 12.20, 79.70, 47.20,
#'                                   41.00, 297.00, 164.00, 401.00, 766.00, 592.00, 769.00,
#'                                   577.00, 834.00, 1160.00, 1160.00)/median(c(2.82, 1.68, 2.87, 
#'                                   8.14, 4.39, 13.30, 12.20, 79.70, 47.20, 41.00, 297.00, 164.00, 
#'                                   401.00, 766.00, 592.00, 769.00, 577.00, 834.00, 1160.00, 
#'                                   1160.00)),
#'                      X        = NULL
#'                      )
#'                      
#' Doses <- c(0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 50.0)
#'
#'
#' # Analysis: Fixed effect model assuming monotone relationship between PK and Doses
#' 
#' fit.BLRM.PK.Fixed.mono <- SingleBLRM.PK(Data           = Data.Dose.PK,
#'                                         Agent          = "DRUG 1",
#'                                         Doses          = Doses,
#'                                         Doses.Label    = as.character(Doses),
#'                                         DoseRef        = 270,
#'                                         Prior.logAB    = list(c(log(0.50), 0, 2, 1, 0)),
#'                                         pmix.logAB     = c(1),
#'                                         monotone.logAB = TRUE,
#'                                         Prior.logAB.PK = list(c(1, 0, 3, 1, 0)),
#'                                         pmix.logAB.PK  = c(1),
#'                                         monotone.PK    = TRUE,
#'                                         Prior.sd.PK    = c(0.25, log(2)/1.96),
#'                                         random.model   = FALSE,
#'                                         Prior.beta     = rbind(NULL,c(0, 10)),
#'                                         Prior.Only     = FALSE,
#'                                         Pcutoffs       = c(0.16, 0.33),
#'                                         New.cohortsize = 3,
#'                                         Outfile        = "out_BLRMPK_Fixed_mono",
#'                                         MCMC           = c(5000, 10000, 4, 1),
#'                                         DIC            = TRUE,
#'                                         Rnd.seed       = 1452671,
#'                                         RSeed.WinBUGS  = 1367
#'                                         )
#' 
#' }
#' 
#' @export


SingleBLRM.PK <- function(Data               = NULL,
                          Agent              = NULL,
                          Doses              = NULL,
                          Doses.Label        = NULL,
                          DoseRef            = NULL,
                          Prior.logAB        = NULL,
                          pmix.logAB         = c(1),
                          monotone.logAB     = TRUE,
                          Prior.logAB.PK     = NULL,
                          pmix.logAB.PK      = c(1),
                          monotone.PK        = TRUE,
                          Prior.sd.PK        = NULL,
                          random.model       = TRUE,
                          random.slope       = FALSE,
                          Prior.re.sd.PK     = NULL,
                          Prior.re.rho       = c(0,1),
                          Xout               = NULL,
                          Xout.label         = NULL,
                          Prior.beta         = c(0,10),
                          Prior.Only         = FALSE,
                          Pcutoffs           = c(0.16, 0.33),
                          New.cohortsize     = 3,
                          Outfile            = "Out",
                          Plot               = TRUE,
                          int.crit           = c(1,1,0.25),
                          pars               = c("logAlphaBeta", "logGamma","pTox","pCat","pred.Tox","Zcat","Zcat.PK", "sd.PK"),
                          MCMC               = c(2500, 7500, 4, 1),
                          Warn               = TRUE,
                          DIC                = TRUE,
                          progress.bar       = "text",
                          Rnd.seed           = 1452671,
                          RSeed.WinBUGS      = 1367
)
{


  version = "15-September-2019"

  set.seed(Rnd.seed)
  
  stopifnot("Mixture weights must add to 1" = sum(pmix.logAB) == 1)
  stopifnot("Mixture weights must add to 1" = sum(pmix.logAB.PK) == 1)

  if(Prior.Only){DIC=FALSE}

  # MCMC parameters
  n.burnin = MCMC[1]
  n.iter = MCMC[2]
  n.chains = MCMC[3]
  n.thin = MCMC[4]

  #Reading Data

  if(!Prior.Only){

    DosesAdm   = Data$DosesAdm
    Dose.Grp   = Data$Dose.Grp
    Ntox       = Data$Ntox*Data$Weight
    Npat       = Data$Npat*Data$Weight
    PK         = Data$PK

    DosesAdm[DosesAdm == 0] = 0.00001

    Ncohorts <- length(DosesAdm)

    DosesAdm<- cbind(NULL, DosesAdm)

    Ntox <- cbind(NULL, Ntox)
    Npat <- cbind(NULL, Npat)
    PK <- cbind(NULL, PK)
    zero <- rep(0,Ncohorts)

  }

  if(Prior.Only){Npat=0; Ntox=0; Ncohorts=0}

  X = Data$X


  if (is.null(X)) {
    # create X matrix with value = 0
    withCov = FALSE
    X = matrix(0,Ncohorts,1)
    Xout=matrix(0,1,1)
    }

  if (is.null(Xout)) {
    # create Xout matrix with value = 0
    Xout=matrix(0,1,1)

  }

  if(is.null(Xout.label)){Xout.label=""}else{Xout.label = Xout.label}

  if (!is.matrix(X))
    X = matrix(X,1)

  if (!is.matrix(Xout))
    Xout = matrix(Xout,1)

  Nout = nrow(Xout)
  Ncov = ncol(X)

  Ndoses <- length(Doses)

  Ncat <- length(Pcutoffs) + 1
  Pcutoffs1 = c(0, Pcutoffs, 1)

  intLabels = paste(Pcutoffs1[1:Ncat], "-", Pcutoffs1[2:(Ncat + 1)], sep = "")

  if(is.null(Doses.Label)){Doses.Label = paste(Agent, "=", Doses, sep = '')}else{Doses.Label =Doses.Label}
  pred.tox.Label = paste("DLT=",0:New.cohortsize, sep='')

  Doses.1 =Doses

  Doses[Doses ==0] = 0.00001

  Doses <- cbind(NULL, Doses)

  if(min(pmix.logAB) < 1){pars <- c(pars, "Zcat")}
  if(min(pmix.logAB.PK) < 1){pars <- c(pars, "Zcat.PK")}

  #Priors

  #Bivariate normal prior for (logAlpha,logBeta)

  Prior.logAB.1 <- c(Prior.logAB, list(c(0,0,0.1,0.1,0)))
  Pmix <- c(pmix.logAB, 0)
  Prior <- do.call(rbind,Prior.logAB.1)

  Pmix <- cbind(NULL, Pmix)
  Nmix1 <- length(Prior.logAB.1)
  monotone.DLT <- ifelse(monotone.logAB == TRUE, 1, 2)

  #Bivariate normal prior for (logGamma1,logGamma2)

  Prior.PK.1 <- c(Prior.logAB.PK, list(c(0,0,0.1,0.1,0)))
  Pmix.PK <- c(pmix.logAB.PK, 0)
  Prior.PK <- do.call(rbind,Prior.PK.1)

  Pmix.PK <- cbind(NULL, Pmix.PK)
  Nmix2 <- length(Prior.PK.1)
  monotone.DE <- ifelse(monotone.PK == TRUE, 1, 2)
  rand.m <- ifelse(random.model == TRUE, 1, 2)
  rand.slope <- ifelse(random.slope == TRUE, 1, 2)
  
  if(rand.m ==2){rand.slope <- 2}
  
  if(rand.m ==2){Prior.re.sd.PK <- rbind(c(0, 0.001), c(0, 0.001))}

  #Normal prior for beta

  Prior.beta <- Prior.beta

  #Prior for sd
  

  Prior.sd.PK <- rbind(NULL,Prior.sd.PK)
  
  if(nrow(Prior.sd.PK)==1){Prior.sd.PK <- matrix(rep(Prior.sd.PK, times=Ndoses), nrow=Ndoses, byrow=T)}

  
  #Prior for random effect sd

  Prior.re.sd.PK <- rbind(NULL, Prior.re.sd.PK)
  if(nrow(Prior.re.sd.PK)==1){Prior.re.sd.PK <- rbind(Prior.re.sd.PK, c(Prior.re.sd.PK[1,1]/2, Prior.re.sd.PK[1,2]))}

  Prior.re.rho <- Prior.re.rho

  #Choice of model

  if(Prior.Only){model= SingleBLRM.PK.Prior.WB; model.index=1}
  if(!Prior.Only & (any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)))){model= SingleBLRM.PK.Weight.WB; model.index=2}
  if(!Prior.Only & !(any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)))){model=SingleBLRM.PK.WB; model.index=3}

  #Data for jags
  if(Prior.Only){

    data= list("Ndoses", "Ncat", "Nmix1", "Nmix2",
               "Ncov", "monotone.DLT", "monotone.DE", "Nout",
               "Prior", "Prior.PK", "Prior.beta", "Prior.sd.PK","Prior.re.sd.PK","Pmix", "Pmix.PK", "Prior.re.rho",
               "Doses", "DoseRef","Xout", 
               "Pcutoffs", "New.cohortsize")
  }

  if(!Prior.Only & model.index== 3){

    data= list("Ncohorts", "Ndoses", "Ncat", "Nmix1", "Nmix2",
               "Ncov", "monotone.DLT", "monotone.DE", "Nout",
               "Prior", "Prior.PK", "Prior.beta", "Prior.sd.PK","Prior.re.sd.PK","Pmix", "Pmix.PK","rand.m", "rand.slope", "Prior.re.rho",
               "Doses", "DoseRef","DosesAdm", "Ntox", "Npat","PK","X", "Xout", "Dose.Grp",
               "Pcutoffs", "New.cohortsize")
  }


  if(!Prior.Only & model.index== 2){

    data= list("Ncohorts", "Ndoses", "Ncat", "Nmix1", "Nmix2",
               "Ncov", "monotone.DLT", "monotone.DE", "Nout",
               "Prior", "Prior.PK", "Prior.beta", "Prior.sd.PK", "Prior.re.sd.PK", "Pmix", "Pmix.PK","rand.m", "rand.slope", "Prior.re.rho",
               "Doses", "DoseRef","DosesAdm", "Ntox", "Npat","PK","X", "Xout", "Dose.Grp","zero",
               "Pcutoffs", "New.cohortsize")
  }


  #Initial values
  inits.fun = function(i){

    logAlphaBetaAll <- matrix(NA, nrow=Nmix1, ncol=2)

    for(i1 in 1:Nmix1){logAlphaBetaAll[i1,1:2] <- rmvnorm(1, sigma = diag(2))}

    Z <- 1

    logGammaAll <- matrix(NA, nrow=Nmix2, ncol=2)

    for(i2 in 1:Nmix2){logGammaAll[i2,1:2] <- rmvnorm(1, sigma = diag(2))}
    Z.PK <- 1

    return(list(logAlphaBetaAll = logAlphaBetaAll,
                Z=Z,
                logGammaAll = logGammaAll,
                Z.PK=Z.PK
    ))
  }

  inits <- lapply(rep(1,n.chains),inits.fun)

  RndSeed.WinBUGS <- RSeed.WinBUGS

  pars = c(pars,"logAlphaEBeta", "logAlphaBeta.XY", "logEGamma", "logGamma.XY")

  if(!Prior.Only){pars=c(pars,"sd.re.PK")}


  if(DIC){rjags::load.module('dic')}

  fit = jags(
    data=data,
    inits= inits,
    parameters.to.save = pars,
    model=model,
    n.chains=n.chains,n.burnin=n.burnin,n.iter=n.iter,n.thin=n.thin,
    jags.seed=RndSeed.WinBUGS,
    DIC= DIC,
    progress.bar=progress.bar
    #bugs.directory="C:/Users/roychs04/Documents/Folder/Software/WinBUGS14"
  )


  #Processing output

  vnames = c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat")
  summary = fit$BUGSoutput$summary
  #fit$BUGSoutput$sims.matrix = NULL
  fit$BUGSoutput$sims.array = NULL

  if(Warn){
    logAlphaBeta.neff <- fit$BUGSoutput$summary[c("logAlphaBeta[1]", "logAlphaBeta[2]"), "n.eff"]

    if (any(logAlphaBeta.neff  < 1000)) {
      cat("\nWARNING: Neff < 1000 for logAlphaBeta  Consider increasing the sample size (number of iterations).\n")
      message("WARNING: Neff < 1000 for logAlphaBeta\ Consider increasing the sample size (number of iterations).")
    }

    logGamma.neff <- fit$BUGSoutput$summary[c("logGamma[1]", "logGamma[2]"), "n.eff"]

    if (any(logGamma.neff  < 1000)) {
      cat("\nWARNING: Neff < 1000 for logGamma  Consider increasing the sample size (number of iterations).\n")
      message("WARNING: Neff < 1000 for logGamma\ Consider increasing the sample size (number of iterations).")
    }

  }


  logAlphaBeta = NULL
  Rhat.logAlphaBeta = NULL
  logAlphaBeta.cor=NULL
  Rhat.logAlphaBeta.cor=NULL
  logAlphaEBeta =NULL
  Rhat.logAlphaEBeta =NULL

  logGamma = NULL
  Rhat.logGamma = NULL
  logGamma.cor=NULL
  Rhat.logGamma.cor=NULL
  logEGamma =NULL
  Rhat.logEGamma =NULL

  beta = NULL
  Rhat.beta = NULL

  sd.PK =NULL
  Rhat.sd.PK = NULL

  sd.re.PK =NULL
  Rhat.sd.re.PK = NULL

  P = NULL
  Rhat.P = NULL

  Pcat = NULL
  PcatP = NULL
  pred.Tox =NULL
  Post.pmix = NULL
  Post.pmix.PK = NULL



  if(is.element("logAlphaBeta", pars)){
    logAlphaBeta = BUGS.Select.Output("logAlphaBeta", summary, cols = vnames)
    Rhat.logAlphaBeta = logAlphaBeta[, "Rhat"]

    logAlphaBeta = BUGS.Table2Array(logAlphaBeta, Labels = list(c("logAlpha", "logBeta")), Dim = 2)

    logAlphaEBeta = BUGS.Select.Output("logAlphaEBeta", summary, cols = vnames)
    Rhat.logAlphaEBeta = logAlphaEBeta[, "Rhat"]

    logAlphaEBeta = BUGS.Table2Array(logAlphaEBeta, Labels = list(c("logAlpha", "beta")), Dim = 2)

  }

  if(is.element("logAlphaBeta.XY", pars)){

    logAlphaBeta.XY = BUGS.Select.Output("logAlphaBeta.XY", summary, cols = vnames)
    Rhat.logAlphaBeta.XY.cor = logAlphaBeta.XY["Rhat"]

    logAlphaBeta.cor <-  (logAlphaBeta.XY[1] - logAlphaBeta[1, 1]*logAlphaBeta[2,1])/logAlphaBeta[1, 2]/logAlphaBeta[2, 2]

  }

  if(is.element("logGamma", pars)){
    logGamma = BUGS.Select.Output("logGamma", summary, cols = vnames)
    Rhat.logGamma = logGamma[, "Rhat"]

    logGamma = BUGS.Table2Array(logGamma, Labels = list(c("logGamma.1", "logGamma.2")), Dim = 2)

    logEGamma = BUGS.Select.Output("logEGamma", summary, cols = vnames)
    Rhat.logEGamma = logEGamma[, "Rhat"]

    logEGamma = BUGS.Table2Array(logEGamma, Labels = list(c("logGamma.1", "Gamma.2")), Dim = 2)

  }

  if(is.element("logGamma.XY", pars)){

    logGamma.XY = BUGS.Select.Output("logGamma.XY", summary, cols = vnames)
    Rhat.logGamma.XY.cor = logGamma.XY["Rhat"]

    logGamma.cor <-  (logGamma.XY[1] - logGamma[1, 1]*logGamma[2,1])/logGamma[1, 2]/logGamma[2, 2]

  }


  if(is.element("beta", pars)){

    beta = BUGS.Select.Output("beta", summary, cols = vnames)
    Rhat.beta = beta["Rhat"]
  }

  if(is.element("sd.PK", pars)){

    sd.PK = BUGS.Select.Output("sd.PK", summary, cols = vnames)
    Rhat.sd.PK = sd.PK["Rhat"]
  }
  

  if(is.element("sd.re.PK", pars)){

    sd.re.PK = BUGS.Select.Output("sd.re.PK", summary, cols = vnames)
    Rhat.sd.re.PK = sd.re.PK["Rhat"]
  }
  

  if (is.element("pTox", pars)) {
    P = BUGS.Select.Output("pTox", summary, cols = vnames)
    Rhat.P = P[, "Rhat"]
    aiw = round(Beta.ms2ab(P[, "mean"], P[, "sd"])$n, 1)
    P = cbind(P, aiw)
    P = BUGS.Table2Array(P, Labels = list(Doses.Label, Xout.label), Dim = c(Ndoses,Ncov))
    P = aperm(P, c(1,3,2))

  }
  


  if (is.element("pCat", pars)) {
    Pcat = BUGS.Select.Output("pCat", summary, cols = vnames)
    Pcat = BUGS.Table2Array(Pcat, Labels = list(Doses.Label, intLabels, Xout.label), Dim = c(Ndoses, Ncat, Nout))
    Pcat = Pcat[, , , "mean", drop = FALSE]
    #Pcat = aperm(Pcat, c(1, 3, 2, 4))

  }
  


  if (is.element("pCat", pars) & is.element("pTox", pars)){

    PcatP <- cbind(rbind(NULL,Pcat[,,,1]), rbind(NULL,P[,,1]))


    vnames1 = c(intLabels,vnames,"aiw")

    dimnames(PcatP) = list(Doses.Label, vnames1)
  }

  if (is.element("pred.Tox", pars)){

    pred.Tox = BUGS.Select.Output("pred.Tox", summary, cols = vnames)
    Rhat.pred.Tox= pred.Tox[, "Rhat"]



    pred.Tox = BUGS.Table2Array(pred.Tox, Labels = list(Doses.Label, pred.tox.Label, Xout.label), Dim = c(Ndoses, New.cohortsize+1, Nout))
    pred.Tox = pred.Tox[, ,,"mean", drop = FALSE]

  }

 
  if (is.element("Zcat", pars)){

    Post.pmix.1 = BUGS.Select.Output("Zcat", summary, cols = "mean")
    Post.pmix <- cbind(NULL,Post.pmix.1[-Nmix1])
    Post.pmix.nms <- paste("mix",1:(Nmix1-1),sep='')
    rownames(Post.pmix) <- Post.pmix.nms

  }


  if (is.element("Zcat.PK", pars)){

    Post.pmix.PK.1 = BUGS.Select.Output("Zcat.PK", summary, cols = "mean")
    Post.pmix.PK <- cbind(NULL,Post.pmix.PK.1[-Nmix2])
    Post.pmix.PK.nms <- paste("mix",1:(Nmix2-1),sep='')
    rownames(Post.pmix.PK) <- Post.pmix.PK.nms

  }


  Rhat = c(Rhat.logAlphaBeta, Rhat.logAlphaBeta.cor, Rhat.logGamma, Rhat.logGamma.cor,Rhat.P)

  Rhat.max = max(Rhat)

  if (Rhat.max > 1.1) {
    warning("\n Warning: at least one of the parameters has a Gelman-Rubin diagnostic Rhat > 1.1")
  }

  Priors = list(Prior.logAB = Prior.logAB,
                pmix.logAB = pmix.logAB,
                Prior.logGamma = Prior.logAB.PK,
                pmix.logGamma = pmix.logAB.PK,
                Prior.beta    = Prior.beta,
                Prior.sd.PK = Prior.sd.PK,
                Prior.re.sd.PK = Prior.re.sd.PK
  )

  labels = list(intLabels = intLabels)
  Agents = Agent

  if(!Prior.Only){
    Data.1 = data.frame(t(do.call(rbind, Data)))[, -2]  # remove Dose.Grp because it jeopardize all simulation calculation
    if(is.null(Data$X)){colnames(Data.1) = c("DoseAdm", "Npat","Ntox", "Weight", "PK")}
     else{colnames(Data.1) =c("DoseAdm", "Npat","Ntox", "Weight", "PK","X")}
    N = c(Ncohorts = Ncohorts, Ndoses = Ndoses, Nmix1=Nmix1-1, Nmix2=Nmix2-1, Ncov=Ncov)

    Data.DLT = subset(Data.1, select=c(DoseAdm, Npat,Ntox, Weight))
    Data11 <- aggregate(Data.DLT, by=list(Data.DLT$DoseAdm, Data.DLT$Weight),sum)[,c(1,2,4,5)]

    Data1 <- Data11[,c(1,3,4,2)]
    colnames(Data1) = c("DoseAdm", "Npat","Ntox", "Weight")

    Data.PK =  subset(Data.1, select=c(DoseAdm, PK))
    Data12.1 <- aggregate(Data.PK, by=list(Data.PK$DoseAdm), mean)[,c(1,3)]
    Data12.2 <- aggregate(Data.PK, by=list(Data.PK$DoseAdm), sd)[,c(1,3)]
    Data2 <- merge( Data12.1, Data12.2, by="Group.1")
    colnames(Data2) = c("DoseAdm", "mean.PK","sd.PK")

  }

  if(Prior.Only){
    Data1 = NULL
    Data2 = NULL
    N = c(Ndoses = Ndoses, Nmix1=Nmix1-1, Nmix2=Nmix2-1, Ncov=Ncov)
  }

  #Collecting outputs

  outlist = list(Data = list(data = list(Data1, Data2), N = N, labels = labels, Agents = Agents),
                 Priors = Priors,
                 logAlphaBeta = logAlphaBeta,
                 logAlphaBeta.cor=logAlphaBeta.cor,
                 logAlphaEBeta = logAlphaEBeta,
                 logGamma = logGamma,
                 logGamma.cor=logGamma.cor,
                 logEGamma = logEGamma,
                 sd.PK     = sd.PK,
                 sd.re.PK  = sd.re.PK,
                 beta      = beta,
                 P = P,
                 Pcat = Pcat,
                 PcatP = PcatP,
                 pred.Tox = pred.Tox,
                 Post.pmix= Post.pmix,
                 Post.pmix.PK= Post.pmix.PK,
                 Rhat = list(Rhat = Rhat, Rhat.max = Rhat.max),
                 MCMC = list(pars.monitored = pars, MCMC = MCMC, RndSeed.R = Rnd.seed, RndSeed.WinBUGS = RndSeed.WinBUGS),
                 R2WB = fit)


  if(Plot){

    for(ll in 1:Nout){

    P1 <- outlist$P[,,ll]
    
    Plot.P(Est = P1, RowNames = Doses.Label)

    pdf(file = paste(Outfile,"_P","Xout_",ll,'.pdf',sep=''), onefile = F)

    Plot.P(Est = P1, RowNames = Doses.Label)

    dev.off()

    }

    #Pca1 = array(NA, dim=c(dim(outlist$Pcat)[1], dim(outlist$Pcat)[2]))
    for(kk in 1:Nout){

    Pcat1 <- outlist$Pcat[,,,kk]
     
    Plot.Pcat(Pcat = Pcat1, crit=int.crit, RowNames = Doses.Label)

    pdf(file = paste(Outfile,"_Pcat","Xout_",kk,'.pdf',sep=''), onefile = F)

    Plot.Pcat(Pcat = Pcat1, crit=int.crit, RowNames = Doses.Label)

    dev.off()

    }
  }


  #Creating Output file with necessary Outputs

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
  cat("Parameters:DLT model  \n\n")
  cat("\n")
  print(outlist$logAlphaBeta)
  cat("\n")
  print(outlist$logAlphaEBeta)
  cat("\n")
  cat("\n")
  cat("\n")
  cat("Parameters:PK model  \n\n")
  cat("\n")
  print(outlist$logGamma)
  cat("\n")
  print(outlist$logEGamma)
  cat("\n")
  cat("\n")
  cat("\n")
  cat("Covariates:PK Model \n\n")
  cat("\n")
  print(outlist$beta)
  cat("\n")
  cat("\n")
  cat("Variace component:PK Model \n\n")
  cat("\n")
  print(outlist$sd.PK)
  cat("\n")
  print(outlist$sd.re.PK)
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
  cat("Mixture Probabilities: DLT Parameters \n\n")
  cat("\n")
  print(outlist$Post.pmix)
  cat("\n")


  cat("\n")
  cat("Mixture Probabilities: PK Parameters \n\n")
  cat("\n")
  print(outlist$Post.pmix.PK)
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
