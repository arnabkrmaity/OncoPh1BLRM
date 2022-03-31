#' @title Meta analytic Predictive Prior
#'
#' @description This function derives Meta-analytic predictive prior for single agent Bayesian logstic linear model (BLRM) parameters for Phase I Oncology trials.
#' @param Agent name of the Agent under investigation.
#' @param Doses a vector of Doses for which inferential results will be provided.
#' @param DoseRef refernce dose for MAP.BLRM.
#' @param Data  external dose-DLT data from different strata: a list object contains DoseAdm, Npat, Ntox, Study and StudyLabel.
#' DoseAdm, Npat, and Ntox represent vectors of doses administered, Npat number of patients treated and the number of
#' dose-limiting toxicitites observed in different strata. Study is a vector indicating different strata. Here strata represents
#' study, subgroups and other source of external data. The lenght of DoseAdm, Npat, Ntox, and Study must be same. StudyLabel is 
#' a vector of strata names. The length of  StudyLabel needs to be unique number of elements in Study.
#' @param Prior.muA prior mean and sd of the normal distribution for the mean \eqn{\mu_{\log(\alpha)}} 
#' of the exchangeability distribution (see details).
#' @param Prior.muB prior mean and sd of the normal distribution for the mean \eqn{\mu_{\log(\beta)}} 
#' of the exchangeability distribution (see details).
#' @param Nexch Number of exchangeable set for between-strata standard deviations. Default value is 1 (full exchangeability model).
#' @param tau.index \eqn{(Nstrata+1)} vector with indices for between-strata standard deviations; indices
#' must be in \eqn{1.\ldots ,Nexch}. tau.index = 1 suffices for the full exchangeability model.
#' @param Prior.tauA \eqn{Nexch \times 2} matrix with mean (1st column) and sd (2nd column) 
#' of the log-normal prior(s) for \eqn{\tau_{\log(\alpha)}} of the exchangeability distribution; 
#' \eqn{Nexch} is the number of between-strata standard deviations (see tau.index and details).
#' @param Prior.tauB \eqn{Nexch \times 2} matrix with mean (1st column) and sd (2nd column) 
#' of the log-normal prior(s) for \eqn{\tau_{\log(\beta)}} of the exchangeability distribution; 
#' \eqn{Nexch} is the number of between-strata standard deviations (see tau.index and details).
#' @param Prior.rho vector with parameters of uniform distribution for correlation \eqn{\rho} of the exchangeability
#' distribution of \eqn{(\log(\alpha), \log(\beta))}(see details). The default is Uniform(-1,1).
#' @param Pcutoffs cutoffs to define the intervals for the true DLT probability. The default is
#' (0.16, 0.33) defining three categories: underdosing ({[}0--0.16{]}), target ((0.16--0.33)), and 
#' overdosing ({[}0.33--1.00{]}).
#' @param New.cohortsize number of subjects in a new cohort, for which the predictive distribution 
#' of the number of DLT's will be computed.
#' @param Outfile name of the output file.
#' @param pars MCMC parameters to save for output summary.
#' @param MCMC four component vector of MCMC parameters n.burnin, n.iterations, n.chains, n.thin.
#' See \link{R2WinBUGS} for details.
#' @param Warn logical: if TRUE, the convergence warning will be printed.
#' @param DIC logical: if TRUE, compute deviance, pD, and DIC. See \link{R2WinBUGS} for details.
#' @param progress.bar progress.bar type of progress bar. Possible values are “text”, “gui”, and “none”. Type “text”
#' is displayed on the R console. See \link{rjags} for details.
#' @param Rnd.seed Random seed for the R random stream. This is important for reproducibility. 
#' @param RSeed.WinBUGS Random seed for the WinBUGS random stream. This is important for reproducibility.
#'
#' @details Let \eqn{r_{d,s}} and \eqn{n_{d,s}} be the number of patients treated and the number of DLT's observed 
#' at dose level \eqn{d} from strata \eqn{s} .The aim is to derive a prior distribution 
#' for the logistic parameters of new stratum \eqn{(\log(\alpha^*),\log(\beta^*))}. The model specifications are 
#' as follows:
#' \deqn{r_{d,s} \sim Bin(\pi_{d,s},n_{d,s}); \; s=1, \ldots H} 
#' \deqn{\mbox{logit}(\pi_{d,s}) = \log(\alpha_s) + \beta_s \log(d/d^*)}
#' \deqn{\mbox{logit}({\pi_d}^*) = \log(\alpha^*) + \beta^* \log(d/d^*)}
#' \deqn{(\log(\alpha_s),\log(\beta_s))  \sim  BVN((\mu_{\log(\alpha)}, \mu_{\log(\beta)}), \Psi)}
#' \deqn{(\log(\alpha^*),\log(\beta^*))  \sim  BVN((\mu_{\log(\alpha)}, \mu_{\log(\beta)}), \Psi)}
#' 
#' where \eqn{\Psi}	is the between-strata covariance matrix with standard deviations \eqn{\tau_{\log(\alpha)}},
#' \eqn{\tau_{\log(\beta)}} and correlation \eqn{\rho}. The MAP prior of new stratum \eqn{(\log(\alpha^*),\log(\beta^*))}
#' is the predictive distribution and approximated by bivariate normal distribution:
#' \deqn{(\log(\alpha^*),\log(\beta^*))|(r_{d,s}, n_{d,s} : s= 1, \ldots, H)}
#' 
#' Furthermore, following prior distributions are used for model parameters:
#
#'  \itemize{
#'   \item normal priors are used for \eqn{\mu_{\log(\alpha)}} and \eqn{\mu_{\log(\beta)}}. 
#'   The expected toxicity at DoseRef can be corroborated by setting the mean appropriately.
#'   For example, when expected toxicity os 0.5 then the mean is log(0.5/(1 - 0.5)) which
#'   is default. 
#'   \item log-normal priors for \eqn{\tau_{\log(\alpha)}} and  \eqn{\tau_{\log(\beta)}} and 
#'   \item a uniform prior for \eqn{\rho}. 
#'}
#'
#'For differential discounting, different \eqn{\tau} can be used (see tau.index). For typical
#'values for the medians of the log-normal prior distributions, see Neuenschwander 2015.
#'
#' @references Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. A Bayesian Industry
#' Approach to Phase I Combination Trials in Oncology. Statistical Methods for Drug Combination
#' Studies. Taylor & Francis, 2015.
#' 
#' 
#' @return The function returns a list of the following main items:
#' \item{Data}{The dataset analyzed}
#' \item{Priors}{summary of prior specifications}
#' \item{MAP.logAB.sim}{MCMC sample of the predictive distribution of \eqn{(\log(\alpha^*),\log(\beta^*))}}
#' \item{MAP.Prior}{Bivariate normal approximation of MAP prior for \eqn{(\log(\alpha^*),\log(\beta^*))}}
#' \item{MAP.P}{summary of MAP prior for DLT probability in a new stratum}
#' \item{MAP.Pcat}{summary of interval probabilities for DLT probability in a new stratum}
#' \item{MAP.PcatP}{MAP.Pcat and MAP.P combined}
#' \item{logAB}{Posterior summaries BLRM parameters for all strata including new stratum (MAP.prior)}
#' \item{logAB.cor}{Estimated correlation for all strata including new stratum (MAP prior)}
#' \item{P}{summary of DLT probabilities for all strata including the new stratum (MAP prior)}
#' \item{Pcat}{summary of interval probabilities for DLT probability for all strata including the new stratum (MAP prior)}
#' \item{pred.Tox}{predicted probabilities for the number of DLT's for a cohort of size New.cohortSize for all strata including new stratum (MAP prior)}
#' \item{mu.tau.rho}{Posterior summary of the exchageability distribution and between-trial standard deviation}
#' \item{Rhat}{Gelman-Rubin MCMC convergence diagnostics}
#' \item{MCMC}{MCMC specifications}
#' \item{R2WB}{WinBUGS object with summaries of MCMC analysis (see \link{R2WinBUGS})}
#' In addition, a text output file (Outfile.txt) with all summaries is produced.
#' 
#' @examples 
#' \dontrun{
#' 
#'##Example 1: Prior derivation for Japanese study using Western data (example from 
#'section 6.3.3.3 in Neuenschwander 2015)
#'
#'# Western first-time-in-human trial is specified as the historical data
#'Data.western <- list(DosesAdm    = c(1, 2, 4, 8, 15, 30),
#'                      Npat       = c(5, 6, 5, 9,  8,  4),
#'                      Ntox       = c(0, 0, 0, 0,  1,  3),
#'                      Study      = c(1, 1, 1, 1,  1,  1),
#'                      StudyLabel = c("Western")
#'                      )
#'                      
#' # MAP Prior for Japanese Phase I trial
#' prior.japan <- MAP.BLRM(Doses = c(1, 2, 4, 8, 12, 15, 30, 40),
#'                DoseRef        = 10,
#'                Data           = Data.western,
#'                Prior.muA      = c(0,2),
#'                Prior.muB      = c(0,1),
#'                Prior.tauA     = c( log(0.25), log(2)/1.96),
#'                Prior.tauB     = c( log(0.125), log(2)/1.96),
#'                Pcutoffs       = c(0.16, 0.33),
#'                New.cohortsize = 3,
#'                Outfile        = "Output_Japan_MAP",
#'                MCMC           = c(2500, 7500, 4, 1),
#'                DIC            = TRUE,
#'                Rnd.seed       = 1452671,
#'                RSeed.WinBUGS  = 1367
#'                )
#'                
#' MAP.prior.japan <- prior.japan$MAP.Prior
#' MAP.prior.japan
#' 
#' ##Example 2: Differential discounting
#' 
#'}
#'
#'
#' @export




MAP.BLRM <- function(Agent              = NULL,
                     Doses              = NULL,
                     DoseRef            = NULL,
                     Data               = NULL,
                     Prior.muA          = c(0,2),
                     Prior.muB          = c(0,1),
                     Nexch              = 1,
                     tau.index          = NULL,
                     Prior.tauA         = NULL,
                     Prior.tauB         = NULL,
                     Prior.rho          = c(-1,1),
                     Pcutoffs           = c(0.16, 0.33),
                     New.cohortsize     = 3,
                     Outfile            = "Out",
                     pars               = c("mu.ex","logAB","pTox","pCat","MAP.pTox","MAP.pCat","pred.Tox","MAP.logAB","tau","rho","MAP.logAB.XY", "logAB.XY"),
                     MCMC               = c(2500, 7500, 4, 1),
                     Warn               = TRUE,
                     DIC                = TRUE,
                     progress.bar       = "text",
                     Rnd.seed           = 1452671,
                     RSeed.WinBUGS      = 1367
                     )
{

  version = "25-November-2017"

  set.seed(Rnd.seed)

  # MCMC parameters
  n.burnin = MCMC[1]
  n.iter = MCMC[2]
  n.chains = MCMC[3]
  n.thin = MCMC[4]

  #Reading the historical data

  DosesAdm   = Data$DosesAdm
  Ntox       = Data$Ntox
  Npat       = Data$Npat
  Study      = Data$Study
  StudyLabel = c(Data$StudyLabel, "MAP")
  
  if((length(DosesAdm) != length(Ntox)) || (length(DosesAdm) != length(Npat)) || 
     (length(DosesAdm) != length(Study)))
    stop("All the list items in Data object must have same length")

  Doses.Label = paste(Agent, "=", Doses, sep = '')

  Nobs <- length(Data$Ntox)
  Nstudy <- length(unique(Data$Study))
  Ncat <- length(Pcutoffs) + 1
  Ndoses <- length(Doses)

  Pcutoffs1 = c(0, Pcutoffs, 1)

  intLabels = paste(Pcutoffs1[1:Ncat], "-", Pcutoffs1[2:(Ncat + 1)], sep = "")
  Doses.Label = paste(Agent, "=", Doses, sep='')
  pred.tox.Label = paste("DLT=",0:New.cohortsize, sep='')

  DosesAdm[DosesAdm == 0] = 0.00001
  Doses[Doses == 0] = 0.00001

  DosesAdm = cbind(NULL, DosesAdm)
  Doses = cbind(NULL, Doses)
  Npat = cbind(NULL, Npat)
  Ntox = cbind(NULL, Ntox)
  Study =cbind(NULL, Study)


  if (Prior.rho[1] == -1)
    Prior.rho[1] <- -0.999
  if (Prior.rho[2] == 1)
    Prior.rho[2] <- 0.999

  #Between trial heterogeneity

  if(Nexch==1){tau.index = rep(1, Nstudy+1)}

  tau.index <- cbind(NULL, tau.index)

  Prior.tauA = rbind(Prior.tauA)
  Prior.tauB = rbind(Prior.tauB)

  model = MAP.BLRM.WB

  data = list("Nstudy", "Nobs", "Ndoses", "Nexch", "Ncat",
              "tau.index", "Doses", "DoseRef", "DosesAdm", "Npat", "Ntox", "Study",
              "Prior.muA","Prior.muB", "Prior.tauA", "Prior.tauB", "Prior.rho",
              "Pcutoffs", "New.cohortsize")

  inits.fun = function(i){
    return(list(rho = 0, mu.unconst = c(rnorm(1, Prior.muA[1], 0.1)/2, rnorm(1, Prior.muB[1], 0.1)),
                tau = cbind(exp(rnorm(Nexch, Prior.tauA[, 1], Prior.tauA[, 2]/2)), exp(rnorm(Nexch, Prior.tauB[, 1], Prior.tauB[, 2]/2)))))
                        }

  inits <- lapply(rep(1,n.chains),inits.fun)

  if(!is.element("MAP.logAB", pars)){
                 pars = c(pars, "MAP.logAB")
                                    }

  if(!is.element("MAP.logAB.XY", pars)){
    pars = c(pars, "MAP.logAB.XY")
  }

  if(!is.element("MAP.pCat", pars)){
    pars = c(pars, "MAP.pCat")
  }

  if(!is.element("MAP.pTox", pars)){
    pars = c(pars, "MAP.pTox")
  }

  if(!is.element("pred.Tox", pars)){
    pars = c(pars, "pred.Tox")
  }



  RndSeed.WinBUGS <- RSeed.WinBUGS

  fit = jags(
    data=data,
    inits=inits,
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
  fit$BUGSoutput$sims.matrix = NULL
  fit$BUGSoutput$sims.array = NULL

  #Collecting MCMC samples for MAP.logAB parameters

  logAB.MAP.sample <- fit$BUGSoutput$sims.list$MAP.logAB
  colnames(logAB.MAP.sample) <- c("logA", "logB")

  if(Warn){

  logAB.MAP.neff <- fit$BUGSoutput$summary[c("MAP.logAB[1]", "MAP.logAB[2]"), "n.eff"]
  if (any(logAB.MAP.neff < 1000)) {
    cat("\nWARNING: Neff < 1000 for logAB.MAP. Consider increasing the sample size (number of iterations).\n")
    message("WARNING: Neff < 1000 for logAB.MAP. Consider increasing the sample size (number of iterations).")
  }
}

#Collecting outputs

  logAB = NULL
  Rhat.logAB = NULL
  logAB.cor=NULL
  Rhat.logAB.cor=NULL
  P = NULL
  Rhat.P = NULL
  Pcat = NULL
  PcatP = NULL
  MAP.P =NULL
  MAP.Pcat =NULL
  MAP.PcatP =NULL
  mu = NULL
  Rhat.mu = NULL
  tau = NULL
  Rhat.tau = NULL
  rho = NULL
  Rhat.rho = NULL
  pred.Tox =NULL
  MAP.logAB.sim <- NULL

  logAB.MAP = BUGS.Select.Output("MAP.logAB", summary, cols = vnames)
  rownames(logAB.MAP) = c("MAP.logAlpha", "MAP.logBeta")
  logAB.MAP.XY = BUGS.Select.Output("MAP.logAB.XY", summary, cols = vnames)
  Rhat.logAB.MAP = logAB.MAP[, "Rhat"]
  MAP.logAB.cor = (logAB.MAP.XY[1] - logAB.MAP[1, 1] * logAB.MAP[2,1])/logAB.MAP[1, 2]/logAB.MAP[2, 2]

  MAP.Prior = c(logAB.MAP[, 1], logAB.MAP[, 2], MAP.logAB.cor)
  names(MAP.Prior)=c("MAP.logAlpha.mean", "MAP.logBeta.mean", "MAP.logAlpha.sd", "MAP.logBeta.sd", "MAP.AlphaBeta.cor")

  MAP.Pcat = BUGS.Select.Output("MAP.pCat", summary, cols = vnames)
  MAP.Pcat = BUGS.Table2Array(MAP.Pcat, Labels = list(Doses,intLabels), Dim = c(Ndoses, Ncat))
  MAP.Pcat = MAP.Pcat[, , "mean", drop = FALSE]
  MAP.Pcat = aperm(MAP.Pcat, c(1, 2, 3))

  MAP.P = BUGS.Select.Output("MAP.pTox", summary, cols = vnames)
  aiw = round(Beta.ms2ab(MAP.P[, "mean"], MAP.P[, "sd"])$n, 1)
  MAP.P = cbind(MAP.P, aiw)
  MAP.P = BUGS.Table2Array(MAP.P, Labels = list(Doses), Dim = c(Ndoses))
  MAP.P = aperm(MAP.P, c(1,2))

  MAP.PcatP <- cbind(MAP.Pcat[,,"mean"], MAP.P)

  if(is.element("logAB", pars)){
    logAB = BUGS.Select.Output("logAB", summary, cols = vnames)
    Rhat.logAB = logAB[, "Rhat"]

    logAB = BUGS.Table2Array(logAB, Labels = list(StudyLabel, c("logAlpha", "logBeta")), Dim = c(Nstudy+1, 2))
    logABl = list()
    for (j in 1:(Nstudy+1)){
      logABl[[j]] <- logAB[j, , ]
    }
    names(logABl) = StudyLabel
    logAB <- logABl

  }

  if(is.element("logAB.XY", pars) & is.element("logAB", pars)){

  logAB.XY = BUGS.Select.Output("logAB.XY", summary, cols = vnames)
  Rhat.logAB.cor = logAB.XY[, "Rhat"]

  logAB.cor <- list()

  for (j in 1:(Nstudy+1)){
  logAB.cor[[j]] <-  (logAB.XY[j, 1] - logAB[[j]][1, 1] * logAB[[j]][2,1])/logAB[[j]][1, 2]/logAB[[j]][2, 2]
  }
  names(logAB.cor) = StudyLabel

  }

  if (is.element("pTox", pars)) {
    P = BUGS.Select.Output("pTox", summary, cols = vnames)
    Rhat.P = P[, "Rhat"]

    aiw = round(Beta.ms2ab(P[, "mean"], P[, "sd"])$n, 1)
    P = cbind(P, aiw)
    P = BUGS.Table2Array(P, Labels = list(StudyLabel, Doses.Label), Dim = c(Nstudy+1, Ndoses))
    P = aperm(P, c(2,3,1))

  }


  if (is.element("pCat", pars)) {
    Pcat = BUGS.Select.Output("pCat", summary, cols = vnames)
    Pcat = BUGS.Table2Array(Pcat, Labels = list(StudyLabel, Doses.Label,intLabels), Dim = c(Nstudy+1, Ndoses, Ncat))
    Pcat = Pcat[, , , "mean", drop = FALSE]
    Pcat = aperm(Pcat, c(2, 3, 4, 1))

  }

  if (is.element("pCat", pars) & is.element("pTox", pars)) {

       PcatP = lapply(1:(Nstudy+1), function(i){cbind(Pcat[,,,i], P[,,i])})
       names(PcatP) <- StudyLabel

  }


  if (is.element("mu.ex", pars)){

    mu = BUGS.Select.Output("mu.ex", summary, cols = vnames)
    Rhat.mu = mu[, "Rhat"]

  }

  if (is.element("pred.Tox", pars)){

    pred.Tox = BUGS.Select.Output("pred.Tox", summary, cols = vnames)
    Rhat.pred.Tox= pred.Tox[, "Rhat"]

    pred.Tox = BUGS.Table2Array(pred.Tox, Labels = list(Doses.Label,pred.tox.Label), Dim = c(Ndoses, New.cohortsize+1))
    pred.Tox = pred.Tox[, , "mean", drop = FALSE]
    pred.Tox = aperm(pred.Tox, c(1, 2, 3))

  }

  if (is.element("tau", pars)) {
    tau = BUGS.Select.Output("tau", summary, cols = vnames)
    if (Nexch > 1)
      rownames(tau) = as.vector(outer(c("logAlpha exch",
                                        "logBeta  exch"), 1:Nexch, paste))
    else rownames(tau) = c("logAlpha", "logBeta")

      Rhat.tau = tau[, "Rhat"]
      names(Rhat.tau) = paste("tau.", rownames(tau), sep = "")

  }

  if (is.element("rho", pars)) {
    rho = BUGS.Select.Output("rho", summary, cols = vnames)
      Rhat.rho = rho["Rhat"]
  }



  Rhat = c(Rhat.logAB.MAP, Rhat.logAB, Rhat.P, Rhat.mu, Rhat.tau, Rhat.rho, Rhat.logAB.cor)
  Rhat.max = max(Rhat)

  if (Rhat.max > 1.1) {
  warning("\n Warning: at least one of the parameters has a Gelman-Rubin diagnostic Rhat > 1.1")
  }

  Priors = list(mu.logA = Prior.muA,
                mu.logB = Prior.muB,
                tau.logA = Prior.tauA,
                tau.logB = Prior.tauB,
                rho = Prior.rho)

  labels = list(StudyLabel = StudyLabel, intLabels = intLabels)
  N = c(Nobs = Nobs, Ndoses = Ndoses, Nstrata = Nstudy, Nexch = Nexch)

  Data1 = data.frame(t(do.call(rbind, Data)))
  colnames(Data1) = c("DoseAdm", "Npat","Ntox", "Study","Study Name")

#Collecting outputs

  outlist = list(Data = list(data = Data1, N = N, labels = labels),
                 Priors = Priors,
                 MAP.logAB.sim = logAB.MAP.sample,
                 MAP.logAB = logAB.MAP,
                 logAB.MAP.XY = logAB.MAP.XY,
                 MAP.Prior = MAP.Prior,
                 MAP.P = MAP.P,
                 MAP.Pcat = MAP.Pcat,
                 MAP.PcatP = MAP.PcatP,
                 logAB = logAB,
                 logAB.cor = logAB.cor,
                 P = P,
                 Pcat = Pcat,
                 PcatP = PcatP,
                 pred.Tox = pred.Tox,
                 mu.tau.rho = list(mu = mu, tau = tau, rho = rho),
                 Rhat = list(Rhat = Rhat, Rhat.max = Rhat.max),
                 MCMC = list(pars.monitored = pars, MCMC = MCMC, RndSeed.R = Rnd.seed, RndSeed.WinBUGS = RndSeed.WinBUGS),
                 R2WB = fit)


#Creating Output file with necessary

sink(paste(Outfile, ".txt", sep = ""), append = F)
  cat("\n --------------------------------------------------------------------------------")
  cat("Meta-Analytic-Predictive (MAP) Logistic Bayesian inference")
  cat("\n", date())
  cat("\n Working directory: ", getwd())
  cat("\n --------------------------------------------------------------------------------\n")
  cat("Data \n\n")
  print(outlist$Data)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("Meta-Analytic-Predictive Distribution \n\n")
  cat("\n")

  cat("Prior Distribution \n\n")
  print(outlist$MAP.Prior)
  cat("\n")
  cat("DLT rates \n\n")
  print(outlist$MAP.P)
  cat("\n")
  cat("Interval probabilities \n\n")
  print(outlist$MAP.Pcat)
  cat("\n")
  cat("Interval probabilities and DLT rates\n\n")
  print(outlist$MAP.PcatP)
  cat("\n")
  cat("Predictive probability for number of DLTs in new cohort \n\n")
  print(outlist$pred.Tox)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("Parameters of Exchangeability Distributions: mu, tau, rho \n\n")
  print(outlist$mu.tau.rho)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("Logistic parameters: by studies and MAP \n\n")
  print(outlist$logAB)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("DLT rates: by studies and MAP \n\n")
  print(outlist$P)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("Interval probabilities: by-studies and MAP \n\n")
  print(outlist$Pcat)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("DLT rates and interval probabilities: by-studies and MAP \n\n")
  print(outlist$PcatP)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("Gelman-Rubin diagnostics \n\n")
  print(outlist$Rhat)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("Prior distributions \n\n")
  print(outlist$Priors)
  cat("\n")
  cat("\n --------------------------------------------------------------------------------\n")
  cat("MCMC parameters \n\n")
  print(outlist$MCMC)
  cat("\n")
  sink()


  return(outlist)

}

