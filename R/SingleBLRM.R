#' @title SingleBLRM Function 
#'
#' @description This function performs analysis of dose limiting toxicity (DLT) data with Bayesian logistic regression (BLRM) for 
#' single agent Phase I Oncology trial
#'
#' @param Data dose-DLT data: a list object contains DoseAdm, Npat, Ntox, and Weight. DoseAdm, Npat, 
#'  and Ntox are vectors of doses administered, Npat number of patients treated and the number of
#' dose-limiting toxicitites observed per cohort. The vector Weight represents the relevance of the cohort
#' for current trial. If Weight remains unspecified, then it will equal 1 for each cohort.
#' @param Agent name of the Agent under investigation.
#' @param Doses a vector of Doses for which inferential results will be provided.
#' @param DoseRef refernce dose for BLRM.
#' @param Prior.logAB prior for the model parameters Bayesian Logistic Regession Model (a list object).
#' each element of the list is vector with five elements. 
#' (\eqn{m_\alpha}, \eqn{m_\beta}, \eqn{s_\alpha}, \eqn{s_\beta}, \eqn{Corr(\alpha,\beta)}). These five 
#' parameters represent prior mean, sd and correlation of the bivariate normal prior for 
#' (\eqn{\log(\alpha)}, \eqn{\log(\beta)}). If no mixture prior is used, the list contains a single vector. 
#' Otherwise, the list  contains a number of five component vectors. The number of vectors are equal to number
#' of mixture components.   
#' @param pmix.logAB  a vector with apriori mixture probabilities. The number of components of pmix.logAB.
#' must be same as number of components in Prior.logAB.See the package vignettes for details.
#' @param Prior.Only logical: TRUE means data is ignored; if FALSE (the default), then data is taken
#' into account.  
#' @param Pcutoffs cutoffs to define the intervals for the true DLT probability. The default is
#' (0.16, 0.33) defining three categories: underdosing ({[}0--0.16{]}), target ((0.16--0.33)), and 
#' overdosing ({[}0.33--1.00{]}).
#' @param New.cohortsize number of subjects in a new cohort, for which the predictive distribution 
#' of the number of DLT's will be computed.
#' @param Outfile name of the output file.
#' @param Plot logical: if TRUE (default) plots will be generated.
#' @param Table logical: if TRUE (default) the output tables will be generated.
#' @param int.crit vector of critical values for each interval defined by Pcutoffs (only used for
#' graphical output).
#' @param pars MCMC parameters to save for output summary.
#' @param MCMC four component vector of MCMC parameters n.burnin, n.iterations, n.chains, n.thin.
#' See \link{R2WinBUGS} for details.
#' @param Warn logical: if TRUE, the convergence warning will be printed.
#' @param DIC logical: if TRUE, compute deviance, pD, and DIC. See \link{R2WinBUGS} for details.
#' @param progress.bar type of progress bar. Possible values are “text”, “gui”, and “none”. Type “text”
#' is displayed on the R console. See \link{rjags} for details.
#' @param Rnd.seed Random seed for the R random stream. This is important for reproducibility. 
#' @param RSeed.WinBUGS Random seed for the WinBUGS random stream. This is important for reproducibility. 
#'
#' @details Statistical model:
#' \deqn{r_j|n_j \sim \mbox{Binomial} (\pi(d_j), n_j)} 
#' \deqn{logit(\pi(d_j)) = \log (\alpha) + \beta \log (d_j/d^*).}
#'
#' Here, \eqn{r_j} and \eqn{n_j} are number of patients treated and number of DLT's observed 
#' at dose level \eqn{d_j}. Moreover, \eqn{\alpha} is the odds, \eqn{\pi/(1 - \pi} of a DLT at \eqn{d^*}, an 
#' arbitrary scaling dose. \eqn{\beta > 0} is the increase in the log-odds of a DLT by a unit increase in 
#' log-dose. For instance, doubling the dose from \eqn{d} to \eqn{2d} will increase the odds of a DLT by 
#' the factor \eqn{2^\beta}, from \eqn{\pi_d/(1 - \pi_d)} to \eqn{2^\beta \pi_d/( 1- \pi_d)}. 
#' Neuenschwander et al. 2008 have suggested bivariate normal prior for \eqn{(\log(\alpha), \log(\beta))}.  
#' 
#' \deqn{(\log(\alpha), \log(\beta)) \sim N_2 (m_\alpha, m_\beta, s_\alpha, s_\beta, Corr(\alpha,\beta))} 
#'
#' Furthermore, mixture priors can be used for \eqn{(\log(\alpha), \log(\beta))},
#' \deqn{(\log(\alpha), \log(\beta)) \sim \sum \omega_i N_2 ({m_\alpha}_i, {m_\beta}_i, {s_\alpha}_i, {s_\beta}_i, Corr(\alpha_i,\beta_i)); \;\; \sum \omega_{i} = 1}
#'
#' Details on prior derivation can be found in the book chapter by Neuenschwander et al 2015.
#'
#'@seealso \code{\link{SingleBLRMSim}}, \code{\link{MAP.BLRM}}.
#'
#'@references 
#'Neuenschwander B, Branson M, Gsponer T. Critical aspects of the Bayesian approach to phase I
#'cancer trials. Stat Med. 2008;27(13):2420-39.
#' 
#'Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. A Bayesian Industry
#'Approach to Phase I Combination Trials in Oncology. Statistical Methods for Drug Combination
#'Studies. Taylor & Francis, 2015.
#'
#'@return 
#'The function returns as a list the following main items:
#' \item{data}{The dataset analyzed}
#' \item{Priors}{summary of prior specifications}
#' \item{logAlphaBeta}{Posterior summary of (\eqn{\log(\alpha)}, \eqn{\log(\beta)})}
#' \item{logAlphaBeta.cor}{Posterior summary of the correlation between (\eqn{\log(\alpha)}, \eqn{\log(\beta)})}
#' \item{logAlphaEBeta}{Posterior summary of (\eqn{\log(\alpha)}, \eqn{\beta})}
#' \item{P}{Posterior summary of DLT probabilities}
#' \item{Pcat}{Posterior summary of DLT interval probabilities}
#' \item{PcatP}{P and Pcat combined}
#' \item{pred.Tox}{Predicted probabilities for the number of DLT's for a cohort of size New.cohortSize}
#' \item{Post.pmix}{Posterior mixture weigths (for mixture prior)}
#' \item{Rhat}{Gelman-Rubin MCMC convergence diagnostics}
#' \item{MCMC}{MCMC specifications}
#' \item{R2WB}{WinBUGS object with summaries of MCMC analysis (see \link{R2WinBUGS})}
#' 
#' In addition the following files will be created:
#' \item{Outfile}{a text output file with summary outputs}
#' \item{Outfile_P}{estimates and 95\% interval plots for DLT rates (in pdf format)}
#' \item{Outfile_Pcat}{interval plots for DLT rates (in pdf format)}
#' 
#' @examples 
#' \dontrun{
#' ## Example 1: Analysis of phase I single agent data: prior only 
#' 
#' fit.BLRM.single.P <- SingleBLRM(Data           = NULL,
#'                                 Agent          = "DRUG 1",
#'                                 Doses          = c(1, 2, 4, 8, 15, 30, 40), 
#'                                 DoseRef        = 10,
#'                                 Prior.logAB    = list(c(log(0.5), 0, 2, 1, 0)),
#'                                 pmix.logAB     = c(1),
#'                                 Prior.Only     = TRUE,
#'                                 Pcutoffs       = c(0.16, 0.33),
#'                                 New.cohortsize = 3,
#'                                 Outfile        = "Output_single_prior",
#'                                 MCMC           = c(5000, 10000, 4, 1),
#'                                 DIC            = FALSE,
#'                                 Rnd.seed       = 1452671,
#'                                 RSeed.WinBUGS  = 1367
#'                                 )
#'    
#'    
#'    
#' ## Example 2: Analysis of phase I single agent data: analysis of 6 cohort data 
#' 
#' fit.BLRM.single <- SingleBLRM(Data = list(DosesAdm = c(1, 2, 4, 8, 15, 30),
#'                                           Npat     = c(5, 6, 5, 9,  8,  4),
#'                                           Ntox     = c(0, 0, 0, 0,  1,  3),
#'                                           Weight   = c(1, 1, 1, 1,  1,  1)
#'                                           ),
#'                               Agent          = "DRUG 1",
#'                               Doses          = c(1, 2, 4, 8, 15, 30, 40), 
#'                               DoseRef        = 10,
#'                               Prior.logAB    = list(c(log(0.5), 0, 2, 1, 0)),
#'                               pmix.logAB     = c(1),
#'                               Prior.Only     = FALSE,
#'                               Pcutoffs       = c(0.16, 0.33),
#'                               New.cohortsize = 3,
#'                               Outfile        = "Output_single",
#'                               MCMC           = c(5000, 10000, 4, 1),
#'                               DIC            = FALSE,
#'                               Rnd.seed       = 1452671,
#'                               RSeed.WinBUGS  = 1367
#'                               )
#'          
#'                               
#'                                                                         
#'
#' ## Example 3: Analysis of phase I single agent data: analysis of 6 cohort data 
#' ##            together with down-weighted data from other trial#' 
#' 
#' fit.BLRM.single.D <- SingleBLRM(Data = list(DosesAdm = c(1, 2, 4, 8, 15, 30,   8,  15),
#'                                             Npat     = c(5, 6, 5, 9,  8,  4,   6,   6),
#'                                             Ntox     = c(0, 0, 0, 0,  1,  3,   0,   2),
#'                                             Weight   = c(1, 1, 1, 1,  1,  1, 0.5, 0.5)
#'                                             ),
#'                                 Agent          = "DRUG 1",
#'                                 Doses          = c(1, 2, 4, 8, 15, 30, 40), 
#'                                 DoseRef        = 10,
#'                                 Prior.logAB    = list(c(log(0.5), 0, 2, 1, 0)),
#'                                 pmix.logAB     = c(1),
#'                                 Prior.Only     = FALSE,
#'                                 Pcutoffs       = c(0.16, 0.33),
#'                                 New.cohortsize = 3,
#'                                 Outfile        = "Output_single_dwght",
#'                                 MCMC           = c(5000, 10000, 4, 1),
#'                                 DIC            = FALSE,
#'                                 Rnd.seed       = 1452671,
#'                                 RSeed.WinBUGS  = 1367
#'                                 )
#'                                 
#'        
#'                                 
#'                                                          
#' ## Example 4: Analysis of phase I single agent data with a mixture prior.
#' ##            A 0.5-0.5 mixture of prior distribution is applied
#' 
#' fit.mixture <- SingleBLRM(Data = list(DosesAdm = c(1, 2, 4, 8, 15, 30),
#'                                       Npat     = c(5, 6, 5, 9,  8,  4),
#'                                       Ntox     = c(0, 0, 0, 0,  1,  3),
#'                                       Weight   = c(1, 1, 1, 1,  1,  1)
#'                                       ),
#'                           Agent          = "DRUG 1",
#'                           Doses          = c(1, 2, 4, 8, 15, 30, 40), 
#'                           DoseRef        = 10,
#'                           Prior.logAB    = list(c(log(0.5), 0, 2, 1, 0),
#'                                                 c(0, 0, 1, 1, 0.5)),
#'                           pmix.logAB     = c(0.5, 0.5),
#'                           Prior.Only     = FALSE,
#'                           Pcutoffs       = c(0.16, 0.33),
#'                           New.cohortsize = 3,
#'                           Outfile        = "Output_mixture",
#'                           MCMC           = c(5000, 10000, 4, 1),
#'                           DIC            = FALSE,
#'                           Rnd.seed       = 1452671,
#'                           RSeed.WinBUGS  = 1367
#'                           )                               
#'}
#'
#'
#' @export


SingleBLRM <- function(Data               = NULL,
                       Agent              = NULL,
                       Doses              = NULL,
                       DoseRef            = NULL,
                       Prior.logAB        = NULL,
                       pmix.logAB         = c(1),
                       Prior.Only         = FALSE,
                       Pcutoffs           = c(0.16, 0.33),
                       New.cohortsize     = 3,
                       Outfile            = "Out",
                       Plot               = TRUE,
                       Table              = TRUE,
                       int.crit           = c(1,1,0.25),
                       pars               = c("logAlphaBeta", "pTox","pCat","pred.Tox"),
                       MCMC               = c(2500, 7500, 4, 1),
                       Warn               = TRUE,
                       DIC                = TRUE,
                       progress.bar       = "text",
                       Rnd.seed           = 1452671,
                       RSeed.WinBUGS      = 1367
)
{


  version = "15-May-2019"
  
  if(is.null(Data) && (Prior.Only == FALSE))
  stop("Either provide data OR set Prior.Only = TRUE")
  
  stopifnot("Mixture weights must add to 1" = sum(pmix.logAB) == 1)
  
  set.seed(Rnd.seed)

  if(Prior.Only){DIC=FALSE}

  # MCMC parameters
  n.burnin = MCMC[1]
  n.iter = MCMC[2]
  n.chains = MCMC[3]
  n.thin = MCMC[4]

  #Reading Data

  if(!Prior.Only){

    DosesAdm   = Data$DosesAdm
    Ntox       = Data$Ntox*Data$Weight
    Npat       = Data$Npat*Data$Weight
    
    if((length(DosesAdm) != length(Ntox)) || (length(DosesAdm) != length(Npat)) || 
       (length(DosesAdm) != length(Data$Weight)))
      stop("All the list items in Data object must have same length")
    
    
    DosesAdm[DosesAdm ==0] = 0.00001

    Ncohorts <- length(DosesAdm)

    DosesAdm<- cbind(NULL, DosesAdm)

    Ntox <- cbind(NULL, Ntox)
    Npat <- cbind(NULL, Npat)
    zero <- rep(0,Ncohorts)

  }

  if(Prior.Only){Npat=0; Ntox=0}

  Ndoses <- length(Doses)

  Ncat <- length(Pcutoffs) + 1
  Pcutoffs1 = c(0, Pcutoffs, 1)

  intLabels = paste(Pcutoffs1[1:Ncat], "-", Pcutoffs1[2:(Ncat + 1)], sep = "")
  Doses.Label = paste(Agent, "=", Doses, sep = '')
  pred.tox.Label = paste("DLT=",0:New.cohortsize, sep='')

  Doses.1 =Doses

  Doses[Doses ==0] = 0.00001

  Doses <- cbind(NULL, Doses)

  if(min(pmix.logAB) < 1){pars <- c(pars, "Zcat")}

  #Priors

  #Bivariate normal prior for (logAlpha,logBeta)

  Prior.logAB.1 <- c(Prior.logAB, list(c(0,0,0.1,0.1,0)))
  Pmix <- c(pmix.logAB, 0)
  Prior <- do.call(rbind,Prior.logAB.1)

  Pmix <- cbind(NULL, Pmix)
  Nmix <- length(Prior.logAB.1)


  #Choice of model

  if(Prior.Only){model= SingleBLRM.Prior.WB; model.index=1}
  if(!Prior.Only & (any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)))){model= SingleBLRM.Weight.WB; model.index=2}
  if(!Prior.Only & !(any(Npat <= 1) || any(Ntox != round(Ntox)) || any(Npat != floor(Npat)))){model= SingleBLRM.WB; model.index=3}

  #Data for jags
  if(Prior.Only){

    data= list("Ndoses", "Ncat", "Nmix",
               "Prior", "Pmix",
               "Doses", "DoseRef",
               "Pcutoffs", "New.cohortsize")
  }
  if(!Prior.Only & model.index== 3){

    data= list("Ncohorts", "Ndoses", "Ncat", "Nmix",
               "Prior", "Pmix",
               "Doses", "DoseRef",
               "DosesAdm", "Ntox", "Npat",
               "Pcutoffs", "New.cohortsize")
  }

  if(!Prior.Only & model.index== 2){

    data= list("Ncohorts", "Ndoses", "Ncat", "Nmix",
               "Prior", "Pmix",
               "Doses", "DoseRef",
               "DosesAdm", "Ntox", "Npat", "zero",
               "Pcutoffs", "New.cohortsize")
  }


  #Initial values
  inits.fun = function(i){

    logAlphaBetaAll <- matrix(NA, nrow=Nmix, ncol=2)

    for(i1 in 1:Nmix){logAlphaBetaAll[i1,1:2] <- mvtnorm::rmvnorm(1, sigma = diag(2))}

    Z <- 1

    return(list(logAlphaBetaAll = logAlphaBetaAll,
                Z=Z
    )
    )
  }

  inits <- lapply(rep(1,n.chains),inits.fun)

  RndSeed.WinBUGS <- RSeed.WinBUGS

  pars = c(pars,"logAlphaEBeta", "logAlphaBeta.XY")


  if(DIC){rjags::load.module('dic')}

  fit = R2jags::jags(
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
    message("WARNING: Neff < 1000 for logAlphaBeta1\ Consider increasing the sample size (number of iterations).")
  }

         }
  logAlphaBeta = NULL
  Rhat.logAlphaBeta = NULL
  logAlphaBeta.cor=NULL
  Rhat.logAlphaBeta.cor=NULL
  logAlphaEBeta =NULL
  Rhat.logAlphaEBeta =NULL


  P = NULL
  Rhat.P = NULL

  Pcat = NULL
  PcatP = NULL
  pred.Tox =NULL
  Post.pmix = NULL


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

  if (is.element("pTox", pars)) {
    P = BUGS.Select.Output("pTox", summary, cols = vnames)
    Rhat.P = P[, "Rhat"]
    aiw = round(Beta.ms2ab(P[, "mean"], P[, "sd"])$n, 1)
    P = cbind(P, aiw)
    P = BUGS.Table2Array(P, Labels = list(Doses.Label), Dim = Ndoses)
    #P = aperm(P, c(1,3,2))

  }


  if (is.element("pCat", pars)) {
    Pcat = BUGS.Select.Output("pCat", summary, cols = vnames)
    Pcat = BUGS.Table2Array(Pcat, Labels = list(Doses.Label, intLabels), Dim = c(Ndoses, Ncat))
    Pcat = Pcat[, , "mean", drop = FALSE]
    #Pcat = aperm(Pcat, c(1, 3, 2, 4))

  }

 if (is.element("pCat", pars) & is.element("pTox", pars)){

       PcatP <- cbind(rbind(NULL,Pcat[,,1]), rbind(NULL,P))


    vnames1 = c(intLabels,vnames,"aiw")

    dimnames(PcatP) = list(Doses.Label, vnames1)
  }


  if (is.element("pred.Tox", pars)){

    pred.Tox = BUGS.Select.Output("pred.Tox", summary, cols = vnames)
    Rhat.pred.Tox= pred.Tox[, "Rhat"]



    pred.Tox = BUGS.Table2Array(pred.Tox, Labels = list(Doses.Label, pred.tox.Label), Dim = c(Ndoses, New.cohortsize+1))
    pred.Tox = pred.Tox[, ,"mean", drop = FALSE]

  }

  if (is.element("Zcat", pars)){

    Post.pmix.1 = BUGS.Select.Output("Zcat", summary, cols = "mean")
    Post.pmix <- cbind(NULL,Post.pmix.1[-Nmix])
    Post.pmix.nms <- paste("mix",1:(Nmix-1),sep='')
    rownames(Post.pmix) <- Post.pmix.nms

  }


  Rhat = c(Rhat.logAlphaBeta, Rhat.logAlphaBeta.cor, Rhat.logAlphaEBeta, Rhat.P)

  Rhat.max = max(Rhat)

  if (Rhat.max > 1.1) {
    warning("\n Warning: at least one of the parameters has a Gelman-Rubin diagnostic Rhat > 1.1")
  }

  Priors = list(Prior.logAB = Prior.logAB,
                pmix.logAB = pmix.logAB
  )

  labels = list(intLabels = intLabels)
  Agents = Agent

  if(!Prior.Only){
    Data.1 = data.frame(t(do.call(rbind, Data)))
    colnames(Data.1) = c("DoseAdm", "Npat","Ntox", "Weight")
    N = c(Ncohorts = Ncohorts, Ndoses = Ndoses, Nmix=Nmix-1)

    Data11 <- stats::aggregate(Data.1, by=list(Data.1$DoseAdm, Data.1$Weight),sum)[,c(1,2,4,5)]
    Data1 <- Data11[,c(1,3,4,2)]
    colnames(Data1) = c("DoseAdm", "Npat","Ntox", "Weight")

  }

  if(Prior.Only){
    Data1 = NULL
    N = c(Ndoses = Ndoses, Nmix=Nmix-1)
  }

  #Collecting outputs

  outlist = list(Data = list(data = Data1, N = N, labels = labels, Agents = Agents),
                 Priors = Priors,
                 logAlphaBeta = logAlphaBeta,
                 logAlphaBeta.cor=logAlphaBeta.cor,
                 logAlphaEBeta = logAlphaEBeta,
                 P = P,
                 Pcat = Pcat,
                 PcatP = PcatP,
                 pred.Tox = pred.Tox,
                 Post.pmix= Post.pmix,
                 Rhat = list(Rhat = Rhat, Rhat.max = Rhat.max),
                 MCMC = list(pars.monitored = pars, MCMC = MCMC, RndSeed.R = Rnd.seed, RndSeed.WinBUGS = RndSeed.WinBUGS),
                 R2WB = fit)

   if(Plot == TRUE){

     P1 <- outlist$P

     grDevices::pdf(file = paste(Outfile,"_P",'.pdf',sep=''), onefile = F)

     Plot.P(Est = P1, RowNames = Doses)

     grDevices::dev.off()
     
     Plot.P(Est = P1, RowNames = Doses)  # this will plot locally

     #Pca1 = array(NA, dim=c(dim(outlist$Pcat)[1], dim(outlist$Pcat)[2]))
     Pcat1 <- outlist$Pcat[,,1]

     grDevices::pdf(file = paste(Outfile,"_Pcat",'.pdf',sep=''), onefile = F)

     Plot.Pcat(Pcat = Pcat1, crit=int.crit, RowNames = Doses)

     grDevices::dev.off()
     
     Plot.Pcat(Pcat = Pcat1, crit=int.crit, RowNames = Doses)  # this will plot locally


   }


  #Creating Output file with necessary Outputs
  
  if(Table == TRUE)
  {
    
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
    cat("Parameters  \n\n")
    cat("\n")
    print(outlist$logAlphaBeta)
    cat("\n")
    print(outlist$logAlphaEBeta)
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
    print(outlist$Post.pmix)
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
    
  }

  
  return(outlist)

}
