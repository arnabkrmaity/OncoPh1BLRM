#options(warn=-1)

#' @title SingleBLRMSim Function
#'
#' @description This function performs Phase I Oncology trial simulation for single agent using Bayesian 
#' Logistic Regression Model (BLRM) and Escalation with Overdose Criteria (EWOC)
#'
#' @param Design   Design specification: contains required parameters for \code{\link{SingleBLRM}}. This incudes provisional
#' dose levels, reference dose, and prior specifications. See \code{\link{SingleBLRM}} for more details.
#' @param nTrials Number of simulated trials.
#' @param TrueProbs a vector of true toxicity  probabilities for different dose levels.
#' @param CohortRnd Cohort size (Random). A list with the first element a vector of the available
#'  cohort size(s), and the second element a vector of the same length with the respective selection
#'  probability of the cohort size. If a fixed cohort size is desired use e.g. CohortRnd = list(c(3), c(1))
#' @param Histdatf logical: if TRUE, historical data is used (optional)
#' @param Histdata a list with named elements  DosesAdm1, DosesAdm2, Npat, Ntox, Weight containing the
#'  historical data (optional)
#' @param Pcutoffs cutoffs to define the intervals for the true DLT probability. The default is
#' (0.16, 0.33) defining three categories: underdosing ({[}0--0.16{]}), target ((0.16--0.33)), and 
#' overdosing ({[}0.33--1.00{]}).
#' @param Nsim four component vector of MCMC parameters n.burnin, n.iterations, n.chains, n.thin.
#' See \link{R2WinBUGS} for details.
#' @param Decision Function which decides after each simulated cohort how to proceed
#' @param Outfile A R object saves all simulations output objects
#' @param Out.summary A txt file with simulation summary
#' @param RndSeed Seed used for the random number generator. Required parameter
#'
#' @details See vignettes. browseVignettes("OncoPh1BLRM")
#' 
#' @references 
#' #'Neuenschwander B, Branson M, Gsponer T. Critical aspects of the Bayesian approach to phase I
#'cancer trials. Stat Med. 2008;27(13):2420-39.
#' 
#' Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. A Bayesian Industry
#' Approach to Phase I Combination Trials in Oncology. Statistical Methods for Drug Combination
#' Studies. Taylor & Francis, 2015.
#' 
#' @return Input arguments are 
#' \item{Design}{The input data and prior information}
#' \item{nTrials}{Number of simulated trials}
#' \item{TrueProbs}{True toxicity vector}
#' \item{StartDose}{Start dose for the escalation trial}
#' \item{Pcutoffs}{Target toxicity interval}
#' \item{CohortRnd}{}
#' \item{Decision}{}
#' \item{Outfile}{Name of the output file}
#' \item{Out.summary}{Name of the summary output file}
#' \item{RndSeed}{The seed which was fixed to reproduce the result}
#' 
#' Output arguments are 
#' \item{DLTMATRX}{the occurrences of DLT's per dose per trial}
#' \item{MTDSIM}{a binary matrix which is only TRUE when MTD has been observed in the trial}
#' \item{NNMATRX}{the cohort sizes per dose per trial}
#' \item{res.sim$simul}{descrioption of each trial -- Trial no. (Trial), Cohort no. (Cohort), 
#' Number of toxicity (Ntox), Number of patients (Npat), Weight, Current dose (DoseAdm), 
#' underdose probability (under), probability of target toxicity (target), 
#' overdose probability (over), probability of the current dose (true.prob), 
#' suggestion for the next dose, if any, (nextDoseAdm), underdose probability of the next 
#' dose (nextunder), target toxicity probability of the next dose (nexttarget), 
#' probability of the next dose (nextprob), a binary variable which is only 1 (TRUE) 
#' if the current dose is declared as the MTD, if the first dose is so toxic that 
#' \eqn{\Pr(\pi_d > 0.33) \ge EWOC} then the trial is stooped and too.toxic variable is set to 1, 
#' the Rhat.max information, crash = 1 if the code is crashed, the seed number for R (Rseed) 
#' and JAGS (WBseed)}
#' \item{res.sim$simsummary}{description of those cohorts within those trials 
#' in which the corresponding dose is the MTD}
#' \item{sim.sum}{simulation summary -- Proportion of patients in target dose, 
#' Proportion of patients in over dose, Proportion of patients in under dose,
#' Proportion of trials with MTD within target dose region, 
#' Proportion of trails with MTD within over dose region, 
#' Proportion of trails with MTD within under dose region, 
#' Proportion of trials stopped due to too toxicity at the start dose, 
#' Average sample size}
#' 
#' In addition, an R object output and a text output are generated containing all of the above. 
#' 
#' @seealso \code{\link{SingleBLRM.PK}}
#' 
#' @examples 
#' \dontrun{
#' 
#' require(boot)   # for logit
#' require(doRNG)  # for parallel running
#' 
#' # Set a PK profile
#' True.PK  <- c(0.05239955, 0.12834821, 0.27734375, 0.67438616,
#'               0.75083705, 1.09784226, 1.14880952, 1.15364583, 1.16316964)
#' # Obtain a true troxicity profile using invers logit tranformation                                              
#' TrueProbs <- inv.logit(True.PK - 2)  # move back to TrueProbs to generate the DLT data
#'
#'
#'
#' # Design Parameters
#' simulation.design <- list(Agent          = "PTK7",
#'                           Doses          = c(0.133333, 0.3333333, 0.8333333, 1.4, 1.8666667, 
#'                                              2.1, 2.4666667, 2.8, 3.2),
#'                           Dose.Grp       = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
#'                           DoseRef        = 3.7, 
#'                           Prior.logAB    = list(c(log(0.5), 1, 1, 1, 0)),  
#'                           pmix.logAB     = c(1),
#'                           PK             = True.PK,
#'                           PK.sd          = rep(0.5, length(True.PK)),
#'                           Prior.logAB.PK = list(c(0, 0, 1, 1, 0)),
#'                           pmix.logAB.PK  = c(1)
#'                           )
#' 
#' 
#'   
#' # Running Simulation
#' sim.args <- list(Design       = simulation.design,
#'                  nTrials      = 10,
#'                  TrueProbs    = TrueProbs,
#'                  StartDose    = 0.133333,
#'                  Pcutoffs     = c(0.16, 0.33),
#'                  CohortRnd    = list(c(3), c(1)),
#'                  Decision     = SingleBLRMPKSimDecisionRule(MaxN = 50, mnMTD = 6, mnTOT = 12,
#'                                                             mnTAR = 0.5, EWOC = 0.25,
#'                                                             Escalation = "maximal", 
#'                                                             doseMaxEsc = c(3)
#'                                                             ),
#'                  Outfile      = "simulation_blrmpk",
#'                  Out.summary  = "summary_simulation_blrmpk.txt",
#'                  RndSeed      = 123   # SEED2: To CONTROL REPRODUCIBILITY
#'                  )
#' 
#' # Setting up back-end for parallel computing (DO NOT CHANGE)
#' library(doParallel)
#' nCores <- detectCores()  # number of cores
#' cl <- makeCluster(nCores) 
#' registerDoParallel(cl)
#' registerDoRNG(seed = 123, once = FALSE) #SEED1 : To initiate enviorment. 
#' # THIS IS THE IMPORTANT STEP FOR PARALLELISM AND REPRODUCIBILITY
#' getDoParWorkers()
#' 
#' simulation.blrmpk <- do.call(SingleBLRMPKSim, sim.args)
#' 
#' # Stop all clusters (DO NOT CHANGE)
#' OncoPh1BLRM:::unregister()
#' stopImplicitCluster() 
#' 
#' simulation.blrmpk$sim.sum  # summary of the simulated trials
#' }
#' 
#' 
#' @export


SingleBLRMPKSim <- function(Design       = NULL,
                            nTrials      = 1000,
                            TrueProbs    = NULL,
                            StartDose    = NULL,
                            CohortRnd    = list(c(  3,   4,    5,    6),
                                                c(0.8, 0.1, 0.05, 0.05)),
                            Histdatf     = 0,
                            Histdata     = NULL,
                            Pcutoffs     = c(0.16, 0.33),
                            Nsim         = c(1000, 2000, 2, 1),
                            Decision     = SingleBLRMSimDecisionRule(MaxN = 60, mnMTD = 6, 
                                                                     mnTOT = 15, mnTAR = 0.5, 
                                                                     EWOC = 0.25, 
                                                                     Escalation = "maximal", 
                                                                     doseMaxEsc = 2),
                            Outfile      = "Sims",
                            Out.summary  = "sims_summary.txt",
                            RndSeed      = 12
)
{

  set.seed(RndSeed)
  SeedParam <- RndSeed

  #SeedParam <- sample.int(10, nTrials)

  suppressPackageStartupMessages(require(foreach, warn.conflicts=FALSE, quiet=TRUE))
  suppressPackageStartupMessages(require(doRNG, warn.conflicts=FALSE, quiet=TRUE))

  if(missing(SeedParam)) {
    stop("ERROR: Random seed must be set for this operation!\nPlease specify SeedParam.")
  }

  forArgs <- list(i = 1:nTrials, .packages=c("methods","R2WinBUGS","rjags","R2jags"), .options.RNG=SeedParam)

  if(getDoParRegistered()) {
    ## for the MPI backend set the chunksize to be at most 10
    if(getDoParName() == "doMPI") {
      forArgs$.options.mpi <- list(chunkSize=min(10, ceiling(nTrials/getDoParWorkers())))
      cat("Setting doMPI backend chunkSize to", forArgs$.options.mpi$chunkSize, "\n")
    }
  }
  

  

  Results <- do.call(foreach, forArgs) %dorng% {
    
    # parallel.sim <- function(index){
      
      
    cat("Running trial", i, "\n")

    
  Ntot       <- 0
  is.MTD     <- FALSE
  stop.trial <- FALSE
  too.toxic  <- FALSE
  
  SimulationTrial <- matrix(nrow = 50, ncol = 25)
    colnames(SimulationTrial) = c("Trial", "Cohort", "Ntox", "Npat", "Weight", "DoseAdm",  
                                  "under",  "target", "over", "true.prob", "truePK", "mean(PK)", "sd(PK)", 
                                  "nextDoseAdm", "nextunder", "nexttarget", "nextover", "nextprob", 
                                  "MTD", "too.toxic", "DoseGroup",
                                  "Rhat.max", "crash", "Rseed", "WBseed")

    Agent          <- Design$Agent
    Doses          <- Design$Doses
    Dose.Grp       <- Design$Dose.Grp
    DoseRef        <- Design$DoseRef
    PK             <- Design$PK
    PK.sd          <- Design$PK.sd
    Prior.logAB    <- Design$Prior.logAB
    pmix.logAB     <- Design$pmix.logAB
    Prior.logAB.PK <- Design$Prior.logAB.PK
    pmix.logAB.PK  <- Design$pmix.logAB.PK
    dose.pk        <- cbind(Doses, PK, Dose.Grp)
    

    pars=c("logAlphaBeta","pTox","pCat")


    DosesAdm.obj <- NULL
    Dose.Grp.obj <- NULL
    Npat.obj     <- NULL
    Ntox.obj     <- NULL
    Weight.obj   <- NULL
    PK.obj       <- NULL
    
    currentDoseAdm <- StartDose
    cohort <- 1
    crashed <- FALSE
    
    while (!stop.trial & !crashed) {
      
      if (cohort > nrow(SimulationTrial)) {
        SimulationTrialL <- matrix(nrow = nrow(SimulationTrial) + 15, ncol = 23)
        colnames(SimulationTrialL) <- colnames(SimulationTrial)
        SimulationTrialL[1:nrow(SimulationTrial), ] <- SimulationTrial
        SimulationTrial <- SimulationTrialL
      }
      
      SimulationTrial[cohort, "Trial"] <- i
      SimulationTrial[cohort, "Cohort"] <- cohort
      SimulationTrial[cohort, "Rseed"] <- i

      if (length(CohortRnd[[1]]) == 1) {
        
        CohortSize <- CohortRnd[[1]]
        
      } else {
        CohortSize <- sample(CohortRnd[[1]], 1, prob = CohortRnd[[2]])
      }
      

      currentDoseGrp <- rep(Dose.Grp[which(Doses == currentDoseAdm)], CohortSize)
      current.truePK <- PK[which(Doses == currentDoseAdm)]
      currentPk.sd   <- PK.sd[which(Doses == currentDoseAdm)]
      current.tox    <- TrueProbs[which(Doses == currentDoseAdm)]
       
      # The easy way to generate toxicity data is above
      # The following is to match the first cohort data which is produced by BLRM simulation
      current.ntox <- rep(0, CohortSize)
      nntox        <- rbinom(n = 1, size = CohortSize, prob = current.tox)
      temp.index   <- sample(x = 1:CohortSize, size = nntox) 
      current.ntox[temp.index] <- 1
      
      # current.ntox   <- try(rbinom(n = CohortSize, size = 1, prob = current.tox))
      currentPK      <- rlnorm(n = CohortSize, meanlog = current.truePK, sdlog = currentPk.sd)
    
      
      
      DosesAdm.obj <- c(DosesAdm.obj, rep(currentDoseAdm, CohortSize))
      Dose.Grp.obj <- c(Dose.Grp.obj, rep(currentDoseGrp[1], CohortSize))
      Npat.obj     <- c(Npat.obj, rep(1, CohortSize))
      Ntox.obj     <- c(Ntox.obj, current.ntox)
      Weight.obj   <- c(Weight.obj, rep(1, CohortSize))
      PK.obj       <- c(PK.obj, sample(x = currentPK, size = CohortSize, replace = TRUE))
      
      
      SimulationTrial[cohort, "Ntox"]      <- sum(current.ntox)
      SimulationTrial[cohort, "Npat"]      <- CohortSize
      SimulationTrial[cohort, "Weight"]    <- 1
      SimulationTrial[cohort, "DoseAdm"]   <- currentDoseAdm
      SimulationTrial[cohort, "true.prob"] <- current.tox
      SimulationTrial[cohort, "truePK"]    <- current.truePK
      SimulationTrial[cohort, "mean(PK)"]  <- mean(currentPK)
      SimulationTrial[cohort, "sd(PK)"]    <- currentPk.sd
      SimulationTrial[cohort, "DoseGroup"] <- currentDoseGrp[1]
      

      Ntot <- Ntot + CohortSize
      
      
      if (Histdatf == 1) {
        current.data <- list(DosesAdm = c(Histdata$DosesAdm, SimulationTrial[, "DoseAdm"][1:cohort]),
                             Dose.Grp = c(Histdata$DosesAdm, SimulationTrial[, "DoseGroup"][1:cohort]),
                             Npat     = c(Histdata$Npat, SimulationTrial[, "Npat"][1:cohort]),
                             Ntox     = c(Histdata$Ntox, SimulationTrial[, "Ntox"][1:cohort]),
                             Weight   = c(Histdata$Weight, SimulationTrial[, "Weight"][1:cohort]),
                             PK       = c(Histdata$DosesAdm, SimulationTrial[, "PK"][1:cohort]),
                             X        = NULL
                             )
      } else {
        current.data <- list(DosesAdm = DosesAdm.obj,
                             Dose.Grp = Dose.Grp.obj,
                             Npat     = Npat.obj,
                             Ntox     = Ntox.obj,
                             Weight   = Npat.obj,
                             PK       = PK.obj, 
                             X        = NULL
        )
      }

      n.att <- 1
      did.run <- FALSE

      sink(paste(Outfile, "_log.txt", sep=""))
      

      while((n.att <= 10) & !did.run) {
        
        simSeeds <- sample.int(.Machine$integer.max, 2)
        currentRun <- try(SingleBLRM.PK(Data           = current.data,
                                        Agent          = Agent,
                                        Doses          = Doses,
                                        DoseRef        = DoseRef,
                                        Prior.logAB    = Prior.logAB,
                                        pmix.logAB     = c(1),
                                        monotone.logAB = TRUE,
                                        Prior.logAB.PK = Prior.logAB.PK,
                                        pmix.logAB.PK  = pmix.logAB.PK,
                                        monotone.PK    = TRUE,
                                        Prior.sd.PK    = c(0.25, log(2)/1.96),
                                        random.model   = FALSE,
                                        Prior.beta     = rbind(NULL,c(0,10)),
                                        Prior.Only     = FALSE,
                                        Pcutoffs       = c(0.16, 0.33),
                                        New.cohortsize = 3,
                                        Outfile        = "PTK7_BLRMPK_Fixed_M",
                                        MCMC           = Nsim,
                                        DIC            = TRUE,
                                        Rnd.seed       = simSeeds[1],
                                        RSeed.WinBUGS  = simSeeds[2]
        ))
                       
          
        did.run <- is.null(attr(currentRun, "class"))
        n.att <- n.att + 1
        SimulationTrial[cohort, "WBseed"] <- simSeeds[2]
        
      }
      
      sink()

      if (!did.run) {
        crashed <- TRUE
        SimulationTrial[, "crash"] <- TRUE
      } else {
        
        
        SimulationTrial[, "crash"]          <- FALSE
        SimulationTrial[cohort, "Rhat.max"] <- currentRun$Rhat$Rhat.max
        SimulationTrial[cohort, "under"]    <- currentRun$Pcat[, 1,, ][which(currentDoseAdm == Doses)]
        SimulationTrial[cohort, "target"]   <- currentRun$Pcat[, 2,, ][which(currentDoseAdm == Doses)]
        SimulationTrial[cohort, "over"]     <- currentRun$Pcat[, 3,, ][which(currentDoseAdm == Doses)]

        #print(SimulationTrial)

        decisionVars <- names(formals(Decision))
        decision <- do.call(Decision, mget(decisionVars))

        too.toxic <- decision$too.toxic
        stop.trial <- decision$stop.trial
        is.MTD <- decision$is.MTD
        nextDoseAdm <- decision$DoseAdm
        
        print(decision)
        print(currentRun$Data$data[[1]])

        if (!too.toxic) {

          SimulationTrial[cohort, "nextunder"]  <- currentRun$Pcat[, 1,, ][which(nextDoseAdm == Doses)]
          SimulationTrial[cohort, "nexttarget"] <- currentRun$Pcat[, 2,, ][which(nextDoseAdm == Doses)]
          SimulationTrial[cohort, "nextover"]   <- currentRun$Pcat[, 3,, ][which(nextDoseAdm == Doses)]
          SimulationTrial[cohort, "nextprob"]   <- TrueProbs[which(Doses == nextDoseAdm)]


          SimulationTrial[cohort, "MTD"]         <- is.MTD
          SimulationTrial[cohort, "too.toxic"]   <- too.toxic
          SimulationTrial[cohort, "nextDoseAdm"] <- nextDoseAdm


        } else {
          
          SimulationTrial[cohort, "nextunder"]   <- NA
          SimulationTrial[cohort, "nexttarget"]  <- NA
          SimulationTrial[cohort, "nextover"]    <- NA
          SimulationTrial[cohort, "MTD"]         <- is.MTD
          SimulationTrial[cohort, "too.toxic"]   <- too.toxic
          SimulationTrial[cohort, "nextprob"]    <- NA
          SimulationTrial[cohort, "nextDoseAdm"] <- nextDoseAdm
          
        }
        
        currentDoseAdm <- nextDoseAdm
        
        
        cat(paste("Cohort: ", cohort, "\n", sep = ""))
        cohort <- cohort + 1
        
      }
      
    }
    

    mtd.done <- subset(SimulationTrial, SimulationTrial[, "MTD"] == 1)
    simsummaryTrial <- mtd.done
    #print(SimulationTrial[1:(cohort - 1), ])
    ncohort <- matrix(NA, nrow = NROW(Doses), ncol = 1)
    dltcohort <- matrix(NA, nrow = NROW(Doses), ncol = 1)
    for (l in 1:NROW(Doses)) {
        ncohort1       <- subset(SimulationTrial, SimulationTrial[, "DoseAdm"] == Doses[l])
        ncohort[l,1]   <- sum(ncohort1[, "Npat"])
        dltcohort[l,1] <- sum(ncohort1[, "Ntox"])
                 }

    NMATRXTrial   <- ncohort/sum(ncohort)
    NNMATRXTrial  <- ncohort
    DLTMATRXTrial <- dltcohort
    cat(paste("Trial ", i, ": finalized.\n", sep = ""))
    
    
    res <- list(Simulation = SimulationTrial, simsummary = simsummaryTrial,
                NMATRX = NMATRXTrial, NNMATRX = NNMATRXTrial, DLTMATRX = DLTMATRXTrial)
    res
  }


    # library(doParallel)
    # nCores <- detectCores()
    # cl <- makeCluster(nCores)
    # registerDoParallel(cl)
    # # registerDoParallel(cores=4)
    # registerDoRNG(seed = 123, once = FALSE) #SEED1 : To initiate enviorment. THIS IS THE IMPORTANT STEP FOR PARALLELISM AND REPRODUCIBILITY
    # getDoParWorkers()
    ######################################################################################################
    
    
    # Results <- foreach(i = 1:nTrials, 
    #                    .packages=c("methods","R2WinBUGS","rjags","R2jags"),
    #                    .options.RNG = SeedParam
    # ) %dorng% parallel.sim(irep)
    
    
    #########################################################################################################
    #Stop all clusters (DO NOT CHANGE)
    
    # OncoPh1BLRM:::unregister()
    # #stopCluster()  #for unix
    # stopImplicitCluster()
    
    
  Doses   <- Design$Doses
  Agent   <- Design$Agent
  

  extract    <- function(struct, elem) lapply(struct, "[[", elem)
  Simulation <- extract(Results, "Simulation")
  simsummary <- extract(Results, "simsummary")
  NMATRX     <- extract(Results, "NMATRX")
  NNMATRX    <- extract(Results, "NNMATRX")
  DLTMATRX   <- extract(Results, "DLTMATRX")
  RNGSEEDS   <- attr(Results, "rng")



  MTDsim <- list()
  simsummarys <- do.call("rbind", simsummary)
  if (NROW(simsummarys) > 0) {
    for (q in 1:nrow(simsummarys)) {
      MTDmatrix <- matrix(NA, nrow = NROW(Doses), ncol = 1)
      for (a in 1:NROW(Doses)) {
          MTDmatrix[a, 1] <- (simsummarys[q, "DoseAdm"] == Doses[a])
                               }
      MTDsim[[q]] <- MTDmatrix
    }
  }

  if (NROW(simsummarys) == 0) {
    MTDmatrix <- matrix(NA, nrow = NROW(Doses), ncol = 1)

    MTDsim[[1]] <- MTDmatrix

  }

  res.sim=list(simul      = Simulation,
               simsummary = simsummarys,
               NMATRX     = NMATRX,
               NNMATRX    = NNMATRX,
               DLTMATRX   = DLTMATRX,
               MTDsim     = MTDsim,
               RNGSEEDS   = RNGSEEDS)


  sim.sum.Vars <- names(formals(SingleBLRMSimSummary))

  sim.sum <- do.call(SingleBLRMSimSummary, mget(sim.sum.Vars))

  save(list = names(res.sim), file = paste(Outfile, ".Rdata", sep = ""), envir = list2env(res.sim))

  res.sim.f <- list(res.sim=res.sim, sim.sum= sim.sum)

  return(res.sim.f)

}
