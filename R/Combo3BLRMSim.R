#options(warn=-1)

#' @title Combo3BLRMSim Function
#'
#' @description This function performs trial simulation with triple combination model
#'
#' @param Design Triple combination design specification: contains required parametrs of \link[OncoPh1BLRM]{Combo3BLRM} (a list objest contains Agent1 (name of agent 1),
#'  Agent1 (name of Agent 2), Doses1 (doses of agent 1), Doses2 (doses of Agent 2), Doses3 (doses of Agent 3),
#'  DoseRef1 (reference dose of agent 1), DoseRef2 (reference dose of agent 2), DoseRef3 (reference dose of agent 3),
#'  Prior.logAB.1 (prior for agent 1), pmix.logAB1 (mixture probability of agent 1), Prior.logAB.2 (prior for agent 2),
#'  pmix.logAB2 (mixture probability of agent 2), Prior.eta (prior for eta)). Dose of Agent 3 is fixed.
#' @param nTrials Number of simulated trials
#' @param TrueProbs True toxicity matrix for each combination dose (a matrix with dimension number of
#'  doses for agent 1 X number of doses for agent 2)
#' @param CohortRnd Cohort size (Random). A list with the first element a vector of the available
#'                  cohort size(s), and the second element a vector of the same length with the 
#'                  respective selection probability of the cohort size. For example, when set to  
#'                  list(c(3, 2), c(0.7, 0.3)), then the cohortsize is chosen as 3 or 2 at random 
#'                  with probability 0.7 or 0.3 respectively. If a fixed cohort size is desired use 
#'                  e.g. CohortRnd = list(c(3), c(1))
#' @param Histdatf logical: if TRUE, historical data is used (optional)
#' @param Histdata a list with named elements  DosesAdm1, DosesAdm2, Npat, Ntox, Weight containing the
#'  historical data (optional)
#' @param Pcutoffs cutoffs to define the intervals for the true event probability. The default is
#'  c(0.16, 0.33) defining three categories, [0.00 -0.16); [0.16-0.33); (0.33-1.00] which correspons to
#'  underdosing, targeted dosing and overdosing
#' @param Nsim 4-vector of MCMC parameters n.burnin, n.iterations, n.chains, n.thin (see jags)
#'  number of chains, thin. Default is c(2500, 5000, 2, 1).
#' @param Decision Function which decides after each simulated cohort how to proceed
#' @param Outfile A R object saves all simulations output objects
#' @param Out.summary A txt file with simulation summary
#' @param RndSeed Seed used for the random number generator. Required parameter
#'
#'
#' @details See vignettes. browseVignettes("OncoPh1BLRM") 
#' 
#' Note: Currently only a fix dose is allowed for the third agent. 
#' 
#' @references Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. A Bayesian Industry
#' Approach to Phase I Combination Trials in Oncology. Statistical Methods for Drug Combination
#' Studies. Taylor & Francis, 2015.
#' 
#' 
#' @return Input arguments are 
#' \item{Design}{The input data and prior information}
#' \item{nTrials}{Number of simulated trials}
#' \item{TrueProbs}{True toxicity vector}
#' \item{StartDose}{Start dose for the escalation trial}
#' \item{Pcutoffs}{Target toxicity interval}
#' \item{CohortRnd}{A list of size two with first quantity is the vector of cohortsize and second quantity is 
#'                  the vector of probabilities assigned to each cohortsize. For example, when set to 
#'                  list(c(3, 2), c(0.7, 0.3)), then the cohortsize is chosen as 3 or 2 at random 
#'                  with probability 0.7 or 0.3 respectively.}
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
#'
#' @seealso \code{\link{SingleBLRMSim}}, \code{\link{Combo3BLRM}}
#' 
#' @examples 
#' \dontrun{
#' # Design Parameters
#' Design.ex <- list(Agent1        = "Agent1",
#'                   Agent2        = "Agent2",
#'                   Agent3        = "Agent3",
#'                   Doses1        = c(0.5, 0.75, 1.03),
#'                   Doses2        = c(22.5, 30, 45),
#'                   Doses3        = c(10),
#'                   DoseRef1      = 1,
#'                   DoseRef2      = 45,
#'                   DoseRef3      = 10,
#'                   Prior.logAB.1 = list(c(-1.5185331, 0.2237979, 0.6732011, 0.7611914, 0.2236197)),
#'                   pmix.logAB1   = c(1),
#'                   Prior.logAB.2 = list(c(-3.6214488, 0.4964869, 0.7296964, 0.9112984, -0.3744258)),
#'                   pmix.logAB2   = c(1),
#'                   Prior.logAB.3 = list(c(-2.6788781, -0.0429276, 0.9714378, 0.8205235, -0.2379236)),
#'                   pmix.logAB3   = c(1),
#'                   Prior.eta     = list(etaPrior12  = c(0.00000,    0.20687),
#'                                        etaPrior13  = c(0.00000,    0.20687),
#'                                        etaPrior23  = c(0.0000000,  0.6391648),
#'                                        etaPrior123 = c(0.00000000, 0.02489294)
#'                  )
#' )
#'
#'
#' # Simulation Scenario
#' TrueProbs.prior <- rbind(c(0.1877346, 0.1949610, 0.2157443),
#'                          c(0.2328766, 0.2399549, 0.2598694),
#'                          c(0.2877317, 0.2942871, 0.3123852)
#'                          )
#'
#' # Running Simulation
#' sim.args <- list(Design      = Design.ex,
#'                  nTrials     = 5,
#'                  TrueProbs   = TrueProbs.prior,
#'                  StartDose   = c(0.5,30,10),
#'                  Pcutoffs    = c(0.16, 0.33),
#'                  Decision    = Combo3BLRMSimDecisionRule(MaxN = 60, mnMTD = 6, mnTOT = 15, mnTAR = 0, 
#'                                                          EWOC = 0.25, Escalation = "maximal"),
#'                  Outfile     = "Sims_triple",
#'                  Out.summary = "Sim_triple_summary.txt",
#'                  RndSeed     = 123
#' )
#'
#' # Setting up back-end for parallel computing (DO NOT CHANGE)
#' library(doParallel)
#' require(doRNG)
#' nCores <- detectCores()  # number of cores
#' cl <- makeCluster(nCores) 
#' registerDoParallel(cl)
#' registerDoRNG(seed = 123, once = FALSE) #SEED1 : To initiate enviorment. 
#' # THIS IS THE IMPORTANT STEP FOR PARALLELISM AND REPRODUCIBILITY
#' getDoParWorkers()
#' 
#' simulation.triple <- do.call(Combo3BLRMSim, sim.args)
#'
#' # Stop all clusters (DO NOT CHANGE)
#' OncoPh1BLRM:::unregister()
#' stopImplicitCluster() 
#' }
#' 
#' 
#' @export


Combo3BLRMSim <- function (Design =NULL,
                           nTrials = 1000,
                           TrueProbs= NULL,
                           StartDose = NULL,
                           CohortRnd=  list(c(  3,   4,    5,    6),
                                            c(0.8, 0.1, 0.05, 0.05)),
                           Histdatf= 0,
                           Histdata = NULL,
                           Pcutoffs= c(0.16, 0.33),
                           Nsim= c(2500, 5000, 2, 1),
                           Decision= Combo3BLRMSimDecisionRule(MaxN = 60, mnMTD = 6, mnTOT = 15, mnTAR = 1, EWOC=0.25, Escalation = "target", doseMaxEsc1 = 2, doseMaxEsc2 = 2),
                           #Stop.rule = list(MaxN = 60, mnMTD = 6, mnTOT = 15, mnTAR = 1, Escalation = "maximal", doseMaxEsc1 = 2, doseMaxEsc2 = 2),
                           Outfile = "Sims",
                           Out.summary="sims_summary.txt",
                           RndSeed=12
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

#set.seed(RndSeed)

Results <- do.call(foreach, forArgs) %dorng% {
    cat("Running trial", i, "\n")

    Ntot <- 0
    is.MTD <- FALSE
    stop.trial <- FALSE
    too.toxic <- FALSE
    SimulationTrial <- matrix(nrow = 15, ncol = 25)
    colnames(SimulationTrial) = c("Trial", "Cohort", "Ntox",
                                  "Npat", "Weight", "DoseAdm1", "DoseAdm2", "DoseAdm3","under",
                                  "target", "over", "true.prob", "nextDoseAdm1", "nextDoseAdm2",
                                  "nextDoseAdm3", "nextunder", "nexttarget", "nextover", "nextprob",
                                  "MTD", "too.toxic", "Rhat.max", "crash", "Rseed",
                                  "WBseed")


    Agent1 = Design$Agent1
    Agent2 = Design$Agent2
    Agent3 = Design$Agent3
    Doses1 = Design$Doses1
    Doses2 = Design$Doses2
    Doses3 = Design$Doses3
    DoseRef1 = Design$DoseRef1
    DoseRef2 = Design$DoseRef2
    DoseRef3 = Design$DoseRef3

    Prior.logAB.1= Design$Prior.logAB.1
    pmix.logAB1= Design$pmix.logAB1
    Prior.logAB.2= Design$Prior.logAB.2
    pmix.logAB2= Design$pmix.logAB2
    Prior.logAB.3= Design$Prior.logAB.3
    pmix.logAB3= Design$pmix.logAB3

    Prior.eta= Design$Prior.eta
    pars=c("logAlphaBeta1", "logAlphaBeta2", "logAlphaBeta3", "eta12", "eta13","eta23", "eta123", "P12","pCat")


    StartDose1 <- StartDose[1]
    StartDose2 <- StartDose[2]
    StartDose3 <- StartDose[3]

    currentDoseAdm1 <- StartDose1
    currentDoseAdm2 <- StartDose2
    currentDoseAdm3 <- StartDose3
    cohort <- 1
    crashed <- FALSE
    while (!stop.trial & !crashed) {
      if (cohort > nrow(SimulationTrial)) {
        SimulationTrialL <- matrix(nrow = nrow(SimulationTrial) + 15, ncol = 25)
        colnames(SimulationTrialL) <- colnames(SimulationTrial)
        SimulationTrialL[1:nrow(SimulationTrial), ] <- SimulationTrial
        SimulationTrial <- SimulationTrialL
      }
      SimulationTrial[cohort, "Trial"] <- i
      SimulationTrial[cohort, "Cohort"] <- cohort
      SimulationTrial[cohort, "Rseed"] <- i

      if (length(CohortRnd[[1]]) == 1) {
        CohortSize <- CohortRnd[[1]]
      }
      else {
        CohortSize <- sample(CohortRnd[[1]], 1, prob = CohortRnd[[2]])
      }
      SimulationTrial[cohort, "Npat"] <- CohortSize
      Ntot <- Ntot + CohortSize
      SimulationTrial[cohort, "DoseAdm1"] <- currentDoseAdm1
      SimulationTrial[cohort, "DoseAdm2"] <- currentDoseAdm2
      SimulationTrial[cohort, "DoseAdm3"] <- currentDoseAdm3

      if (NROW(Doses1) > 1 & NROW(Doses2) > 1) {
        SimulationTrial[cohort, "Ntox"] <- rbinom(1, CohortSize, TrueProbs[which(Doses1 == currentDoseAdm1), which(Doses2 == currentDoseAdm2)])
        SimulationTrial[cohort, "true.prob"] <- TrueProbs[which(Doses1 == currentDoseAdm1), which(Doses2 == currentDoseAdm2)]
      }
      if (NROW(Doses1) == 1) {
        SimulationTrial[cohort, "Ntox"] <- rbinom(1, CohortSize, TrueProbs[which(Doses2 == currentDoseAdm2)])
        SimulationTrial[cohort, "true.prob"] <- TrueProbs[which(Doses2 == currentDoseAdm2)]
      }
      if (NROW(Doses2) == 1) {
        SimulationTrial[cohort, "Ntox"] <- rbinom(1, CohortSize, TrueProbs[which(Doses1 == currentDoseAdm1)])
        SimulationTrial[cohort, "true.prob"] <- TrueProbs[which(Doses1 == currentDoseAdm1)]
      }
      SimulationTrial[cohort, "Weight"] <- 1
      if (Histdatf == 1) {
        current.data <- list(DosesAdm1 = c(Histdata$DosesAdm1, SimulationTrial[, "DoseAdm1"][1:cohort]),
                             DosesAdm2 = c(Histdata$DosesAdm2, SimulationTrial[, "DoseAdm2"][1:cohort]),
                             DosesAdm3 = c(Histdata$DosesAdm3, SimulationTrial[, "DoseAdm3"][1:cohort]),
                             Npat = c(Histdata$Npat, SimulationTrial[, "Npat"][1:cohort]),
                             Ntox = c(Histdata$Ntox, SimulationTrial[, "Ntox"][1:cohort]),
                             Weight = c(Histdata$Weight, SimulationTrial[, "Weight"][1:cohort])
                             )
                       }
      else {
        current.data <- list(DosesAdm1 = SimulationTrial[,"DoseAdm1"][1:cohort],
                             DosesAdm2 = SimulationTrial[,"DoseAdm2"][1:cohort],
                             DosesAdm3 = SimulationTrial[,"DoseAdm3"][1:cohort],
                             Npat = SimulationTrial[, "Npat"][1:cohort],
                             Ntox = SimulationTrial[, "Ntox"][1:cohort],
                             Weight = SimulationTrial[, "Weight"][1:cohort]
                             )
      }

      n.att <- 1
      did.run <- FALSE

sink(paste(Outfile, "_log.txt", sep=""))

      while ((n.att <= 10) & !did.run) {
        simSeeds <- sample.int(.Machine$integer.max, 2)
        currentRun <- try(Combo3BLRM(Data          = current.data,
                                     Agent1        = Agent1,
                                     Doses1        = Doses1,
                                     DoseRef1      = DoseRef1,
                                     Prior.logAB.1 = Prior.logAB.1,
                                     pmix.logAB1   = pmix.logAB1,
                                     Agent2        = Agent2,
                                     Doses2        = Doses2,
                                     DoseRef2      = DoseRef2,
                                     Prior.logAB.2 = Prior.logAB.2,
                                     pmix.logAB2   = pmix.logAB2,
                                     Agent3        = Agent3,
                                     Doses3        = Doses3,
                                     DoseRef3      = DoseRef3,
                                     Prior.logAB.3 = Prior.logAB.3,
                                     pmix.logAB3   = pmix.logAB3,
                                     Prior.eta     = Prior.eta,
                                     Pcutoffs      = Pcutoffs,
                                     MCMC          = Nsim,
                                     Warn          = F,
                                     pars          = pars,
                                     Prior.Only    = F,
                                     Plot          = F,
                                     progress.bar  = "none",
                                     Rnd.seed      = simSeeds[1],
                                     RSeed.WinBUGS = simSeeds[2]
                                     ))
        did.run <- is.null(attr(currentRun, "class"))
        n.att <- n.att + 1
        SimulationTrial[cohort, "WBseed"] <- simSeeds[2]
      }
sink()

      if (!did.run) {
        crashed <- TRUE
        SimulationTrial[, "crash"] <- TRUE
      }
      else {
        SimulationTrial[, "crash"] <- FALSE
        SimulationTrial[cohort, "Rhat.max"] <- currentRun$Rhat$Rhat.max

        if (NROW(Doses1) > 1 & NROW(Doses2) > 1) {
          SimulationTrial[cohort, "under"] <- currentRun$Pcat[,1,,,][which(currentDoseAdm1 == Doses1), which(currentDoseAdm2 == Doses2)]
          SimulationTrial[cohort, "target"] <- currentRun$Pcat[,2,,,][which(currentDoseAdm1 == Doses1), which(currentDoseAdm2 == Doses2)]
          SimulationTrial[cohort, "over"] <- currentRun$Pcat[,3,,,][which(currentDoseAdm1 == Doses1), which(currentDoseAdm2 == Doses2)]
          }
        if (NROW(Doses1) == 1) {
          SimulationTrial[cohort, "under"] <- currentRun$Pcat[,1,,,][which(currentDoseAdm2 == Doses2)]
          SimulationTrial[cohort, "target"] <- currentRun$Pcat[,2,,,][which(currentDoseAdm2 == Doses2)]
          SimulationTrial[cohort, "over"] <- currentRun$Pcat[,3,,,][which(currentDoseAdm2 == Doses2)]
        }
        if (NROW(Doses2) == 1) {
          SimulationTrial[cohort, "under"] <- currentRun$Pcat[,1,,,][which(currentDoseAdm1 == Doses1)]
          SimulationTrial[cohort, "target"] <- currentRun$Pcat[,2,,,][which(currentDoseAdm1 == Doses1)]
          SimulationTrial[cohort, "over"] <- currentRun$Pcat[,3,,,][which(currentDoseAdm1 == Doses1)]
        }


        decisionVars <- names(formals(Decision))
        decision <- do.call(Decision, mget(decisionVars))


        too.toxic <- decision$too.toxic
        stop.trial <- decision$stop.trial
        is.MTD <- decision$is.MTD
        nextDoseAdm1 <- decision$DoseAdm1
        nextDoseAdm2 <- decision$DoseAdm2
        nextDoseAdm3 <- decision$DoseAdm3


        if (!too.toxic) {
          if (NROW(Doses1) > 1 & NROW(Doses2) > 1) {
            SimulationTrial[cohort, "nextunder"] <- currentRun$Pcat[,1,,,][which(nextDoseAdm1 == Doses1), which(nextDoseAdm2 == Doses2)]
            SimulationTrial[cohort, "nexttarget"] <- currentRun$Pcat[,2,,,][which(nextDoseAdm1 == Doses1), which(nextDoseAdm2 == Doses2)]
            SimulationTrial[cohort, "nextover"] <- currentRun$Pcat[,3,,,][which(nextDoseAdm1 == Doses1), which(nextDoseAdm2 == Doses2)]
            SimulationTrial[cohort, "nextprob"] <- TrueProbs[which(Doses1 == nextDoseAdm1), which(Doses2 == nextDoseAdm2)]
          }
          if (NROW(Doses1) == 1) {
            SimulationTrial[cohort, "nextunder"] <- currentRun$Pcat[, 1,,,][which(nextDoseAdm2 == Doses2)]
            SimulationTrial[cohort, "nexttarget"] <- currentRun$Pcat[,2,,,][which(nextDoseAdm2 == Doses2)]
            SimulationTrial[cohort, "nextover"] <- currentRun$Pcat[,3,,,][which(nextDoseAdm2 == Doses2)]
            SimulationTrial[cohort, "nextprob"] <- TrueProbs[which(Doses2 == nextDoseAdm2)]
          }
          if (NROW(Doses2) == 1) {
            SimulationTrial[cohort, "nextunder"] <- currentRun$Pcat[,1,,,][which(nextDoseAdm1 == Doses1)]
            SimulationTrial[cohort, "nexttarget"] <- currentRun$Pcat[,2,,,][which(nextDoseAdm1 == Doses1)]
            SimulationTrial[cohort, "nextover"] <- currentRun$Pcat[,3,,,][which(nextDoseAdm1 == Doses1)]
            SimulationTrial[cohort, "nextprob"] <- TrueProbs[which(Doses1 == nextDoseAdm1)]
          }
          SimulationTrial[cohort, "MTD"] <- is.MTD
          SimulationTrial[cohort, "too.toxic"] <- too.toxic
          SimulationTrial[cohort, "nextDoseAdm1"] <- nextDoseAdm1
          SimulationTrial[cohort, "nextDoseAdm2"] <- nextDoseAdm2
          SimulationTrial[cohort, "nextDoseAdm3"] <- nextDoseAdm3
        }
        else {
          SimulationTrial[cohort, "nextunder"] <- NA
          SimulationTrial[cohort, "nexttarget"] <- NA
          SimulationTrial[cohort, "nextover"] <- NA
          SimulationTrial[cohort, "MTD"] <- is.MTD
          SimulationTrial[cohort, "too.toxic"] <- too.toxic
          SimulationTrial[cohort, "nextprob"] <- NA
          SimulationTrial[cohort, "nextDoseAdm1"] <- nextDoseAdm1
          SimulationTrial[cohort, "nextDoseAdm2"] <- nextDoseAdm2
          SimulationTrial[cohort, "nextDoseAdm3"] <- nextDoseAdm3
        }
        currentDoseAdm1 <- nextDoseAdm1
        currentDoseAdm2 <- nextDoseAdm2
        currentDoseAdm3 <- nextDoseAdm3
        cat(paste("Cohort: ", cohort, "\n", sep = ""))
        cohort <- cohort + 1
      }
    }

    mtd.done <- subset(SimulationTrial, SimulationTrial[, "MTD"] == 1)
    simsummaryTrial <- mtd.done
    #print(SimulationTrial[1:(cohort - 1), ])
    ncohort <- matrix(NA, nrow = NROW(Doses1), ncol = NROW(Doses2))
    dltcohort <- matrix(NA, nrow = NROW(Doses1), ncol = NROW(Doses2))
    for (l in 1:NROW(Doses1)) {
      for (h in 1:NROW(Doses2)) {
        ncohort1 <- subset(SimulationTrial, SimulationTrial[, "DoseAdm1"] == Doses1[l] & SimulationTrial[,"DoseAdm2"] == Doses2[h])
        ncohort[l, h] <- sum(ncohort1[, "Npat"])
        dltcohort[l, h] <- sum(ncohort1[, "Ntox"])
      }
    }
    NMATRXTrial <- ncohort/sum(ncohort)
    NNMATRXTrial <- ncohort
    DLTMATRXTrial <- dltcohort
    cat(paste("Trial ", i, ": finalized.\n", sep = ""))
    res <- list(Simulation = SimulationTrial, simsummary = simsummaryTrial,
                NMATRX = NMATRXTrial, NNMATRX = NNMATRXTrial, DLTMATRX = DLTMATRXTrial)
    res
}


   Agent1 = Design$Agent1
   Agent2 = Design$Agent2
   Agent3 = Design$Agent3
   Doses1 = Design$Doses1
   Doses2 = Design$Doses2
   Doses3 = Design$Doses3

  extract <- function(struct, elem) lapply(struct, "[[", elem)
  Simulation <- extract(Results, "Simulation")
  simsummary <- extract(Results, "simsummary")
  NMATRX <- extract(Results, "NMATRX")
  NNMATRX <- extract(Results, "NNMATRX")
  DLTMATRX <- extract(Results, "DLTMATRX")
  RNGSEEDS <- attr(Results, "rng")
  MTDsim <- list()
  simsummarys <- do.call("rbind", simsummary)
  if (NROW(simsummarys) > 0) {
    for (q in 1:nrow(simsummarys)) {
      MTDmatrix <- matrix(NA, nrow = NROW(Doses1), ncol = NROW(Doses2))
      for (a in 1:NROW(Doses1)) {
        for (b in 1:NROW(Doses2)) {
          MTDmatrix[a, b] <- (simsummarys[q, "DoseAdm1"] == Doses1[a] & simsummarys[q, "DoseAdm2"] == Doses2[b])
        }
      }
      MTDsim[[q]] <- MTDmatrix
    }
  }

  if (NROW(simsummarys) == 0) {
    MTDmatrix <- matrix(NA, nrow = NROW(Doses1), ncol = NROW(Doses2))

    MTDsim[[1]] <- MTDmatrix

  }
    res.sim=list(simul = Simulation,
              simsummary = simsummarys,
              NMATRX = NMATRX,
              NNMATRX = NNMATRX,
              DLTMATRX = DLTMATRX,
              MTDsim = MTDsim,
              RNGSEEDS = RNGSEEDS)


    sim.sum.Vars <- names(formals(Combo3BLRMSimSummary))

    sim.sum <- do.call(Combo3BLRMSimSummary, mget(sim.sum.Vars))

  save(list = names(res.sim), file = paste(Outfile, ".Rdata", sep = ""), envir = list2env(res.sim))

  res.sim.f <- list(res.sim=res.sim, sim.sum= sim.sum)

  return(res.sim.f)

}
