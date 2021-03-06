% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Combo2BLRMSim.R
\name{Combo2BLRMSim}
\alias{Combo2BLRMSim}
\title{Combo2BLRMSim Function}
\usage{
Combo2BLRMSim(
  Design = NULL,
  nTrials = 1000,
  TrueProbs = NULL,
  StartDose = NULL,
  CohortRnd = list(c(3, 4, 5, 6), c(0.8, 0.1, 0.05, 0.05)),
  Histdatf = 0,
  Histdata = NULL,
  Pcutoffs = c(0.16, 0.33),
  Nsim = c(2500, 5000, 2, 1),
  Decision = Combo2BLRMSimDecisionRule(MaxN = 60, mnMTD = 6, mnTOT = 15, mnTAR = 1,
    EWOC = 0.25, Escalation = "target", doseMaxEsc1 = 2, doseMaxEsc2 = 2),
  Outfile = "Sims",
  Out.summary = "sims_summary.txt",
  RndSeed = 12
)
}
\arguments{
\item{Design}{Dual combination design specification: contains required parametrs of \link[OncoPh1BLRM]{Combo2BLRM} (a list objest contains Agent1 (name of agent 1),
Agent1 (name of Agent 2), Doses1 (doses of agent 1), Doses2 (doses of Agent 2), DoseRef1
(reference dose of agent 1), DoseRef2 (reference dose of agent 2), Prior.logAB.1 (prior for agent 1),
pmix.logAB1 (mixture probability of agent 1), Prior.logAB.2 (prior for agent 2), pmix.logAB2
(mixture probability of agent 2), Prior.eta (prior for eta), pmix.eta (mixture probability of eta))}

\item{nTrials}{Number of simulated trials}

\item{TrueProbs}{True toxicity matrix for each combination dose (a matrix with dimension number of
doseS for agent 1 X number of doses for agent 2)}

\item{CohortRnd}{Cohort size (Random). A list with the first element a vector of the available
cohort size(s), and the second element a vector of the same length with the 
respective selection probability of the cohort size. For example, when set to  
list(c(3, 2), c(0.7, 0.3)), then the cohortsize is chosen as 3 or 2 at random 
with probability 0.7 or 0.3 respectively. If a fixed cohort size is desired use 
e.g. CohortRnd = list(c(3), c(1))}

\item{Histdatf}{logical: if TRUE, historical data is used (optional)}

\item{Histdata}{a list with named elements  DosesAdm1, DosesAdm2, Npat, Ntox, Weight containing the
historical data (optional)}

\item{Pcutoffs}{cutoffs to define the intervals for the true event probability. The default is
c(0.16, 0.33) defining three categories, [0.00 -0.16); [0.16-0.33); (0.33-1.00] which correspons to
underdosing, targeted dosing and overdosing}

\item{Nsim}{4-vector of MCMC parameters n.burnin, n.iterations, n.chains, n.thin (see jags)
number of chains, thin. Default is c(2500, 5000, 2, 1).}

\item{Decision}{Function which decides after each simulated cohort how to proceed}

\item{Outfile}{A R object saves all simulations output objects}

\item{Out.summary}{A txt file with simulation summary}

\item{RndSeed}{Seed used for the random number generator. Required parameter}
}
\value{
Input arguments are 
\item{Design}{The input data and prior information}
\item{nTrials}{Number of simulated trials}
\item{TrueProbs}{True toxicity vector}
\item{StartDose}{Start dose for the escalation trial}
\item{Pcutoffs}{Target toxicity interval}
\item{CohortRnd}{A list of size two with first quantity is the vector of cohortsize and second quantity is 
                 the vector of probabilities assigned to each cohortsize. For example, when set to 
                 list(c(3, 2), c(0.7, 0.3)), then the cohortsize is chosen as 3 or 2 at random 
                 with probability 0.7 or 0.3 respectively.}
\item{Decision}{}
\item{Outfile}{Name of the output file}
\item{Out.summary}{Name of the summary output file}
\item{RndSeed}{The seed which was fixed to reproduce the result}

Output arguments are 
\item{DLTMATRX}{the occurrences of DLT's per dose per trial}
\item{MTDSIM}{a binary matrix which is only TRUE when MTD has been observed in the trial}
\item{NNMATRX}{the cohort sizes per dose per trial}
\item{res.sim$simul}{descrioption of each trial -- Trial no. (Trial), Cohort no. (Cohort), 
Number of toxicity (Ntox), Number of patients (Npat), Weight, Current dose (DoseAdm), 
underdose probability (under), probability of target toxicity (target), 
overdose probability (over), probability of the current dose (true.prob), 
suggestion for the next dose, if any, (nextDoseAdm), underdose probability of the next 
dose (nextunder), target toxicity probability of the next dose (nexttarget), 
probability of the next dose (nextprob), a binary variable which is only 1 (TRUE) 
if the current dose is declared as the MTD, if the first dose is so toxic that 
\eqn{\Pr(\pi_d > 0.33) \ge EWOC} then the trial is stooped and too.toxic variable is set to 1, 
the Rhat.max information, crash = 1 if the code is crashed, the seed number for R (Rseed) 
and JAGS (WBseed)}
\item{res.sim$simsummary}{description of those cohorts within those trials 
in which the corresponding dose is the MTD}
\item{sim.sum}{simulation summary -- Proportion of patients in target dose, 
Proportion of patients in over dose, Proportion of patients in under dose,
Proportion of trials with MTD within target dose region, 
Proportion of trails with MTD within over dose region, 
Proportion of trails with MTD within under dose region, 
Proportion of trials stopped due to too toxicity at the start dose, 
Average sample size}

In addition, an R object output and a text output are generated containing all of the above.
}
\description{
This function performs trial simulation with dual combination model
}
\details{
See vignettes. browseVignettes("OncoPh1BLRM")
}
\examples{
\dontrun{
# Design Parameters
Design.ex <- list(Agent1        = "Agent1",
                  Agent2        = "Agent2",
                  Doses1        = c(3, 4.5, 6),
                  Doses2        = c(400, 600, 800),
                  DoseRef1      = 3,
                  DoseRef2      = 800,
                  Prior.logAB.1 = list(c(-2.50878483, 0.02742675, 1.00494466, 
                  0.87341511, -0.53973581)),
                  pmix.logAB1   = c(1),
                  Prior.logAB.2 = list(c(-2.1273926, 0.7083885, 0.6155286, 
                  0.6569739, -0.2044163)),
                  pmix.logAB2   = c(1),
                  Prior.eta     = list(c(0.09531018, 1.07240530)),
                  pmix.eta      = c(1)
                  )

# Simulation Scenario 1: using the mean from prior
TrueProbs.prior <- rbind(c(0.11944100, 0.15527464, 0.2090764),
                         c(0.17593587, 0.21362416,  0.2714988),
                         c(0.23071554, 0.27402758,  0.3340971))

#Running Simulation
sim.args <- list(Design = Design.ex, 
                 nTrials = 5,
                 TrueProbs = TrueProbs.prior,
                 StartDose = c(3, 400),
                 Pcutoffs  = c(0.16, 0.33),
                 Decision  = Combo2BLRMSimDecisionRule(MaxN = 60, mnMTD = 6, mnTOT = 15, 
                                                       mnTAR = 0, EWOC = 0.25, 
                                                       Escalation = "maximal"),
                 Outfile = "Sims_dual",
                 Out.summary= "Sim_dual_summary.txt",
                 RndSeed = 123
)

# Setting up back-end for parallel computing (DO NOT CHANGE)
library(doParallel)
require(doRNG)
nCores <- detectCores()  # number of cores
cl <- makeCluster(nCores) 
registerDoParallel(cl)
registerDoRNG(seed = 123, once = FALSE) #SEED1 : To initiate enviorment. 
# THIS IS THE IMPORTANT STEP FOR PARALLELISM AND REPRODUCIBILITY
getDoParWorkers()

simulation.dual <- do.call(Combo2BLRMSim, sim.args)

# Stop all clusters (DO NOT CHANGE)
OncoPh1BLRM:::unregister()
stopImplicitCluster() 
}


}
\references{
Neuenschwander B, Matano A, Tang Z, Roychoudhury S,Wandel S, Bailey S. 
A Bayesian Industry Approach to Phase I Combination Trials in Oncology. 
Statistical Methods for Drug Combination Studies. Taylor & Francis, 2015.
}
\seealso{
\code{\link{SingleBLRMSim}}, \code{\link{Combo2BLRM}}
}
