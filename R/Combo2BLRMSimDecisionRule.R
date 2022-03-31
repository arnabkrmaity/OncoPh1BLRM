#' @title Combo2BLRMSimDecisionRule Function
#'
#' @description This function performs escalation and decision rule trial simulation with dual 
#' combination model
#' 
#' @param MaxN Maximum number of patients in each trial. The default is 60. See Details
#' @param mnMTD At least these many patients have been treated at MTD. The default is 6. See Details
#' @param mnTOT A minimum number of these many patients should have already treated in the trial. 
#' The default is 15. See Details
#' @param mnTAR Probability thresold of target toxicity. The default is 0.5. See details
#' @param EWOC Escalation with overdose control criterion. The default is 0.25. See details
#' @param Escalation Escalation The escalation method to be used. Possible methods -- "maximal":
#' recomended, maxiumum dose which satisfies EWOC criteria; "target": dose which has 
#' highest target toxicity probability; "dose_max": maximum dose, doses will not be skipped  
#' @param doseMaxEsc1 Either a scaler or a vector of length of number of provisional doses for the first agent. 
#' If a scaler \eqn{k} then it used assumed that for all doses escaltion upto \eqn{k}-fold is allowed. If a vector
#' then for each dose the corresponding fold is allowed to escalate upto that dose. The first element is set to 2.
#' @param doseMaxEsc2 Either a scaler or a vector of length of number of provisional doses for the second agent. 
#' If a scaler \eqn{k} then it used assumed that for all doses escaltion upto \eqn{k}-fold is allowed. If a vector
#' then for each dose the corresponding fold is allowed to escalate upto that dose. The first element is set to 2.
#' 
#' @details At any trial dose escalation is achieved satisfying the escalation with overdose control 
#' (EWOC) principle. This is implemented by setting the target toxicity interval to Pcutoffs, with
#' an overdose control criterion \eqn{P(\pi_d > Pcutoffs[2]) < EWOC}. 
#' 
#' Dose escalations continue until declaration of the maximum tolerated dose (MTD) \eqn{d = (d_1, d_2)}. 
#' This dose must meet the following conditions:
#' 
#' 1. At least mnMTD patients have been treated at dose \eqn{d}.
#' 
#' 2. This dose satisfies one of the following conditions:
#' 
#' i. The probability of targeted toxicity at \eqn{d} exceeds mnTAR.
#' 
#' ii. A minimum number of mnTOT pateints should have already been treated in the trial. 
#' 
#' 
#' @seealso \code{\link{Combo2BLRMSim}}
#' 
#' 
#' @seealso \code{\link{SingleBLRMSimDecisionRule}}, \code{\link{Combo2BLRMSim}}
#' 
#'
#' @export



Combo2BLRMSimDecisionRule <- function (MaxN,
                                       mnMTD,
                                       mnTOT,
                                       mnTAR,
                                       EWOC,
                                       Escalation = c("dose_max", "maximal", "target"),
                                       doseMaxEsc1 = 2,
                                       doseMaxEsc2 = 2
                                       )
{

  Escalation <- match.arg(Escalation)
  if (any(doseMaxEsc1 <= 1) | any(doseMaxEsc2 <= 1))
    warning("Some dose escalation factors are smaller than or equal to 1.")
  function(currentRun, Ntot, currentDoseAdm1, currentDoseAdm2,Doses1, Doses2){

    CDoseAdm1 <- currentDoseAdm1
    CDoseAdm2 <- currentDoseAdm2
    datsummary <- currentRun$Data$data
    names(datsummary) <- toupper(names(datsummary))



    Pall <- currentRun$Pcat
    xrng <- max(which(Doses1 >= CDoseAdm1))
    yrng <- max(which(Doses2 >= CDoseAdm2))
    totNPAT <- sum(datsummary[, "NPAT"][datsummary[, "DOSEADM1"] ==  CDoseAdm1 & datsummary[, "DOSEADM2"] == CDoseAdm2 & datsummary[, "WEIGHT"] == 1])
    totNTOX <- sum(datsummary[, "NTOX"][datsummary[, "DOSEADM1"] ==  CDoseAdm1 & datsummary[, "DOSEADM2"] == CDoseAdm2 & datsummary[, "WEIGHT"] == 1])


    odprob <- Pall[, 3, ,]
    tgprob <- Pall[, 2, ,]
    dose.eligible <- odprob < EWOC #0.25
    lDeltaDose1 <- log(if (length(doseMaxEsc1) == 1) rep(doseMaxEsc1, length(Doses1)) else doseMaxEsc1)
    lDeltaDose2 <- log(if (length(doseMaxEsc2) == 1) rep(doseMaxEsc2, length(Doses2)) else doseMaxEsc2)
    ldm1 <- log(Doses1/CDoseAdm1) <= lDeltaDose1
    ldm2 <- log(Doses2/CDoseAdm2) <= lDeltaDose2
    possible.doses <- outer(ldm1, ldm2) == 1

    DM1 <- matrix(rep(Doses1, length(Doses2)), ncol = length(Doses2))
    DM2 <- t(matrix(rep(Doses2, length(Doses1)), ncol = length(Doses1)))
    ndset <- (dose.eligible + possible.doses) == 2
    DM1.H <- DM1 >= CDoseAdm1
    DM2.H <- DM2 >= CDoseAdm2
    ndsetup <- (DM1.H + DM2.H + dose.eligible + possible.doses) == 4



    if (sum(ndset) == 0) {
      too.toxic <- TRUE
      is.MTD <- FALSE
      DoseAdm1 <- NA
      DoseAdm2 <- NA
    }
    else {
      too.toxic <- FALSE
      if (sum(ndsetup) > 0) {
        if (Escalation == "dose_max") {
          dose1_cur <- which(CDoseAdm1 == Doses1)
          dose2_cur <- which(CDoseAdm2 == Doses2)
          while (all(!ndset[dose1_cur, ])) {
            dose1_cur <- dose1_cur - 1
          }
          while (dose2_cur != 1 && !ndset[dose1_cur, dose2_cur]) {
            dose2_cur <- dose2_cur - 1
          }
          dose1_next <- dose1_cur
          while (dose1_next != length(Doses1) && ndset[dose1_next + 1, dose2_cur]) {
            dose1_next <- dose1_next + 1
          }
          dose2_next <- dose2_cur
          while (dose2_next != length(Doses2) && ndset[dose1_next, dose2_next + 1]) {
            dose2_next <- dose2_next + 1
          }
          dose.pair <- matrix(c(dose1_next, dose2_next), 1, 2)
          if (NROW(Doses1) == 1)
            dose.pair <- dose2_next
          if (NROW(Doses2) == 1)
            dose.pair <- dose1_next
        }
        if (Escalation == "maximal") {
          odprob1 <- odprob
          odprob1[!ndset] <- NA
          dose.pair <- which(odprob1 == max(odprob1, na.rm = T), arr.ind = T)
        }
        if (Escalation == "target") {
          tgprob1 <- tgprob
          tgprob1[!ndsetup] <- NA
          dose.pair <- which(tgprob1 == max(tgprob1, na.rm = T), arr.ind = T)
        }
        if (NROW(Doses1) > 1 & NROW(Doses2) > 1) {
          nextDose1 <- dose.pair[1, 1]
          nextDose2 <- dose.pair[1, 2]
          next.equal <- (Doses1[nextDose1] == CDoseAdm1) & (Doses2[nextDose2] == CDoseAdm2)
        }
        if (NROW(Doses1) == 1) {
          nextDose1 <- 1
          nextDose2 <- dose.pair[1]
          next.equal <- (Doses1[nextDose1] == CDoseAdm1) & (Doses2[nextDose2] == CDoseAdm2)
        }
        if (NROW(Doses2) == 1) {
          nextDose1 <- dose.pair[1]
          nextDose2 <- 1
          next.equal <- (Doses1[nextDose1] == CDoseAdm1) & (Doses2[nextDose2] == CDoseAdm2)
        }
      }
      if (sum(ndsetup) == 0) {
        tgprob[ndset == FALSE] <- NA
        dose.pair <- which(tgprob == max(tgprob, na.rm = T), arr.ind = T)

        if (NROW(Doses1) > 1 & NROW(Doses2) > 1) {
          nextDose1 <- dose.pair[1, 1]
          nextDose2 <- dose.pair[1, 2]
          next.equal <- (Doses1[nextDose1] == CDoseAdm1) & (Doses2[nextDose2] == CDoseAdm2)
        }
        if (NROW(Doses1) == 1) {
          nextDose1 <- 1
          nextDose2 <- dose.pair[1]
          next.equal <- (Doses1[nextDose1] == CDoseAdm1) & (Doses2[nextDose2] == CDoseAdm2)
        }
        if (NROW(Doses2) == 1) {
          nextDose1 <- dose.pair[1]
          nextDose2 <- 1
          next.equal <- (Doses1[nextDose1] == CDoseAdm1) & (Doses2[nextDose2] == CDoseAdm2)
        }
      }
      boundD1 <- CDoseAdm1 == max(Doses1)
      boundD2 <- CDoseAdm2 == max(Doses2)
      mnMTD.check <- (totNPAT >= mnMTD)
      if (NROW(Doses1) > 1 & NROW(Doses2) > 1) {
        if (!boundD1 & !boundD2) {
          current.eligible <- odprob[which(Doses1 ==  CDoseAdm1), which(Doses2 == CDoseAdm2)] < EWOC
          current.highestd <- min(odprob[(which(Doses1 == CDoseAdm1) + 1):xrng, (which(Doses2 == CDoseAdm2) + 1):yrng]) >= EWOC
          current.highestx <- min(odprob[(which(Doses1 == CDoseAdm1) + 1):xrng, which(Doses2 == CDoseAdm2)]) >= EWOC
          current.highesty <- min(odprob[which(Doses1 == CDoseAdm1), (which(Doses2 == CDoseAdm2) + 1):yrng]) >= EWOC
          current.highest <- all(current.highestd, current.highestx, current.highesty)
        }
        if (boundD2 & !boundD1) {
          current.eligible <- odprob[which(Doses1 == CDoseAdm1), which(Doses2 == CDoseAdm2)] < EWOC
          current.highestx <- min(odprob[(which(Doses1 == CDoseAdm1) + 1):xrng, which(Doses2 == CDoseAdm2)]) >= EWOC
          current.highest <- current.highestx
        }
        if (boundD1 & !boundD2) {
          current.eligible <- odprob[which(Doses1 == CDoseAdm1), which(Doses2 == CDoseAdm2)] < EWOC
          current.highesty <- min(odprob[which(Doses1 == CDoseAdm1), (which(Doses2 == CDoseAdm2) + 1):yrng]) >= EWOC
          current.highest <- current.highesty
        }
        if (boundD1 & boundD2) {
          current.eligible <- odprob[which(Doses1 == CDoseAdm1), which(Doses2 == CDoseAdm2)] < EWOC
          current.highest <- odprob[which(Doses1 == CDoseAdm1), which(Doses2 == CDoseAdm2)] > EWOC
        }
        MTD.rule <- (Ntot >= mnTOT || tgprob[which(Doses1 == CDoseAdm1), which(Doses2 == CDoseAdm2)] >= mnTAR)
      }
      if (NROW(Doses1) == 1) {
        if (!boundD2) {
          current.eligible <- odprob[which(Doses2 == CDoseAdm2)] < EWOC
          current.highesty <- min(odprob[(which(Doses2 == CDoseAdm2) + 1):yrng]) >= EWOC
          current.highest <- current.highesty
        }
        if (boundD2) {
          current.eligible <- odprob[which(Doses2 == CDoseAdm2)] < EWOC
          current.highest <- odprob[which(Doses2 == CDoseAdm2)] > EWOC
        }
        MTD.rule <- (Ntot >= mnTOT || tgprob[which(Doses2 == CDoseAdm2)] >= mnTAR)
      }
      if (NROW(Doses2) == 1) {
        if (!boundD1) {
          current.eligible <- odprob[which(Doses1 == CDoseAdm1)] < EWOC
          current.highestx <- min(odprob[(which(Doses1 == CDoseAdm1) + 1):xrng]) >= EWOC
          current.highest <- current.highestx
        }
        if (boundD1) {
          current.eligible <- odprob[which(Doses1 == CDoseAdm1)] < EWOC
          current.highest <- odprob[which(Doses1 == CDoseAdm1)] > EWOC
        }
        MTD.rule <- (Ntot >= mnTOT || tgprob[which(Doses1 == CDoseAdm1)] >= mnTAR)
      }


      if (mnMTD.check & current.eligible & MTD.rule & (current.highest | next.equal)) {
        is.MTD <- TRUE
        DoseAdm1 <- CDoseAdm1
        DoseAdm2 <- CDoseAdm2
      }
      else {
        is.MTD <- FALSE
        DoseAdm1 <- Doses1[nextDose1]
        DoseAdm2 <- Doses2[nextDose2]
      }
    }
    stop.trial <- too.toxic | is.MTD | Ntot >= MaxN

    return(list(too.toxic = too.toxic, stop.trial = stop.trial,
                is.MTD = is.MTD, DoseAdm1 = DoseAdm1, DoseAdm2 = DoseAdm2))
  }
}
