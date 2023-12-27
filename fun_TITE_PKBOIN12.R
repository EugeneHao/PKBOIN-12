fun_TITE_PKBOIN12 <-
  function(index, rlist, plist, qlist, trueOBD, pT, qE, pqcorr,
           lambda_e, lambda_d, zeta1, CV, g_P, csize, cN, decisionM,
           ub, u11 = 60, u00 = 40, current = 1, doselimit = Inf,
           accrual, susp, tox_win, eff_win, tox_dist, eff_dist, tox_dist_hyper, eff_dist_hyper,
           use_susp = TRUE, accrual_random = FALSE, considerPK = TRUE
)
{

  dN <- length(plist)
  idlist <- 1:dN
  time_current <- 0                # the current time for assigning the next cohort
  time_next <- NA
  set.seed(index)

  # patDT initialization (update at the beginning of each cohort assignment)
  patDT <- data.frame(cid = rep(1:cN, each = csize), pid = rep(1:csize, cN), id = 1:(csize * cN),
                      PK = NA, tox = NA, eff = NA,
                      enroll = NA, toxend = NA, effend = NA, toxobs = -1, effobs = -1,
                      toxdat = 9999, effdat = 9999, tox_confirm = NA, eff_confirm = NA, d = NA)

  # doseDT initialization (update at the beginning of each cohort assignment)
  doseDT <- data.frame(id = 1:dN, PK = rlist, tox = plist, eff = qlist, n = rep(0, dN),
                       x = rep(0,dN), y = rep(0, dN), keep = rep(1, dN),
                       ESS_t = -1, ESS_e = -1, pi_t_hat = pT/2, pi_e_hat = qE, x_d = 0, # pT/2 and qE are default initial values in the paper
                       r_d = 0, r_sd = 0)

  earlystop <- FALSE
  record <- rep(-1, cN)          # record dose selection

  for(i in 1:cN)
  {
    # update patient information
    TITE_one <- fun_TITE_PK_update(cid = i, CV, g_P, patDT, doseDT, current, time_current, tox_win,
                                eff_win, csize, tox_dist, eff_dist, accrual, susp, u11, u10 = 0, u01 = 100, u00,   # u01, u10: default values
                                tox_dist_hyper, eff_dist_hyper, use_susp, accrual_random, pqcorr)

    patDT <- TITE_one$patDT
    doseDT <- TITE_one$doseDT
    time_current <- TITE_one$time_next           # time for making decision for the next cohort and also ready to assign the next patient

    # to decide the dose level for the next cohort, we only need doseDT
    TITE_PKBOIN12_onestep <- TITE_PKBOIN12_one_V2(doseDT, current, pT, qE, lambda_e, lambda_d, zeta1, csize, decisionM, ub, u11, u00)

    doseDT <- TITE_PKBOIN12_onestep$doseDT
    record[i] <- current
    current <- TITE_PKBOIN12_onestep$newdose   # next cohort dose

    # check the next dose exists or not
    if(is.na(current))
    {
      earlystop <- TRUE
      break
    } else if(doseDT$n[current] + csize > doselimit){  # check the next dose exceeds the dose limit or not
      break   # not early stop
    }
  }

  # after the trial stop, derive the doseDT_final for OBD selection
  doseDT_final <- doseDT %>% select(id, PK, tox, eff, r_d, n, keep) %>%
    left_join(., patDT %>% filter(!is.na(d)) %>% group_by(d) %>%
                summarise(x = sum(toxobs == 1), y = sum(effobs == 1)) %>% "colnames<-"(c("id", "x", "y")))
  doseDT_final$x[is.na(doseDT_final$x)] <- 0
  doseDT_final$y[is.na(doseDT_final$y)] <- 0


  # select optimal dose if more than two doses left after the trial
  if(earlystop == FALSE)
  {
    OBD <- fun_TITE_PKBOIN12_OBD(doseDT_final, pT, u11, u00, zeta1)
    DTremain <- doseDT[doseDT_final$keep == 1 & doseDT_final$n > 0, ]  # need to select doses applied to patients
    rN <- nrow(DTremain)   # how many number of doses remain for selection
  } else
  {
    OBD <- 99   # 99: early stop
    rN <- 0
  }

  trueOBDone <- trueOBD[1]   # the best one

  # with one OBD
  select_OBD <- ifelse(((trueOBDone %in% idlist) & (OBD %in% idlist) & (OBD == trueOBDone))  | ((!trueOBDone %in% idlist) & (!OBD %in% idlist)), 1, 0)
  num_at_OBD <- ifelse(trueOBDone %in% idlist, doseDT_final[trueOBDone, "n"], NA)
  risk_allocate <- ifelse(trueOBDone %in% idlist,
                          ifelse(doseDT_final[trueOBDone, "n"] < csize * cN/5, 1, 0), NA)

  num_overdose_OBD <- ifelse(trueOBDone %in% idlist, ifelse(max(plist) > (pT + 0.1), sum(doseDT_final[which(plist > (pT+0.1)), "n"]), 0), NA)
  num_overdose_nOBD <- ifelse(!trueOBDone %in% idlist, ifelse(max(plist) > (pT + 0.1), sum(doseDT_final[which(plist > (pT+0.1)), "n"]), 0), NA)

  # TITE statistics
  patDT <- patDT %>% filter(!is.na(tox_confirm))
  patDT$tox_confirm[patDT$cid == max(patDT$cid)] <- patDT$toxend[patDT$cid == max(patDT$cid)]  # end the trial until all the patient finish their assessments
  patDT$eff_confirm[patDT$cid == max(patDT$cid)] <- patDT$effend[patDT$cid == max(patDT$cid)]

  duration <- max(c(patDT$tox_confirm, patDT$eff_confirm), na.rm = T)   # the end of the trial

  return(data.frame(earlystop = earlystop, OBD = OBD, rN = rN, trueOBD = trueOBDone,
                    select_OBD = select_OBD, num_at_OBD = num_at_OBD,
                    num_overdose_OBD = num_overdose_OBD,
                    num_overdose_nOBD = num_overdose_nOBD,
                    risk_allocate = risk_allocate,
                    duration = duration) %>% cbind(., t(data.frame(doseDT_final$n)))
  )
}

TITE_PKBOIN12_one <- function(doseDT, current, pT, qE, lambda_e, lambda_d, zeta1, csize, decisionM, ub, u11, u00)
{

  n <- doseDT$n[current]
  phat <- doseDT$pi_t_hat[current]
  qhat <- doseDT$pi_e_hat[current]
  rhat <- doseDT$r_d[current]        # observed PK value
  cn <- round(n/csize)               # number of cohort assigned to the dose level
  dN <- nrow(doseDT)

  above <- ifelse(current == dN, NA, current + which(doseDT$keep[(current+1):dN] == 1)[1]) # NA or number
  below <- ifelse(current == 1, NA, current - which(doseDT$keep[(current-1):1] == 1)[1])   # NA or number

  PKmin <-  suppressWarnings(min(which(doseDT$r_d > zeta1)))   # Inf means no dose has PK value > zeta1

  lower_exist <- FALSE
  if(!is.infinite(PKmin) & !is.na(below))
  {
    if(PKmin < below)
    {
      lower_exist <- TRUE      # TRUE means the admissible set need to include doses {PKmin, ... d-2}
    }
  }

  if(phat * n >= decisionM$DU_T[cn])    # DU_T  (phat * n = x in BOIN12)
  {
    doseDT$keep[current:dN] = 0
    newdose <- below
    return(list(doseDT = doseDT, newdose = newdose))    # skip the following part
  }

  if(is.na(decisionM$DU_E[cn]) == FALSE)   # futility
  {
    if(qhat * n <= decisionM$DU_E[cn])    # (qhat * n = y in BOIN12)
      doseDT$keep[current] = 0
  }

  if(current == dN & n >= 6)    # concentration
  {
    if(pnorm(1.25 * zeta1, mean = doseDT$r_d[current], sd = doseDT$r_sd[current]/sqrt(n)) > 0.95)
    {
      doseDT$keep <- 0
      newdose <- NA
      return(list(doseDT = doseDT, newdose = newdose))
    }
  }

  if(sum(doseDT$keep == 0) == dN)   # no dose left
  {
    newdose <- NA
    return(list(doseDT = doseDT, newdose = newdose))
  }

  if(is.na(above) == FALSE)
  {
    if(n >= 9 & phat < lambda_d & doseDT$n[above] == 0)  # dose exploration
    {
      newdose <- above
      return(list(doseDT = doseDT, newdose = newdose))    # skip the following part
    }
  }

  # posterior prob, Inf means eliminated
  if(doseDT$keep[current] == 0)   # eliminate current dose
  {
    post_current <- -Inf
  } else {
    xd_current <- doseDT$x_d[current]            # different from BOIN12
    post_current <- 1 - pbeta(ub, 1 + xd_current, n + 1 - xd_current)
  }

  if(is.na(below))
  {
    post_below <- -Inf
  } else {
    xd_below <- doseDT$x_d[below]               # different from BOIN12
    post_below <- 1 - pbeta(ub, 1 + xd_below, doseDT$n[below] + 1 - xd_below)
  }

  if(is.na(above))
  {
    post_above <- -Inf
  } else{
    xd_above <-  doseDT$x_d[above]              # different from BOIN12
    post_above <- 1 - pbeta(ub, 1 + xd_above, doseDT$n[above] + 1 - xd_above)
  }

  lower_list <- NULL
  if(lower_exist == TRUE)
  {
    lower_list <- (PKmin:(below - 1))[doseDT$keep[PKmin:(below - 1)] == 1]
    if(length(lower_list) > 0)
    {
      xd_lower <- doseDT$x_d[lower_list]
      post_lower <- 1 - pbeta(ub, 1 + xd_lower, doseDT$n[lower_list] + 1 - xd_lower)
    }
  }

  if(lambda_d <= phat)  # de-escalate
  {
    if(!is.na(below))
    {
      if(rhat > zeta1 & length(lower_list) > 0)  # {PKmin, ... , d-1}
      {
        remain <- c(lower_list, below)
        post <- c(post_lower, post_below)
        newdose <- remain[which.max(post + (1:length(post)) * 1e-6)]
      } else {                                   # d-1
        newdose <- below
      }
      return(list(doseDT = doseDT, newdose = newdose))
    } else {
      if(is.na(current))                         # NA
      {
        newdose <- NA
        return(list(doseDT = doseDT, newdose = newdose))
      } else                                     # current
      {
        newdose <- current    # if current exist, then we use current
        return(list(doseDT = doseDT, newdose = newdose))
      }
    }
  } else   # not de-escalate
  {
    post <- c(post_below, post_current)
    remain <- c(below, current)
    if((n < 6 & !is.na(above)) | (n >= 6 & phat <= lambda_e & !is.na(above)))
    {
      post <- c(post, post_above)
      remain <- c(remain, above)
    }
    if(rhat > zeta1 & length(lower_list) > 0)  # include {PKmin, ... , d-2}
    {
      remain <- c(lower_list, remain)
      post <- c(post_lower, post)
    }

    newdose <- remain[which.max(post + (1:length(post)) * 1e-6)]
    if(max(post) == -Inf)
    {
      return(list(doseDT = doseDT, newdose = NA))
    } else
    {
      return(list(doseDT = doseDT, newdose = newdose))
    }
  }

}

