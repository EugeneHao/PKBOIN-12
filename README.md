# Simulation Example for PKBOIN-12 and TITE-PKBOIN-12

This document provides guidance on how to use our code to conduct simulations with PKBOIN-12 and TITE-PKBOIN-12.

```r
# library all required packages: 
library(dplyr)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(inline)
require(truncnorm)
library(knitr)

# source functiosn: 
source("fun_findOBD.R")
source("fun_pava.R")
source("fun_TITE_PKBOIN12dec.R")
source("fun_TITE_PK_core_para.R")
source("fun_TITE_PK_fixsimu.R")
source("fun_TITE_PK_update.R")
source("fun_TITE_PKBOIN12_OBD.R")
source("fun_TITE_PKBOIN12.R")
source("fun_TITE_PKBOIN12_one.R")

# set design parameters: 
dN = 6;         # number of dose levels 
pT = 0.35;      # maximum acceptable toxicity rate 
qE = 0.25;      # minimum acceptable efficacy rate
pqcorr = 0;     # correlation between toxicity and efficacy outcomes  
psi0PK = 6000;  # target PK value
CV = 0.25;      # coefficient of variation 
g_P = 1;        # ratio that measures the relationship of PK outcome with toxcity and efficacy probabilities
csize = 3;      # number of patients per cohort 
cN = 15;        # maximum number of cohorts 
u11 = 60; u00 = 40;    # utility scores for BOIN12 and TITE-BOIN12
cutoff_tox = 0.95;     # toxicty early termination cutoff 
cutoff_eff = 0.9;      # efficacy early termination cutoff 
repsize = 2000;        # replication size 
n_cores = 10;          # number of cores used for computation 
accrual = 10;          # arrural rate 
tox_win = 30;          # toxicty assessment window
eff_win = 60;          # efficacy assessment window

# Specify the true values: 
rV = c(1000, 1500, 2500, 3600, 4800, 6500)    # PK 
pV = c(0.01, 0.03, 0.05, 0.10, 0.18, 0.24)    # toxicity
qV = c(0.05, 0.10, 0.20, 0.30, 0.45, 0.55)    # efficacy

# identify the true OBD with RDS 
trueOBD <- findOBD_RDS(pV, qV, pT, qE, u11, u00)

```

The main function for running the simulation is called `fun_TITE_PK_fixsimu`. Below are some important arguments for this function: 

+ `design`: specifies the model-assisted design used for the simulation. The options are PKBOIN-12 and TITE-PKBOIN-12
+ `current`: initial dose level.
+ `doselimit`: determines the maximum number of patients assigned to a single dose before terminating the trial. If this limit is not reached, the trial continues until the total sample size reaches csize * cN. By default, doselimit = Inf.
+ `susp`: cutoff point for suspension. For TITE-BOIN12 and TITE-PKBOIN-12, `susp = 0.5`, meaning accrual is suspended to wait for more data if the ascertained toxicity or efficacy outcomes of current dose do not reach 50%. 
+ `tox_dist`: distribution of time to toxicity outcome.
+ `eff_dist`: distribution of time to efficacy outcome.
+ `use_susp`: if true, the suspension rule is applied.
+ `accrual_random`: if true, the accrual time for enrolling the next patient is randomly generated from $U(0, 2 * accrual)$.
+ `considerPK`: if true, consider the PK outcome. This is set to TRUE for PKBOIN-12 and TITE-PKBOIN-12
 
```r
# PKBOIN-12: 
PKBOIN12_result <- 
  fun_TITE_PK_fixsimu(dN = dN,
                rV = rV, pV = pV, qV = qV, 
                pT = pT, qE = qE, pqcorr = 0, psi0PK = psi0PK, CV = CV, g_P = g_P, 
                csize = csize, cN = cN, design = "PKBOIN-12",
                current = 1, doselimit = Inf,
                u11 = u11, u00 = u00, cutoff_tox = cutoff_tox, cutoff_eff = cutoff_eff,
                repsize = repsize, n_cores = n_cores,
                accrual = accrual, susp = 1, tox_win = tox_win, eff_win = eff_win, 
                tox_dist = "Uniform", eff_dist = "Uniform",
                tox_dist_hyper = NULL, eff_dist_hyper = NULL, use_susp = F,
                accrual_random = F, considerPK = T)

# TITE-PKBOIN-12 
TITE_PKBOIN12_result <- 
  fun_TITE_PK_fixsimu(dN = dN,
                rV = rV, pV = pV, qV = qV, 
                pT = pT, qE = qE, pqcorr = 0, psi0PK = psi0PK, CV = CV, g_P = g_P, 
                csize = csize, cN = cN, design = "TITE-PKBOIN-12",
                current = 1, doselimit = Inf,
                u11 = u11, u00 = u00, cutoff_tox = cutoff_tox, cutoff_eff = cutoff_eff,
                repsize = repsize, n_cores = n_cores,
                accrual = accrual, susp = 0.5, tox_win = tox_win, eff_win = eff_win, 
                tox_dist = "Uniform", eff_dist = "Uniform",
                tox_dist_hyper = NULL, eff_dist_hyper = NULL, use_susp = T,
                accrual_random = F, considerPK = T)
```

The output of running `fun_TITE_PK_fixsimu` is a `data.frame` with `repsize` rows. The variables included in the output are as follows:

+ `n`: maximum sample size
+ `design`: the specified design method  
+ `trueOBD`: the true OBD based on the true probability values given $u_{11}$ and $u_{00}$. A value of -1 indicates that none of the dose levels have a true OBD
+ `OBD`: the selected OBD, where `OBD = 99` means early termination without selecting any dose as the OBD
+ `rN`: the number of remaining dose levels for OBD selection 
+ `select_OBD`: indicates whether the true OBD was correctly identified
+ `num_at_OBD`: the number of patients assigned to the OBD
+ `risk_allocate`: a value of 1 indicates that the trial assigned less than 20% of patients to the true OBD
+ `num_overdose_OBD`: the number of patients assigned to the toxic doses if the true OBD exists
+ `num_overdose_nOBD`: the number of patients assigned to the toxic doses if the true OBD does not exist
+ `duration`: the trial duration in days
+ `1`, ..., `6`: the number of patients assigned to dose level 1 through 6 


Next we combine the above results of all designs and summarize the outputs: 

```r
result <- rbind(PKBOIN12_result, TITE_PKBOIN12_result)

result_summary <- 
   result %>% 
   group_by(design) %>% 
   summarise(OBD99 = mean(OBD == 99, na.rm = T) * 100,
             OBDother = mean(OBD == trueOBD, na.rm = T) * 100,
             p1 = mean(OBD == 1, na.rm = T) * 100, p2 = mean(OBD == 2, na.rm = T) * 100,
             p3 = mean(OBD == 3, na.rm = T) * 100, p4 = mean(OBD == 4, na.rm = T) * 100, 
             p5 = mean(OBD == 5, na.rm = T) * 100, p6 = mean(OBD == 6, na.rm = T) * 100, 
             n1 = mean(`1`, na.rm = T), n2 = mean(`2`, na.rm = T), 
             n3 = mean(`3`, na.rm = T), n4 = mean(`4`, na.rm = T), 
             n5 = mean(`5`, na.rm = T), n6 = mean(`6`, na.rm = T),
             dur = mean(duration, na.rm = T)/30) %>% 
   mutate(p_OBD = ifelse(rep(trueOBD, 2) <= 0, OBD99, OBDother), 
          p_rej = 100 - p1 - p2 - p3 - p4 - p5 - p6) %>% 
   select(-OBD99, -OBDother) %>% "["(,c(1, 15, 2:7, 16, 8:14))
   
result_summary
```

```
# A tibble: 2 × 16
  design         p_OBD    p1    p2    p3    p4    p5    p6 p_rej    n1    n2    n3    n4    n5    n6   dur
  <chr>          <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
1 PKBOIN-12       53.8     0     0  0.1    3.3  40.9  53.8  1.90  3.05  3.32  5.05  8.20  12.6  12.6  38.1
2 TITE-PKBOIN-12  52.8     0     0  0.05   3.2  41.9  52.8  2     3.21  3.58  5.39  8.49  12.7  11.4  25.0
```

The variables included in the summarized output are as follows:

+ `p_OBD`: the probability of correctly identifying the true OBD or correctly terminating the trial without selecting any dose as the OBD if the trueOBD does not exist
+ `pi`: the probability of selecting dose level $i$ as the OBD 
+ `p_rej`: the probability of terminating the trial without selecting any dose as the OBD
+ `ni`: the averaged number of patients assigned to dose level $i$ 
+ `dur`: the trial duration in months 




