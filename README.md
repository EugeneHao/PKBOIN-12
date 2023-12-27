# PKBOIN-12

Simulation Example: 

```r
library(dplyr)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(inline)
require(truncnorm)
library(knitr)

dN = 6;
pT = 0.35; qE = 0.25; pqcorr = 0;
psi0PK = 6000; CV = 0.25; g_P = 1; 
csize = 3; cN = 15;
u11 = 60; u00 = 40;
cutoff_tox = 0.95; cutoff_eff = 0.9;
repsize = 2000; n_cores = 10;
accrual = 10; tox_win = 30; eff_win = 60;

rV = c(1000, 1500, 2500, 3600, 4800, 6500)    # PK
pV = c(0.01, 0.03, 0.05, 0.10, 0.18, 0.24)    # tox
qV = c(0.05, 0.10, 0.20, 0.30, 0.45, 0.55)    # eff

trueOBD <- findOBD_RDS(pV, qV, pT, qE, u11 = 60, u00 = 40)

PKBOIN12_result <- 
  fun_TITE_PK_fixsimu(dN = dN,
                rV = rV, pV = pV, qV = qV, 
                pT = pT, qE = qE, pqcorr = 0, psi0PK = psi0PK, CV = CV, g_P = g_P, 
                csize = csize, cN = cN, design = "PKBOIN-12",
                alphaT = 1, betaT = 1, alphaE = 1, betaE = 1,
                alphaTO = 0.5, betaTO = 0.5, alphaEO = 0.5, betaEO = 0.5, 
                current = 1, doselimit = Inf,
                u11 = u11, u00 = u00, cutoff_tox = cutoff_tox, cutoff_eff = cutoff_eff,
                repsize = repsize, n_cores = n_cores,
                accrual = accrual, susp = 1, tox_win = tox_win, eff_win = eff_win, 
                tox_dist = "Uniform", eff_dist = "Uniform",
                tox_dist_hyper = NULL, eff_dist_hyper = NULL, use_susp = T,
                accrual_random = F, considerPK = T)

TITE_PKBOIN12_result <- 
  fun_TITE_PK_fixsimu(dN = dN,
                rV = rV, pV = pV, qV = qV, 
                pT = pT, qE = qE, pqcorr = 0, psi0PK = psi0PK, CV = CV, g_P = g_P, 
                csize = csize, cN = cN, design = "TITE-PKBOIN-12",
                alphaT = 1, betaT = 1, alphaE = 1, betaE = 1,
                alphaTO = 0.5, betaTO = 0.5, alphaEO = 0.5, betaEO = 0.5, 
                current = 1, doselimit = Inf,
                u11 = u11, u00 = u00, cutoff_tox = cutoff_tox, cutoff_eff = cutoff_eff,
                repsize = repsize, n_cores = n_cores,
                accrual = accrual, susp = 0.5, tox_win = tox_win, eff_win = eff_win, 
                tox_dist = "Uniform", eff_dist = "Uniform",
                tox_dist_hyper = NULL, eff_dist_hyper = NULL, use_susp = T,
                accrual_random = F, considerPK = T)
```
