#Replication Package for "Conditional inference and bias reduction for partial
#effects estimation of fixed-effects logit models"

#Francesco Bartolucci, Claudia Pigini and Francesco Valentini

## Empirical Application
emp_app:
	cd application/; Rscript ape_appli.R


## Simulations

## Static Logit

## Table 1
simu_stat:
	cd simulations/static_logit_hk/static ; Rscript simula_ape_stat_hk.R

simu_stat_12:
	cd  simulations/static_logit_hk/static12; Rscript simula_ape_stat_hk.R

## Dynamic Logit

## Tables 2 and A1
simu_dyn_025:
	cd simulations/dynamic_logit_hk/gamma_025/dynamic_025 ; Rscript simula_ape_dyn_hk.R

simu_dyn_025_12:
	cd simulations/dynamic_logit_hk/gamma_025/dynamic12_025 ; Rscript simula_ape_dyn_hk.R

## Tables 3 and A2
simu_dyn_05:
	cd simulations/dynamic_logit_hk/gamma_05/dynamic ; Rscript simula_ape_dyn_hk.R

simu_dyn_05_12:
	cd simulations/dynamic_logit_hk/gamma_05/dynamic12 ; Rscript simula_ape_dyn_hk.R

## Tables 4 and A3
simu_dyn_075:
	cd simulations/dynamic_logit_hk/gamma_075/dynamic_075 ; Rscript simula_ape_dyn_hk.R

simu_dyn_075_12:
	cd simulations/dynamic_logit_hk/gamma_075/dynamic12_075 ; Rscript simula_ape_dyn_hk.R

## Table 5
bc_vs_nobc:
	cd simulations/correction_effect ; Rscript simula_ape_stat_hk.R; Rscript simula_ape_dyn_hk.R

###### Robustness #########

## Table A4
time-dummies:
	cd simulations/time-dummies ; Rscript simula_ape_dyn_hk.R

## Table A5
nerlove:
	cd simulations/nerlove_process/nerlove ; Rscript simula_ape_stat_hn.R

nerlove12:
	cd simulations/nerlove_process/nerlove12 ; Rscript simula_ape_stat_hn.R


## Table A6
chi2alpha:
	cd simulations/chi2_alpha ; Rscript simula_ape_dyn_hk.R


## Table A7
rare:
	cd simulations/rare-events ; Rscript simula_ape_dyn_hk.R
