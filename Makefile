#Replication Package for "Conditional inference and bias reduction for partial
#effects estimation of fixed-effects logit models"

#Francesco Bartolucci, Claudia Pigini and Francesco Valentini

emp_app:
	cd application/; Rscript ape_appli.R

simu_stat:
	cd simulations/static/ ; Rscript simula_ape_stat_hn.R

simu_stat_12:
	cd simulations/static12/; Rscript simula_ape_stat_hn.R

simu_dyn:
	cd simulations/dynamic/; Rscript simula_ape_dyn_hk.R

