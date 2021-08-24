# replication_ape_logit

Replication Package for "Conditional inference and bias reduction for partial
effects estimation of fixed-effects logit models"

Francesco Bartolucci, Claudia Pigini and Francesco Valentini

email: f.valentini@univpm.it

SOFTWARE: R, version 4.0.5

Folders:

- simulations/: files for the replication of the Monte Carlo experiments in Tables 1 - 3;
- packages/: extra R packages, source files;
- application/: files for the replication of the empirical application in Tables 4 and 5.

Folders include a README file with details.

########## Makefile #########

simu_stat: reporduce results in Table 1, for T = 4 and T = 8;

simu_stat_12: reporduce results in Table 1, for T = 12 (this set was run separately to speed up computations);

simu_dyn: reporduce results in Table 2 and 3;

emp_app: reproduce results in Tables 4 and 5.

##############################
