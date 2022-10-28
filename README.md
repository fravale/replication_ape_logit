# replication_ape_logit

Replication Package for "Conditional inference and bias reduction for partial
effects estimation of fixed-effects logit models"

Francesco Bartolucci, Claudia Pigini and Francesco Valentini

Empirical Economics (2022), https://doi.org/10.1007/s00181-022-02313-6

email: f.valentini@univpm.it

SOFTWARE: R, version 4.0.5

Folders:

- simulations/: files for the replication of the Monte Carlo experiments in Tables 1 - 5 and Appendices.
- packages/: extra R packages, source files;
- application/: files for the replication of the empirical application in Tables 6 and 7.

Folders include a README file with details.

########## Makefile #########

simu_stat: replicates results in Table 1, for T = 4 and T = 8;

simu_stat_12: replicates results in Table 1, for T = 12 (this set was
run separately to speed up computations);

simu_dyn_025: replicates results in Tables 2 and 8, for T = 4 and T = 8;

simu_dyn_025_12: replicates results in Tables 2 and 8, for T = 12
(this set was run separately to speed up computations);

simu_dyn_05: replicates results in Tables 3 and 9, for T = 4 and T = 8;

simu_dyn_05_12: replicates results in Tables 3 and 9, for T = 12
(this set was run separately to speed up computations);

simu_dyn_075: replicates results in Tables 4 and 10, for T = 4 and T = 8;

simu_dyn_075_12: replicates results in Tables 4 and 10, for T = 12
(this set was run separately to speed up computations);

bc_vs_nobc: replicates results repoted in Table 5 (Sequence: static,
dynamic)

time-dummies: replicates results repoted in Table 11

nerlove: replicates results in Table 12, for T = 4 and T = 8;

nerlove12: replicates results in Table 12, for T = 12 (this set was
run separately to speed up computations);

chi2alpha: replicates results in Table 13;

rare: replicates results in Table 14.

##############################
