## --------------------------------------------------------------------
## Replication Package for "Conditional inference and bias reduction for partial
## effects estimation of fixed-effects logit models"

## Francesco Bartolucci, Claudia Pigini and Francesco Valentini
##--------------------------------------------------------------------

rm(list=ls())


require(cquad)
source("jackknife_logit.R")
source("cml_bc_ape.R")
require(bife)


## DATASET
dati = read.csv("psid_emp_app.csv") ## load data


### Getting variables
id = dati$id2
year = dati$year

lfp = dati$part
l_lfp = dati$part1

nc02 = dati$kd2
nc35 = dati$kd5
nc617 = dati$kd17
hus_l_inc = dati$tempinc


#### Prepare lists

X = cbind(nc02,nc35,nc617,hus_l_inc)

y = lfp

Xp = X
Xp = Xp[year!=79,]

idp = id[year!=79]
yp = y[year!=79]


### STATIC LOGIT MODEL

## ML
outml =  jackknife_logit(idp, yp, Xp, jk = TRUE, dyn = FALSE)
mlbe = outml$be
mlse = outml$se_be

## JK
jkbe = outml$be_jk

##CML
mod0 = cquad_basic(idp,yp,Xp)
be0 = mod0$coefficients
se0 = mod0$se
scv = mod0$scv
J0 = mod0$J
## CML - BC APE
apef = cquad_bc_ape(idp,yp,Xp,be0,scv,J0,dyn=FALSE)

## Fernandez-Vàl BC
datibife = data.frame(dati[dati$year!=79,])

bife_stat  <- bife(part ~ kd2 + kd5 + kd17 + tempinc | id2, data=datibife)
bife_stat_bc  <- bias_corr(bife_stat)

apes_stat_bc <- get_APEs(bife_stat_bc)


## Print results

rowsn = c("nc02","nc35","nc617","hus_l_inc")
coln = c("coeff","se","pval")

cat ("STATIC LOGIT MODEL: PARAMETERS\n")
cat("ML\n")
mlc = round(cbind(mlbe,mlse,2*pnorm(abs(mlbe/mlse),lower.tail=FALSE)),digits=3)
rownames(mlc) <- rowsn
colnames(mlc) <- coln
print(mlc)

cat("JK\n")
jkc = round(cbind(jkbe,mlse,2*pnorm(abs(jkbe/mlse),lower.tail=FALSE)),digits=3)
rownames(jkc) <- rowsn
colnames(jkc) <- coln
print(jkc)


cat("CML\n")
cmlc = round(cbind(be0,se0,2*pnorm(abs(be0/se0),lower.tail=FALSE)),digits=3)
rownames(cmlc) <- rowsn
colnames(cmlc) <- coln
print(cmlc)

cat("BC\n")
print(summary(bife_stat_bc))
cat("\n")



cat("STATIC LOGIT MODEL: APE\n")
cat("ML\n")
ml = round(cbind(outml$ape,outml$se_ape,2*pnorm(abs(outml$ape/outml$se_ape),lower.tail=FALSE)),digits=3)
rownames(ml) <- rowsn
colnames(ml) <- coln
print(ml)

cat("JK\n")
jk = round(cbind(outml$ape_jk,outml$se_ape_jk,2*pnorm(abs(outml$ape_jk/outml$se_ape_jk),lower.tail=FALSE)),digits=3)
rownames(jk) <- rowsn
colnames(jk) <- coln
print(jk)

cat("CML-BC\n")
cml = round(cbind(apef$ape,apef$se,2*pnorm(abs(apef$ape/apef$se),lower.tail=FALSE)),digits=3)
rownames(cml) <- rowsn
colnames(cml) <- coln
print(cml)


cat("BC\n")
print(summary(apes_stat_bc))



### DYNAMIC LOGIT MODEL

Xp = cbind(X,l_lfp)
Xp = Xp[year!=79,]


## ML
outml =  jackknife_logit(idp, yp, Xp, jk = TRUE, dyn = TRUE)
mlbe = outml$be
mlse = outml$se_be

## JK
jkbe = outml$be_jk


## CML
mod0 = cquad_pseudo(id,y,X)
be0 = mod0$coefficients
se0 = mod0$se
scv = mod0$scv
J0 = mod0$J
## CML - BC APE
apef = cquad_bc_ape(idp,yp,Xp,be0,scv,J0,dyn=TRUE)


## Fernandez-Vàl
datibife = data.frame(dati[dati$year!=79,])

bife_dyn  <- bife(part ~ kd2 + kd5 + kd17 + tempinc + part1 | id2, data=datibife)
bife_dyn_bc  <- bias_corr(bife_dyn, L = 1L)

apes_dyn_bc <- get_APEs(bife_dyn_bc)

## Print results

drowsn = c("nc02","nc35","nc617","hus_l_inc","l_lfp")
dcoln = c("coeff","se","pval")

cat ("DYNAMIC LOGIT MODEL: PARAMETERS\n")
cat("ML\n")
dmlc = round(cbind(mlbe,mlse,2*pnorm(abs(mlbe/mlse),lower.tail=FALSE)),digits=3)
rownames(dmlc) <- drowsn
colnames(dmlc) <- dcoln
print(dmlc)

cat("JK\n")
djkc = round(cbind(jkbe,mlse,2*pnorm(abs(jkbe/mlse),lower.tail=FALSE)),digits=3)
rownames(djkc) <- drowsn
colnames(djkc) <- dcoln
print(djkc)
      
cat("CML-BC\n")
dcmlc = round(cbind(be0,se0,2*pnorm(abs(be0/se0),lower.tail=FALSE)),digits=3)
rownames(dcmlc) <- drowsn
colnames(dcmlc) <- dcoln
print(dcmlc)


cat("BC\n")
print(summary(bife_dyn_bc))
cat("\n")


cat ("DYNAMIC LOGIT MODEL: APE\n")
cat("ML\n")
dml = round(cbind(outml$ape,outml$se_ape,2*pnorm(abs(outml$ape/outml$se_ape),lower.tail=FALSE)),digits=3)
rownames(dml) <- drowsn
colnames(dml) <- dcoln
print(dml)

cat("JK\n")
djk = round(cbind(outml$ape_jk,outml$se_ape_jk,2*pnorm(abs(outml$ape_jk/outml$se_ape_jk),lower.tail=FALSE)),digits=3)
rownames(djk) <- drowsn
colnames(djk) <- dcoln
print(djk)

cat("CML-BC\n")
dcml = round(cbind(apef$ape,apef$se,2*pnorm(abs(apef$ape/apef$se),lower.tail=FALSE)),digits=3)
rownames(dcml) <- drowsn
colnames(dcml) <- dcoln
print(dcml)


cat("BC\n")
print(summary(apes_dyn_bc))

