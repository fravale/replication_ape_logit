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


set.seed(20210512)
## set simulation parameters
size = c(100, 500)
time = 12 
nit = 1000
be = 1
sdx = pi/sqrt(3)

for(n in size){
    for(T in time){
        filename = sprintf("APEbc_stat_N%g_T%g_HK_IT%g",n, T, nit)
            APEML =  matrix(0,nit,1); SEML = matrix(0,nit,1); INFBML = INFAML = matrix(0,nit,1)
            APEJK   = matrix(0,nit,1);   SEJK = matrix(0,nit,1)
            APEBC  = matrix(0,nit,1); SEBC = matrix(0,nit,1)
            APECMLF = matrix(0,nit,1); SECMLF = matrix(0,nit,1);
            BCPAR = CMLPAR = matrix(0,nit,1)
            TMP  = matrix(0,nit,1)


        for(it in 1:nit){## generate data
                id = (1:n)%x%rep(1,T)
                label = unique(id)
                X = y = eta = alpha = PE =  indr = NULL

                for(i in 1:n){
                    Xi = rnorm(T,0,sdx)
                    
                    X = c(X,Xi)
                    alphai = mean(Xi[1:4])*rep(1,T)
                    alpha = c(alpha, alphai)
                    
                    
                    etai = Xi*be + alphai
                    eta = c(eta, etai)        
                    p = exp(etai)/(1+exp(etai))
                    yi = 1*(runif(T)<p); y = c(y,yi)                
                    pe_i = mean(p*(1-p)*be)
                    
                                       
                    PE = c(PE,pe_i)
                }

  
                
                TMP[it] = mean(PE)
                X = as.matrix(X)
                print(mean(y))
                ## Estimation
                                      
                ## ML
                mod_jk = jackknife_logit(id, y, X, jk = TRUE, dyn = FALSE)
                APEML[it] = mod_jk$ape
                ;  SEML[it] = mod_jk$se_ape
                beml = mod_jk$be; alml = mod_jk$AL
                
                ## JACKKNIFE
                APEJK[it] = mod_jk$ape_jk;  SEJK[it] = mod_jk$se_ape_jk



                ##ANALYTICAL BIAS CORRECTION
                datibife = data.frame(id,y,X)
                bife_stat  <- bife(y ~ X | id, data=datibife, control=list(dev_tol = 10^(-12)))
                bife_stat_bc  <- bias_corr(bife_stat)

                apes_stat_bc <- get_APEs(bife_stat_bc)
                
 
                APEBC[it] = apes_stat_bc$delta               
                SEBC[it] = sqrt(diag(apes_stat_bc$vcov))


                ## CML + PENALIZED ML (FIRTH)
                mod0 = cquad_basic(id,y,X)
                be0 = mod0$coefficients
                scv = mod0$scv
                J0 = mod0$J
                apef = cquad_bc_ape(id,y,X,be0,scv,J0,dyn=FALSE)
                APECMLF[it] = apef$ape; SECMLF[it] = apef$se


                if(it>1){
                    print("-----------------------------------------")
                    cat(c("it  =",it,"\n"))                        
                    cat(c("n =", n,"\n"))
                    cat(c("T =", T,"\n"))
                    
                    RES = matrix(rep(NA),4,6)
                    
                    rows = c("ape ml", "ape jk", "ape bc", "ape cml_bc")
                    cols = c("mean ratio", "med ratio", "sd", "p .05", "p .10",
                             "se/sd")
                    
                    rownames(RES) = rows
                    colnames(RES) = cols


                    RES[1,1] = mean(APEML[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[1,2] = median(APEML[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[1,3] = sd(APEML[1:it,1])
                    
                    CIl1 = APEML[1:it,1] - SEML[1:it,1]*1.96
                    CIh1 = APEML[1:it,1] + SEML[1:it,1]*1.96
                    RES[1,4] = mean(TMP[1:it,1] > CIh1 | TMP[1:it,1] < CIl1)
                    
                    CIl1 = APEML[1:it,1] - SEML[1:it,1]*1.65
                    CIh1 = APEML[1:it,1] + SEML[1:it,1]*1.65
                    RES[1,5] = mean(TMP[1:it,1] > CIh1 | TMP[1:it,1] < CIl1)

                    RES[1,6] = mean(SEML[1:it,1])/sd(APEML[1:it,1])


                    RES[2,1] = mean(APEJK[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[2,2] = median(APEJK[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[2,3] = sd(APEJK[1:it,1])
                    
                    CIl1 = APEJK[1:it,1] - SEJK[1:it,1]*1.96
                    CIh1 = APEJK[1:it,1] + SEJK[1:it,1]*1.96
                    RES[2,4] = mean(TMP[1:it,1] > CIh1 | TMP[1:it,1] < CIl1)
                    
                    CIl1 = APEJK[1:it,1] - SEJK[1:it,1]*1.65
                    CIh1 = APEJK[1:it,1] + SEJK[1:it,1]*1.65
                    RES[2,5] = mean(TMP[1:it,1] > CIh1 | TMP[1:it,1] < CIl1)

                    RES[2,6] = mean(SEJK[1:it,1])/sd(APEJK[1:it,1])


                    RES[3,1] = mean(APEBC[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[3,2] = median(APEBC[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[3,3] = sd(APEBC[1:it,1])
                    
                    CIl1 = APEBC[1:it,1] - SEBC[1:it,1]*1.96
                    CIh1 = APEBC[1:it,1] + SEBC[1:it,1]*1.96
                    RES[3,4] = mean(TMP[1:it,1] > CIh1 | TMP[1:it,1] < CIl1)
                    
                    CIl1 = APEBC[1:it,1] - SEBC[1:it,1]*1.65
                    CIh1 = APEBC[1:it,1] + SEBC[1:it,1]*1.65
                    RES[3,5] = mean(TMP[1:it,1] > CIh1 | TMP[1:it,1] < CIl1)

                    RES[3,6] = mean(SEBC[1:it,1])/sd(APEBC[1:it,1])

                    

                    RES[4,1] = mean(APECMLF[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[4,2] = median(APECMLF[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[4,3] = sd(APECMLF[1:it,1])
                    
                    CIl1 = APECMLF[1:it,1] - SECMLF[1:it,1]*1.96
                    CIh1 = APECMLF[1:it,1] + SECMLF[1:it,1]*1.96
                    RES[4,4] = mean(TMP[1:it,1] > CIh1 | TMP[1:it,1] < CIl1)
                    
                    CIl1 = APECMLF[1:it,1] - SECMLF[1:it,1]*1.65
                    CIh1 = APECMLF[1:it,1] + SECMLF[1:it,1]*1.65
                    RES[4,5] = mean(TMP[1:it,1] > CIh1 | TMP[1:it,1] < CIl1)

                    RES[4,6] = mean(SECMLF[1:it,1])/sd(APECMLF[1:it,1])

                    



                    print(round(RES, digits=3), na.print = "" , quote = FALSE)



                    
                }
            }
            save.image(filename)
       
    }
}

