## --------------------------------------------------------------------
## Replication Package for "Conditional inference and bias reduction for partial
## effects estimation of fixed-effects logit models"

## Francesco Bartolucci, Claudia Pigini and Francesco Valentini
##--------------------------------------------------------------------

rm(list=ls())

require(cquad)
#source("jackknife_logit.R") 
source("cml_bc_ape.R")
require(bife) 


set.seed(20210512)
## set simulation parameters
size = c(100, 500)
time = c(4,8) 
nit = 1000
be = 1
sdx = pi/sqrt(3)

for(n in size){
    for(T in time){
        filename = sprintf("APE_cml_compare_stat_N%g_T%g_HK_IT%g",n, T, nit)

            APECMLF = APECML = matrix(0,nit,1)
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
                    
                                        #if(sum(yi)==0 | sum(yi)==T) pe_i = 0
                    PE = c(PE,pe_i)
                }

                #browser()
                
                TMP[it] = mean(PE)
                X = as.matrix(X)


                ## CML + BC
                mod0 = cquad_basic(id,y,X)
                be0 = mod0$coefficients
                scv = mod0$scv
                J0 = mod0$J
                apef = cquad_bc_ape(id,y,X,be0,scv,J0,dyn=FALSE)
                APECMLF[it] = apef$ape
                APECML[it] =  apef$apeu


                if(it>1){
                    print("-----------------------------------------")
                    cat(c("it  =",it,"\n"))                        
                    cat(c("n =", n,"\n"))
                    cat(c("T =", T,"\n"))
                    
                    RES = matrix(rep(NA),2,2)
                    
                    rows = c("ape cml", "ape cml_bc")
                    cols = c("mean ratio", "med ratio")
                    
                    rownames(RES) = rows
                    colnames(RES) = cols


                    

                    RES[1,1] = mean(APECML[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[1,2] = median(APECML[1:it,1]/TMP[1:it,1], na.rm=T)

                    RES[2,1] = mean(APECMLF[1:it,1]/TMP[1:it,1], na.rm=T)
                    RES[2,2] = median(APECMLF[1:it,1]/TMP[1:it,1], na.rm=T)


                    print(round(RES, digits=3), na.print = "" , quote = FALSE)



                    
                }
            }
            save.image(filename)
       
    }
}

