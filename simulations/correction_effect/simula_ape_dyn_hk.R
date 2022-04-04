## --------------------------------------------------------------------
## Replication Package for "Conditional inference and bias reduction for partial
## effects estimation of fixed-effects logit models"

## Francesco Bartolucci, Claudia Pigini and Francesco Valentini
##--------------------------------------------------------------------

rm(list=ls())

require(cquad)

#source("jackknife_logit_bife.R")

source("cml_bc_ape.R")                              
require(bife)

set.seed(1022022)
## set simulation parameters
size = c(100, 500)
time = c(5,9)
nit = 1000
ga = 0.5
ncov = 1
be = 1
sdx = pi/sqrt(3)

for(n in size){
    for(T in time){

        filename = sprintf("APE_cml_compare_dyn_N%g_T%g_ga%g_HK_IT%g",n, T, ga, nit)

        APECMLF = APECML = matrix(0,nit,2) 
        
        TMP = matrix(0,nit,2)

                                       
        it = 1
        while(it<=nit){
                                        # generate data
            id = (1:n)%x%rep(1,T); idp = (1:n)%x%rep(1,(T-1))
            label = unique(id)
            
            X = y = Xp  = yp = eta  = PEg = PEb = sel = alphavec =  ap  = NULL
            for(i in 1:n){
                
                xi = rnorm(T,0,sdx)
                al = mean(xi[1:4])
                alphavec = c(alphavec,al)
                eta = be*xi[1] + al
                pi = exp(eta)/(1+exp(eta))
                yi = 1*runif(1)<pi
                
                for(t in 2:T){
                    eta =  al + ga*yi[t-1] + be*xi[t] 
                    p = exp(eta)/(1+exp(eta)); pi = c(pi,p)
                    yi = c(yi, 1*runif(1)<p)
                }

                pebe_i = mean(pi[2:T]*(1-pi[2:T]))*be
                                       
                eta = al + be*xi[2:T]
                p0 = exp(eta)/(1+exp(eta))
                p1 = exp(eta + ga)/(1+exp(eta + ga))                
                pega_i = mean(p1 - p0)
                                     
                PEg = c(PEg,pega_i); PEb = c(PEb,pebe_i);

                X = c(X, xi); y = c(y, yi)
                Xp = rbind(Xp, cbind(xi[2:T], yi[1:(T-1)]))
                yp = c(yp,yi[2:T])
                ap = c(ap,al)
                
            }
            ap = ap%x%rep(1,T-1)
            
                                        
            TMP[it,1] = mean(PEb); TMP[it,2] = mean(PEg)
            X = as.matrix(X); Xp = as.matrix(Xp)
            
            ## Estimation
            
            cat("CMLF\n")
            ##  ## CML + ML
            mod0 = cquad_pseudo(id,y,X)
            be0 = mod0$coefficients
            se0 = mod0$ser
            
            scv = mod0$scv
            J0 = mod0$J

            ##CML -BC

            apef = cquad_bc_ape(idp,yp,Xp,be0,scv,J0,dyn=TRUE)
            scarti = 0
            APECMLF[it,1] = apef$ape[1]
            APECMLF[it,2] = apef$ape[2]

            APECML[it,1] = apef$apeu[1]
            APECML[it,2] = apef$apeu[2]

            
            if(1){
                
                
                
                

                
                if(it>1){
                                       
                    
                    print("-----------------------------------------")
                    cat(c("it  =",it,"\n"))                        
                    cat(c("n =", n,"\n"))
                    cat(c("T =", T,"\n"))
                  

                    for(j in 1:2){
                        RES = matrix(rep(NA),2,2)
                        
                        rows = c("ape cml", "ape cml_bc")
                        cols = c("mean ratio", "med ratio")
                        
                        rownames(RES) = rows
                        colnames(RES) = cols
                        
                       
                        RES[1,1] = mean(APECML[1:it,j]/TMP[1:it,j], na.rm=T)
                        RES[1,2] = median(APECML[1:it,j]/TMP[1:it,j], na.rm=T)

                        RES[2,1] = mean(APECMLF[1:it,j]/TMP[1:it,j], na.rm=T)
                        RES[2,2] = median(APECMLF[1:it,j]/TMP[1:it,j], na.rm=T) 

                        



                        
                        cat("APE\n")
                        print(round(RES, digits=3), na.print = "" , quote = FALSE)
                       
                        
                    }
                }
                it = it + 1
            }else{scarti = scarti+1}
            
        }
        save.image(filename)
        
    }
}
