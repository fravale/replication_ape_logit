## --------------------------------------------------------------------
## Replication Package for "Conditional inference and bias reduction for partial
## effects estimation of fixed-effects logit models"

## Francesco Bartolucci, Claudia Pigini and Francesco Valentini
##--------------------------------------------------------------------

rm(list=ls())

require(cquad)

source("jackknife_logit_bife.R")

source("cml_bc_ape.R")                              
require(bife)

set.seed(314052021)
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

        filename = sprintf("APEbc_dyn_N%g_T%g_ga%g_HK_IT%g",n, T, ga, nit)
        APEML =  matrix(0,nit,2); SEML = matrix(0,nit,2); 
        APEJK   = matrix(0,nit,2); SEJK = matrix(0,nit,2)
        APEBC  = matrix(0,nit,2); SEBC = matrix(0,nit,2)
        APECMLF = matrix(0,nit,2); SECMLF = matrix(0,nit,2);
        
        
        TMP = matrix(0,nit,2)
        PARML = PARBC = SEPARML = SEPARBC = PARPCML = SEPARPCML = matrix(0,nit,2)

        
        
        for(it in 1:nit){
            
                                        # generate data
            id = (1:n)%x%rep(1,T); idp = (1:n)%x%rep(1,(T-1))
            idt = rep(1,n)%x%seq(1,T); idtp = rep(1,n)%x%seq(1,T-1)
            label = unique(id)
            
            X = y = Xp  = yp = eta  = PEg = PEb = sel = alphavec =  ap  = NULL

            trend_i = seq(-1,1,by=2/(T-2))
             
            for(i in 1:n){
                
                xi = rnorm(T,0,sdx)
                al = mean(xi[1:4])
                alphavec = c(alphavec,al)
                eta = be*xi[1] + al
                pi = exp(eta)/(1+exp(eta))
                yi = 1*runif(1)<pi
                
                for(t in 2:T){
                    eta =  al + ga*yi[t-1] + be*xi[t] + trend_i[t-1] 
                    p = exp(eta)/(1+exp(eta)); pi = c(pi,p)
                    yi = c(yi, 1*runif(1)<p)
                }

                pebe_i = mean(pi[2:T]*(1-pi[2:T]))*be

                
                eta = al + be*xi[2:T] + trend_i
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

            ## Include Time-dummies among covariates
            for(g in 3:T){
                X = cbind(X,as.numeric(idt==g))
            }

            for(g in 2:(T-1)){
                Xp = cbind(Xp,as.numeric(idtp==g))
            }

            
            ## Estimation
            
            cat("CMLF\n")
            ##  ## CML + ML
            mod0 = cquad_pseudo(id,y,X)
            be0 = mod0$coefficients
            se0 = mod0$ser
            #PARPCML[it,1] = be0[1]; SEPARPCML[it,1] = se0[1]
            #PARPCML[it,2] = be0[length(be0)]; SEPARPCML[it,1] = se0[length(be0)]

            #browser()
            
            scv = mod0$scv
            J0 = mod0$J

            ##CML -BC
            #Adjust order in be0
            be0a = c(be0[1],be0[length(be0)],be0[-c(1,length(be0))])


            apef = cquad_bc_ape(idp,yp,Xp,be0a,scv,J0,dyn=TRUE)


            
            APECMLF[it,1] = apef$ape[1]; SECMLF[it,1] = apef$se[1]
            APECMLF[it,2] = apef$ape[2]; SECMLF[it,2] = apef$se[2]
            
            
            
            
            cat("ML & JK\n")
            
                mod_jk = jackknife_logit(idp, yp, Xp, jk = FALSE, dyn = TRUE)
                APEML[it,1] = mod_jk$ape[1];  SEML[it,1] = mod_jk$se_ape[1]
                APEML[it,2] = mod_jk$ape[2];  SEML[it,2] = mod_jk$se_ape[2]


            cat("BC\n")
            
            ## ANALYTICAL BIAS CORRECTION
            ## FV - by bife()

            datmat = cbind(idp,yp,Xp)
            datibife = as.data.frame(datmat)

            bform = as.formula(paste("yp ~ ", paste(paste0("V", 3:(ncol(datibife))), collapse = " + "),paste("| idp")))


            bife_dyn <- bife(bform, data=datibife, control=list(dev_tol = 10^(-12)))

            
            bife_dyn_bc <- bias_corr(bife_dyn, L = 1L)
            
            
            apes_dyn_bc <- get_APEs(bife_dyn_bc)
            apes_se = sqrt(diag(apes_dyn_bc$vcov))
            
            APEBC[it,1] = apes_dyn_bc$delta[1]               
            SEBC[it,1] = apes_se[1]              
            APEBC[it,2] = apes_dyn_bc$delta[2]
            SEBC[it,2] = apes_se[2]              
            #PARBC[it,1:2] = bife_dyn_bc$coefficients[1:2]; SEPARBC[it,1:2] = sqrt(diag(solve(bife_dyn_bc$Hessian)))
            
            
            
            
            if(it>1){
                
                
                print("-----------------------------------------")
                cat(c("it  =",it,"\n"))                        
                cat(c("n =", n,"\n"))
                cat(c("T =", T,"\n"))
                

                for(j in 1:2){
                    RES = matrix(rep(NA),4,6)
                    
                    rows = c("ape ml", "ape jk", "ape bc", "ape cml_bc")
                    cols = c("mean ratio", "med ratio", "sd",
                             "p .05", "p .10","se/sd")
                    
                    rownames(RES) = rows
                    colnames(RES) = cols
                    
                    RES[1,1] = mean(APEML[1:it,j]/TMP[1:it,j], na.rm=T)
                    RES[1,2] = median(APEML[1:it,j]/TMP[1:it,j], na.rm=T)
                    RES[1,3] = sd(APEML[1:it,j])
                    
                    CIl1 = APEML[1:it,j] - SEML[1:it,j]*1.96
                    CIh1 = APEML[1:it,j] + SEML[1:it,j]*1.96
                    RES[1,4] = mean(TMP[1:it,j] > CIh1 | TMP[1:it,j] < CIl1)
                    
                    CIl1 = APEML[1:it,j] - SEML[1:it,j]*1.65
                    CIh1 = APEML[1:it,j] + SEML[1:it,j]*1.65
                    RES[1,5] = mean(TMP[1:it,j] > CIh1 | TMP[1:it,j] < CIl1)
                    
                    RES[1,6] = mean(SEML[1:it,j])/sd(APEML[1:it,j])
                    

                    RES[2,1] = mean(APEJK[1:it,j]/TMP[1:it,j], na.rm=T)
                    RES[2,2] = median(APEJK[1:it,j]/TMP[1:it,j], na.rm=T)
                    RES[2,3] = sd(APEJK[1:it,j])
                    
                    CIl1 = APEJK[1:it,j] - SEJK[1:it,j]*1.96
                    CIh1 = APEJK[1:it,j] + SEJK[1:it,j]*1.96
                    RES[2,4] = mean(TMP[1:it,j] > CIh1 | TMP[1:it,j] < CIl1)
                    
                    CIl1 = APEJK[1:it,j] - SEJK[1:it,j]*1.65
                    CIh1 = APEJK[1:it,j] + SEJK[1:it,j]*1.65
                    RES[2,5] = mean(TMP[1:it,j] > CIh1 | TMP[1:it,j] < CIl1)
                    
                    RES[2,6] = mean(SEJK[1:it,j])/sd(APEJK[1:it,j])
                    
                    
                    RES[3,1] = mean(APEBC[1:it,j]/TMP[1:it,j], na.rm=T)
                    RES[3,2] = median(APEBC[1:it,j]/TMP[1:it,j], na.rm=T)
                    RES[3,3] = sd(APEBC[1:it,j])
                    
                    CIl1 = APEBC[1:it,j] - SEBC[1:it,j]*1.96
                    CIh1 = APEBC[1:it,j] + SEBC[1:it,j]*1.96
                    RES[3,4] = mean(TMP[1:it,j] > CIh1 | TMP[1:it,j] < CIl1)
                    
                    CIl1 = APEBC[1:it,j] - SEBC[1:it,j]*1.65
                    CIh1 = APEBC[1:it,j] + SEBC[1:it,j]*1.65
                    RES[3,5] = mean(TMP[1:it,j] > CIh1 | TMP[1:it,j] < CIl1)
                    
                    RES[3,6] = mean(SEBC[1:it,j])/sd(APEBC[1:it,j])
                    
                    
                    RES[4,1] = mean(APECMLF[1:it,j]/TMP[1:it,j], na.rm=T)
                    RES[4,2] = median(APECMLF[1:it,j]/TMP[1:it,j], na.rm=T)
                    RES[4,3] = sd(APECMLF[1:it,j])
                    
                    CIl1 = APECMLF[1:it,j] - SECMLF[1:it,j]*1.96
                    CIh1 = APECMLF[1:it,j] + SECMLF[1:it,j]*1.96
                    RES[4,4] = mean(TMP[1:it,j] > CIh1 | TMP[1:it,j] < CIl1)
                    
                    CIl1 = APECMLF[1:it,j] - SECMLF[1:it,j]*1.65
                    CIh1 = APECMLF[1:it,j] + SECMLF[1:it,j]*1.65
                    RES[4,5] = mean(TMP[1:it,j] > CIh1 | TMP[1:it,j] < CIl1)

                    RES[4,6] = mean(SECMLF[1:it,j])/sd(APECMLF[1:it,j])


                    
                    cat("APE\n")
                    print(round(RES, digits=3), na.print = "" , quote = FALSE)
                    
                    
                }
            }
            
            
            
        }
        save.image(filename)
        
    }
}
