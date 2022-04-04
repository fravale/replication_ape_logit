## --------------------------------------------------------------------
## Replication Package for "Conditional inference and bias reduction for partial
## effects estimation of fixed-effects logit models"

## Francesco Bartolucci, Claudia Pigini and Francesco Valentini
##--------------------------------------------------------------------

jackknife_logit <- function(id, y, X, jk = TRUE, dyn = FALSE){

    ## Perform MLE & jackknife correction in fixed-effects logit
    ## HAHN & NEWEY (2004): static logit
    ## DHAENE & JOCHMANS (2015): dynamic logit (half panel)
    ## (coefficients + apes)

    ## boolean jk = FALSE does not perform jackknife correction
    ## (for recursion)
    ## if dyn=TRUE, only jackknife correction and ape standard errors are different
    ## If dynamic, remember to include lagged y in X
    
    ## WARNING: balanced panels only!!!
    
    label = unique(id); n = length(label) 
    X = as.matrix(X); k = ncol(X)
    Tv = rep(0, n)
    for (i in 1:n) Tv[i] = sum(id == label[i])
    TT = max(Tv)
    
    ## Prepare data
    Xal = diag(n)%x%rep(1,TT)
    indc = indr = NULL
    for(i in 1:n){
        y_i = y[id==i]
        sui = sum(y_i)
        if(sui==0 | sui==TT){
            indc = c(indc, i)
            indr = c(indr, rep(1,Tv[i]))
        }else{
            indr = c(indr, rep(0,Tv[i]))
        }
    }
    Xal = Xal[!indr, ]
    if(!is.null(indc)) Xal = Xal[,-indc]
    NT = sum(Tv); NTq = length(id[!indr])
    
     ## MAXIMUM LIKELIHOOD ESTIMATION
    y1 = y[!indr]; X1 = X[!indr,]
    X1 = cbind(X1, Xal); npar = ncol(X1)
    data = as.data.frame(cbind(y1,X1))
    mod = glm(y1 ~ . -1, data = data, family = binomial(link="logit"))
    be = mod$coefficients[1:k]; al = mod$coefficients[(k+1):npar]    

    ## MLE variance
    AL = rep(NA,NT)
    scb = matrix(0,n,k)
    sca = matrix(0,n,ncol(Xal)); j = 0
    for(i in 1:n){
        il = label[i]
        y_i = y[id == il]; sui = sum(y_i)
        if(sui>0 & sui<Tv[i]){
            j = j+1
            x_i = as.matrix(X[id == il,])
            xb = x_i%*%be + al[j]
            F_i = as.vector(exp(xb)/(1 + exp(xb)))
            scb[i,] = t(y_i - F_i)%*%x_i
            sca[i,j] = sum(y_i - F_i)
            AL[id == i] = al[j]
        }
    }
    scv = cbind(scb,sca)
    J = summary(mod)$cov.scaled
    Vr = J%*%(t(scv)%*%scv)%*%J
    se_be = sqrt(diag(as.matrix(Vr[1:k,1:k])))

   ## AVERAGE PARTIAL EFFECTS
    ape = rep(0,k); APE = matrix(0,n,k)
    for(j in 1:k){
        x_j = X[,j]
        if(sum((x_j)^2) == sum(x_j)){
            eta = AL + as.matrix(X[,-j])%*%be[-j]
            F0 = exp(eta)/(1 + exp(eta))
            F1 = exp(eta + be[j])/(1 + exp(eta + be[j]))
            pe = F1 - F0
        }else{
            eta = AL + X%*%be
            F = exp(eta)/(1+exp(eta))
            pe = F*(1-F)*be[j]
        }
        pe[is.na(pe)] = 0
        ape[j] = mean(pe) #,na.rm=TRUE)
        for(i in 1:n) APE[i,j] = mean(pe[id==i])
    }

    ## MLE variance for APE  (Delta method)
    Vape = rep(0,k)
    for(j in 1:k){
        x_j = X[,j]
        if(sum((x_j)^2) == sum(x_j)){
            X2 = as.matrix(X[,-j]); X2 = as.matrix(X2[!indr,])
            eta = AL[!indr] + X2%*%be[-j]
            F0 = exp(eta)/(1 + exp(eta))
            F1 = exp(eta + be[j])/(1 + exp(eta + be[j]))
            R = matrix(0,npar,1)
            v1 = (1/(NT^2))*NTq*var(F1 - F0)                
            f1 = F1*(1-F1); f0 = F0*(1-F0)
            r1 = as.matrix(sum(f1))
            R[j] = r1
            if(k>1){
                r2 = colSums(f1%*%matrix(1,1,(k-1))*X2)-
                    colSums(f0%*%matrix(1,1,(k-1))*X2)
                r2 = as.matrix(r2)
                if(j>1 & j<k){
                    R[1:(j-1)] = r2[1:(j-1)]
                    R[(j+1):k] = r2[(j+1):k]
                }else{
                    if(j==1 & k>1) R[2:k] = r2[2:k]
                    if(j==k & k>1) R[1:(k-1)] = r2[1:(k-1)]                }
            }
            r3 = t(Xal)%*%(F1 - F0); r3 = as.matrix(r3)
            R[(k+1):npar] = r3
            R = R/NT
            v2 = t(R) %*% Vr %*% R
            Vape[j] = v1 + v2
        }else{
            th = c(be, al)
            xb = X1%*%th
            F =  exp(xb)/(1+exp(xb))
            f =  F*(1-F)
            R = matrix(0,npar,1)
            v1 = (1/(NT^2))*NTq*cov(f*be[j])
            g = F*(1-F)*(1-2*F)
            X2 = as.matrix(X1[,1:k])
            r1 = sum(f) + (t(g)%*%X2[,j])*be[j]
            R[j] = r1
            if(k>1){
                r2 = as.matrix(t(g)%*%X2)
                if(j>1 & j<k){
                    R[1:(j-1)] = r2[1:(j-1)]
                    R[(j+1):k] = r2[(j+1):k]
                }else{
                    if(j==1 & k>1) R[2:k] = r2[2:k]
                    if(j==k & k>1) R[1:(k-1)] = r2[1:(k-1)]
                }
            }
            r3 = (t(Xal)%*%g)*be[j]
            R[(k+1):npar] = r3
            R = R/NT
            v2 = t(R) %*% Vr %*% R
            Vape[j] = v1 + v2
        }
    }
    se_ape = sqrt(Vape)
    
    ## JACKKNIFE CORRECTION FOR BETA AND APE
    be_jk = ape_jk = APE_jk = NULL
    if(jk){
        if(dyn){
            bej = apej = matrix(0,k,2)
            APEJ = array(0,c(n,k,2))
            half = c(rep(1,(TT/2)), rep(2,(TT/2)))
            Time = rep(1,n)%x%half
            G = 2
        }else{
            bej = apej = matrix(0,k,TT)
            APEJ = array(0,c(n,k,TT))
            Time = rep(1,n)%x%(1:TT)
            G = TT
        }
        
        for(t in 1:G){
            yj = y[Time!=t]
            Xj = X[Time!=t,]
            idj = id[Time!=t]
            modj = jackknife_logit(idj,yj,Xj,jk=FALSE,dyn=FALSE)
            bej[,t] = modj$be
            apej[,t] = modj$ape         
            APEJ[,,t] = modj$APE
        }

        be_jk = G*be - (G-1)*rowMeans(bej)
        ape_jk = G*ape - (G-1)*rowMeans(apej)
        APE_jk = G*APE - (G-1)*apply(APEJ,1,mean)

                                                   
    }

    ## CORRECTION FOR ALPHA
    ## Warning: (as a(be_jk))
    ## here jk means only that alpha is a function of be_jk
    ## NO CORRECTION for the incidental parameters problem for alpha

    if(jk){
        al_jk = rep(NA,n)
        AL_jk = rep(NA,NT)
        for (i in 1:n){
            il = label[i]; y_i = y[id == il]; sui = sum(y_i)
            if(sui>0 & sui<Tv[i]){
                x_i = as.matrix(X[id == il,])
                int = x_i %*% be_jk
                al_i = 0
                q_i = exp(int)
                q_i = q_i/(1 + q_i)
                lk1 = as.vector(y_i %*% log(q_i) + (1 - y_i) %*% 
                                log(1 - q_i))
                lk1o = -Inf
                while (abs(lk1 - lk1o) > 10^-6) {
                    lk1o = lk1                        
                    dal = sum(y_i - q_i)/sum(q_i * (1 - q_i))
                    mdal = abs(dal)
                    if (mdal > 0.5) dal = dal/mdal * 0.5
                    al_i = al_i + dal
                    q_i = exp(al_i + int)
                    q_i = q_i/(1 + q_i)
                    lk1 = as.vector(y_i %*% log(q_i) + (1 - y_i) %*% 
                                    log(1 - q_i))
                }
                AL_jk[id == i] = al_i #NT
                al_jk[i] = al_i #n
            }
        }
    }


    
    ## MLE variance for APE using Jackknife corrected estimates
    ## (Delta method)
    ## ONLY CROSS_SECTIONAL VARIATION FOR DHAENE & JOCHMANS
    se_ape_jk = NULL
    if(jk){
        Vape = rep(0,k)
        for(j in 1:k){
            x_j = X[,j]
            if(sum((x_j)^2) == sum(x_j)){
                X2 = as.matrix(X[,-j]); X2 = as.matrix(X2[!indr,])
                eta = AL_jk[!indr] + X2%*%be_jk[-j]
                F0 = exp(eta)/(1 + exp(eta))
                F1 = exp(eta + be_jk[j])/(1 + exp(eta + be_jk[j]))
                R = matrix(0,npar,1)
                v1 = (1/(NT^2))*NTq*var(F1 - F0)                    
                f1 = F1*(1-F1); f0 = F0*(1-F0)
                r1 = as.matrix(sum(f1))
                R[j] = r1
                r2 = colSums(f1%*%matrix(1,1,(k-1))*X2)-
                    colSums(f0%*%matrix(1,1,(k-1))*X2)
                r2 = as.matrix(r2)
                if(k>1){
                    r2 = colSums(f1%*%matrix(1,1,(k-1))*X2)-
                        colSums(f0%*%matrix(1,1,(k-1))*X2)
                    r2 = as.matrix(r2)
                    if(j>1 & j<k){
                        R[1:(j-1)] = r2[1:(j-1)]
                        R[(j+1):k] = r2[(j+1):k]
                    }else{
                        if(j==1 & k>1) R[2:k] = r2[2:k]
                        if(j==k & k>1) R[1:(k-1)] = r2[1:(k-1)]
                    }
                }
                r3 = t(Xal)%*%(F1 - F0); r3 = as.matrix(r3)
                R[(k+1):npar] = r3
                R = R/NT
                v2 = t(R) %*% Vr %*% R
                Vape[j] = v1 + v2
            }else{
                th_jk = c(be_jk, al_jk[!is.na(al_jk)])
                xb = X1%*%th_jk
                F =  exp(xb)/(1+exp(xb))
                f =  F*(1-F)
                R = matrix(0,npar,1)
                v1 = (1/(NT^2))*NTq*cov(f*be_jk[j])
                g = F*(1-F)*(1-2*F)
                X2 = as.matrix(X1[,1:k])
                r1 = sum(f) + (t(g)%*%X2[,j])*be_jk[j]
                R[j] = r1
                if(k>1){
                    r2 = as.matrix(t(g)%*%X2)
                    if(j>1 & j<k){
                        R[1:(j-1)] = r2[1:(j-1)]
                        R[(j+1):k] = r2[(j+1):k]
                    }else{
                        if(j==1 & k>1) R[2:k] = r2[2:k]
                        if(j==k & k>1) R[1:(k-1)] = r2[1:(k-1)]
                    }
                }
                r3 = (t(Xal)%*%g)*be_jk[j]
                R[(k+1):npar] = r3
                R = R/NT
                v2 = t(R) %*% Vr %*% R
                Vape[j] = v1 + v2
            }
        }
        
        se_ape_jk = sqrt(Vape)
    }

    out = list(be = be, se_be = se_be, ape = ape, se_ape = se_ape,
           be_jk = be_jk, ape_jk = ape_jk, se_ape_jk = se_ape_jk,
           APE = APE,APE_jk= APE_jk, AL=AL)
}
