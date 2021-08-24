## --------------------------------------------------------------------
## Replication Package for "Conditional inference and bias reduction for partial
## effects estimation of fixed-effects logit models"

## Francesco Bartolucci, Claudia Pigini and Francesco Valentini
##--------------------------------------------------------------------

ape_correct  <- function(pid,yv,X,be,al,dyn=FALSE,ape){
    label = unique(pid)
    n = length(label)
    X = as.matrix(X)
    k = ncol(X)
    
    Tv = rep(0, n)
    ind = pid
    for (i in 1:n) Tv[i] = sum(ind == label[i])
    TT = max(Tv)
    nobs = 0
    NT = sum(Tv)

    D = rep(0,k)
    SPE = rep(0,k)
    be = as.vector(be)
    for(i in 1:n){
        il = label[i]
        y_i = yv[pid == il]; sui = sum(y_i)
        if(sui>0 & sui<Tv[i]){
            x_i = as.matrix(X[pid == il,])
            for(j in 1:k){
                x_ij = x_i[,j]
                if(sum((x_ij)^2) == sum(x_ij)){
                    xb = as.matrix(x_i[,-j])%*%be[-j] + al[i]
                    F_i0 = exp(xb)/(1 + exp(xb))
                    F_i1 = exp(xb + be[j])/(1 + exp(xb + be[j]))
                    f_i0 = F_i0*(1 - F_i0); f_i1 = F_i1*(1 - F_i1)
                    g_i0 = f_i0*(1 - 2*F_i0); g_i1 = f_i1*(1 - 2*F_i1)
                    dmeff = f_i1 - f_i0
                    ddmeff = g_i1 - g_i0

                    xb = x_i%*%be + al[i]
                    F_i = as.vector(exp(xb)/(1 + exp(xb)))                
                    f_i = F_i*(1 - F_i); f_i = matrix(f_i, Tv[i],1)
                    g_i = f_i*(1 - 2*F_i);  g_i = matrix(g_i, Tv[i],1)
                }else{
                    xb = x_i%*%be + al[i]
                    F_i = as.vector(exp(xb)/(1 + exp(xb)))                
                    f_i = F_i*(1 - F_i); f_i = matrix(f_i, Tv[i],1)
                    g_i = f_i*(1 - 2*F_i);  g_i = matrix(g_i, Tv[i],1)
                    h_i = f_i*(((1-2*F_i)^2) - 2*f_i); h_i = matrix(h_i, Tv[i],1)
                    
                    dmeff = g_i*be[j]
                    ddmeff = h_i*be[j]
                    #browser()
                }
                s2 = 1/mean(f_i)
                psi_i = (y_i - F_i)*s2
            
                c = -((s2^2)/(2*Tv[i]))*sum(g_i)                    
                if(dyn) c = c - (s2/((Tv[i]-1))*sum(f_i[-1]*psi_i[-Tv[i]]))
                        
                Di = (c/Tv[i])*sum(dmeff) + (s2/(2*Tv[i]))*sum(ddmeff)
                D[j] = D[j] + Di

                Corri = Di
                
                if(dyn){
                    SPEi = ((1/(Tv[i] -1))*sum(dmeff[-1]*psi_i[-Tv[i]]))
                    SPE[j] = SPE[j] + SPEi
                    Corri = Corri + Spei
                }
            }
            
        }
    }
    
    ape_bc = ape - ((D+SPE)/NT)

    out = list(ape = ape_bc)
}

ape_uncorrect  <- function(pid,yv,X,be,AL,dyn=FALSE){
    label = unique(pid)
    n = length(label)
    X = as.matrix(X)
    k = ncol(X)
    
    Tv = rep(0, n)
    ind = pid
    for (i in 1:n) Tv[i] = sum(ind == label[i])
    TT = max(Tv)
    nobs = 0
    NT = sum(Tv)

  
   
    
    ape = rep(0,k)
    PE = matrix(0,sum(Tv),k)

    for(j in 1:k){
        x_j = X[,j]
        if(sum((x_j)^2) == sum(x_j)){
            eta = AL + as.matrix(X[,-j])%*%be[-j]
            p0 = exp(eta)/(1 + exp(eta))
            p1 = exp(eta + be[j])/(1 + exp(eta + be[j]))
            pe = p1 - p0
            PE[,j] = pe
        }else{
            eta = AL + X%*%be
            p = exp(eta)/(1+exp(eta))
            pe = p*(1-p)*be[j]
            PE[,j] = pe
        }
        pe[is.na(pe)] = 0
        ape[j] = mean(pe)
    }

    out = list(ape = ape, PE = PE)
}


cquad_bc_ape <- function(id, yv, X, be, scv, J, dyn=FALSE){

    pid = id
    r = length(pid)
    label = unique(pid)
    n = length(label)
    X = as.matrix(X)
    k = ncol(X)

    Tv = rep(0, n)
    ind = id
    for (i in 1:n) Tv[i] = sum(ind == label[i])
    TT = max(Tv)
    nobs = 0
    NT = sum(Tv)


    AL = rep(NA, length(yv))
    alp = rep(NA,n)
    for (i in 1:n) {
        il = label[i]
        y_i = yv[pid == il]
        sui = sum(y_i)
        if (sui > 0 & sui < Tv[i]) {
            nobs = nobs+Tv[i]
            x_i = as.matrix(X[pid == il, ])
            int = x_i %*% be
            al = 0
            q_i = exp(int)
            q_i = q_i/(1 + q_i)
            lk1 = as.vector(y_i %*% log(q_i) + (1 - y_i) %*% 
                            log(1 - q_i))
            lk1o = -Inf
            while (abs(lk1 - lk1o) > 10^-6) {
                lk1o = lk1

                    dal = sum(y_i - q_i)/sum(q_i * (1 - q_i))
                
                mdal = abs(dal)
                 if (mdal > 0.5) 
                     dal = dal/mdal * 0.5
                al = al + dal
                q_i = exp(al + int)
                q_i = q_i/(1 + q_i)
                lk1 = as.vector(y_i %*% log(q_i) + (1 - y_i) %*% 
                                log(1 - q_i))
              
            }
            AL[pid == il] = al
            alp[i] = al
        }
        
    }
    

    ## Compute Uncorrected APE
    ua = ape_uncorrect(pid,yv,X,be,AL,dyn)
    PE = ua$PE
    ape = ua$ape

    ## APE correction
    abc = ape_correct(pid,yv,X,be,alp,dyn=FALSE,ape)
    ape_bc = abc$ape 

    
############ Compute GMM standard errors
    
    S2 = matrix(0,n,k)
    J11 = J; j22 = 0
    
    for(i in 1:n){
        il = label[i]
        y_i = yv[pid == il]
        sui = sum(y_i)
        if (sui > 0 & sui < Tv[i]) {
            j22 = j22 + 2
            for(j in 1:k){
                pe_i = PE[pid == il, j]
                S2[i,j] = - (2/Tv[i])*sum(pe_i - ape_bc[j]) #
            }
        }
    }
    J22 = diag(k)*j22
    S = cbind(scv, S2)

    J21 = matrix(0,k,k)
    for(h in 1:k){
        be1 = be
        be1[h] = be1[h] + 10^(-4)
        
        AL1 = rep(NA, length(yv))
        alp1 = rep(NA,n)
        for (i in 1:n) {
            il = label[i]
            y_i = yv[pid == il]
            sui = sum(y_i)
            if (sui > 0 & sui < Tv[i]) {
                x_i = as.matrix(X[pid == il, ])
                int = x_i %*% be1
                al = 0
                q_i = exp(al + int)
                q_i = q_i/(1 + q_i)
                lk1 = as.vector(y_i %*% log(q_i) + (1 - y_i) %*% 
                                log(1 - q_i))
                lk1o = -Inf
                while (abs(lk1 - lk1o) > 10^-10) {
                    lk1o = lk1

                        dal = sum(y_i - q_i)/sum(q_i * (1 - q_i))
                
                 mdal = abs(dal)
                 if (mdal > 0.5) 
                     dal = dal/mdal * 0.5
                    al = al + dal
                    q_i = exp(al + int)
                    q_i = q_i/(1 + q_i)
                    lk1 = as.vector(y_i %*% log(q_i) + (1 - y_i) %*% 
                                log(1 - q_i))

                }
                AL1[pid == il] = al
                alp1[i] = al
            }
        }

        PE_1 = matrix(0,sum(Tv),k)
        ape1 = rep(0,k)
        sc2_1 = rep(0,k)
        sc2 = colSums(S2)


        ## Compute PE
        ua1 = ape_uncorrect(pid,yv,X,be1,AL1,dyn)
        PE1 = ua1$PE
        
     
        ## Numerical Differentiation
        for(j in 1:k){
            for(i in 1:n){
                il = label[i]
                y_i = yv[pid == il]
                sui = sum(y_i)
                if (sui > 0 & sui < Tv[i]) {
                    pe_i1 = PE1[pid == il, j]
                    sc2_1[j] = sc2_1[j] - (2/Tv[i])*sum(pe_i1 - ape_bc[j]) 
                }
            }
        }
        J21[,h] = (sc2_1 - sc2) *  10^4       
        
    }
    PE[is.na(PE)] = 0

    H = rbind(cbind(J11, matrix(0,k,k)), cbind(J21, J22))
    iH = solve(H)
    SS = (t(S)%*%S)
    Va =  iH %*% SS %*% t(iH)
    Va[(k+1):(2*k), (k+1):(2*k)] = Va[(k+1):(2*k), (k+1):(2*k)]
    se = sqrt(diag(Va))
    se = se[(k+1):(2*k)]

   ############################################

    export = !is.na(al)
    alex = al[export]

    
    out = list(ape = ape_bc, al = alex , PE = PE, se = se, Va = Va, nobs = nobs, AL = AL)
    
}

