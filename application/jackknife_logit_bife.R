resampling  <-  function(id,y,X){

    label = unique(id)
    n = length(label)
    T = length(y)/n

    datamat = cbind(y,X)
    
    labeln = sample(label,replace=TRUE)
    id_tmp = labeln%x%rep(1,T)
    #browser()

    datamatb = matrix(0,n*T,ncol(datamat))

    ini = 1
    fin = T
    for(i in 1:n){
        datamatb[ini:fin,] = datamat[id==labeln[i],]
        ini = ini + T
        fin = (i+1)*T
    }

    yb = datamatb[,1]
    Xb = datamatb[,(-1)]

    out = list(y = yb, X = Xb, idb = id_tmp)
    return(out)
}


jackknife_logit <- function(id, y, X, jk = TRUE, dyn = FALSE, boot=FALSE){

    
    datmat = cbind(id,y,X)
    datibife = as.data.frame(datmat)

    if(boot){colnames(datibife) = c("id","y","nc02","nc35","nc617","hus_l_inc","l_lfp")}


    
    bform = as.formula("y~nc02 + nc35 + nc617 + hus_l_inc + l_lfp | id")
    
    bife_dyn <- bife(bform , data=datibife, control=list(dev_tol = 10^(-12)))

    bife_ape  <- get_APEs(bife_dyn)

    be = bife_dyn$coefficients
    se_be = sqrt(diag(solve(bife_dyn$Hessian))) 
    
    ape = bife_ape$delta

    se_ape = sqrt(diag(bife_ape$vcov)) 
    
    
    
    
    
    ## JACKKNIFE CORRECTION FOR BETA AND APE
    k = ncol(X)
    n = length(unique(id))
    TT = length(id)/length(unique(id))
    be_jk = ape_jk = APE_jk = NULL
    if(jk){
        if(dyn){
            bej = apej = matrix(0,k,2)
            #APEJ = array(0,c(n,k,2))
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
            datjmat = cbind(idj,yj,Xj)
            datibifejk = as.data.frame(datjmat)

            bjkform = as.formula("yj ~ nc02 + nc35 + nc617 + hus_l_inc + l_lfp | idj")


            if(boot){colnames(datibifejk) = c("idj","yj","nc02","nc35","nc617","hus_l_inc","l_lfp")}
           
            modj = bife(bjkform, data=datibifejk, control=list(dev_tol = 10^(-12)))
            modjape = get_APEs(modj)
            bej[,t] =  modj$coefficients
            apej[,t] = modjape$delta         
            #APEJ[,,t] = modj$APE
        }

        be_jk = G*be - (G-1)*rowMeans(bej)
        ape_jk = G*ape - (G-1)*rowMeans(apej)
        #APE_jk = G*APE - (G-1)*apply(APEJ,1,mean)

                                                   
    }


    
    out = list(be = be, se_be = se_be, ape = ape, se_ape=se_ape,
           be_jk = be_jk, ape_jk = ape_jk)
}
