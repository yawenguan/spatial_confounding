SP_vecc<- function(Y,X,d=NULL,L=1,nugget=TRUE,nn=30,zhat = NULL,
                   iters=10000,update=iters,thin=1,
                   dist_=rdist,nures=0.5){
  ti = proc.time()
  library(fields)
  n     <- length(Y)
  p     <- ncol(X)
  tau   <- c(var(Y)/2) # total var
  betaX  <- rep(0,p)
  betaL  <- rep(0,L)
  distmat   <- dist_(d)
  
  rhores <- 0.1
  r      <- ifelse(nugget,0.5,1) # var/totalvar
  
  keep.betaX <- matrix(0,iters,p)
  keep.betaL <- matrix(0,iters,L)
  keepers   <- matrix(0,iters,4)
  dev       <- rep(0,iters)
  colnames(keepers) <- c("sig","r","rho","sigb")
  nn_ordered  <- find_ordered_nn(d,nn)
  #
  Li  <- vecc(r,rhores,nures,d,nn_ordered)
  curdet <- vecc_logdet(Li) 
  # parameterization for GpGp
  # variance, range, smoothness, nugget/variance
  
  for(iter in 1:iters){
    for(thinthin in 1:thin){
      # no vecchia
      # SigmaRes = exp(-distmat/rhores)
      # Update beta
      x.inv.x  <- vecc_qf_mat(Li,X,nn_ordered) # (t(cbind(1,X))%*%solve(SigmaRes,cbind(1,X)))
      x.inv.y  <- vecc_mult(X,Li,Y,nn_ordered) # t(cbind(1,X))%*%solve(SigmaRes,Yc)
      VVV  <- solve(1/tau*x.inv.x + 0.01) # beta prior N(0,10)
      MMM  <- 1/tau*x.inv.y
      betaX<- VVV%*%MMM+t(chol(VVV))%*%rnorm(p)
      res  <- Y-X%*%betaX
      
      SSE  <- vecc_qf_vec(Li,res,nn_ordered) #t(res)%*%solve(SigmaRes,res)
      tau  <- 1/rgamma(1,n/2+0.1,SSE/2+0.1) # tau prior invg(0.1,0.1)
      
      block = F
      if(!block){
      # Update rhores
      MH     <- 0.25
      canrho <- exp(rnorm(1,log(rhores),MH))
      canLi  <- vecc(r,canrho,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
      candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
      canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
      # curdet <- vecc_logdet(Li) # -determinant(SigmaRes)$mod
      
      # unif prior[0,1.5] for rhores, with log normal proposal to get positive values. 
      canl <- 1/2*candet - 1/(2*tau)*canSSE + dunif(canrho, 0, 1.5, log = T)
      curl <- 1/2*curdet - 1/(2*tau)*SSE + dunif(rhores, 0, 1.5, log = T)
      R    <- canl-curl-log(rhores)+log(canrho)
      
      if(log(runif(1))<R){
        rhores   <- canrho
        Li       <- canLi # SigmaRes <- canSigmaRes
        SSE      <- canSSE# t(res)%*%solve(SigmaRes,res)
        curdet   <- candet
      }
      
      # Update ratio r, prior unif[0,1]
      if(nugget){
        MH   <- 0.25
        canr <- pnorm(rnorm(1,qnorm(r),MH))
        canLi  <- vecc(canr,rhores,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(canr,canrho,nures,1/canr-1))
        candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
        canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
        
        canl <- 1/2*candet - 1/(2*tau)*canSSE
        curl <- 1/2*curdet - 1/(2*tau)*SSE
        
        R    <- canl-curl-dnorm(qnorm(r),log=TRUE)+dnorm(qnorm(canr),log=TRUE)
        
        if(log(runif(1))<R){
          r      <- canr
          Li     <- canLi 
          SSE    <- canSSE
          curdet <- candet
        }
      } 
      }else{
        # Update rhores
        MH     <- 0.25
        canrho <- exp(rnorm(1,log(rhores),MH))
        canr   <- pnorm(rnorm(1,qnorm(r),MH))
        canLi  <- vecc(canr,canrho,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
        candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
        canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
        # curdet <- vecc_logdet(Li) # -determinant(SigmaRes)$mod
        
        canl <- 1/2*candet - 1/(2*tau)*canSSE + dunif(canrho, 0, 1.5, log = T)
        curl <- 1/2*curdet - 1/(2*tau)*SSE + dunif(rhores, 0, 1.5, log = T)
        R    <- canl-curl-log(rhores)+log(canrho)-dnorm(qnorm(r),log=TRUE)+dnorm(qnorm(canr),log=TRUE)
        
        if(log(runif(1))<R){
          rhores   <- canrho
          r        <- canr
          Li       <- canLi # SigmaRes <- canSigmaRes
          SSE      <- canSSE# t(res)%*%solve(SigmaRes,res)
          curdet   <- candet
        }
      }
    } # end thin
    
    if(iter%%update==0){
      set.panel(2,2)
      plot(keep.betaX[,p],main="betaX",type="l")
      plot(keepers[,1],main="Var",type="l")
      plot(keepers[,2],main="Ratio",type="l")
      plot(keepers[,3],main="Range",type="l")
    }
    
    # dev[iter]        <- -2*sum(dnorm(res,0,1/sqrt(1/tau*w),log=TRUE))
    keep.betaX[iter,] <- betaX
    keep.betaL[iter,] <- betaL #empty
    keepers[iter,]   <- c(tau,r,rhores,NA)
  } # end iter
  
  # mn    <- colMeans(keepers[burn:iters,])
  # sig   <- mn[1]/sqrt(make.prec.gp(mn[2],d))
  # Dbar  <- mean(dev[burn:iters])
  # if(L==1){
  #   b    <- mean(keep.beta[burn:iters])
  #   Dhat <- -2*sum(dnorm(Ystar,Xstar*b,sig,log=TRUE))
  # }
  # if(L>1){
  #   b    <- colMeans(keep.beta[burn:iters,])
  #   Dhat <- -2*sum(dnorm(Ystar,Xstar%*%b,sig,log=TRUE))
  # }
  # pD   <- Dbar-Dhat
  # DIC  <- Dbar + pD
  # DIC  <- list(DIC=DIC,Dbar=Dbar,pD=pD)
  ti = proc.time()-ti
  time = ti[3]
  out <- list(betaX=keep.betaX,betaL=keep.betaL,keepers=keepers,time=time)
  return(out)}

