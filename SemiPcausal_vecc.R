# fit semiparametric matern 

SemiPcausal_vecc<- function(Y,X,d=NULL,L=100,nugget=TRUE,nn=30,zhat = NULL,
                    iters=10000,burn=1000,update=1000,thin=1,Xgrid = NULL){
  library(fields)
  set.panel(1,2)
  ti = proc.time()
  Y = matrix(Y,nc=1)
  if(is.null(ncol(X))) X = matrix(X,nc=1)
  
  # initial param for MCMC
  tau    <- c(var(Y)/2) # total var
  betaX  <- rep(0,1)
  betaL  <- rep(0,L)
  rhozhat<- 0.1 # for nc matern
  nuzhat <- 1.5
  rhores <- 0.1
  nures  <- 0.5
  r      <- ifelse(nugget,0.5,1) # var/totalvar
  tauprecb <- 1 # precision for ICAR prior for b_l
  Q      <- rdist(1:L)==1 
  Q      <- diag(rowSums(Q))-Q # ICAR prior for b_l, coefficients of zhat
  # Q      <- diag(rowSums(Q))-0.99*Q # proper prior for b_l, coefficients of zhat
  
  # place holder
  keep.betaX <- matrix(0,iters,2)
  keep.betaL <- matrix(0,iters,L)
  keepers    <- matrix(0,iters,4)
  dev        <- rep(0,iters)
  colnames(keepers) <- c("sig","r","rho","sigb")
  
  
  # precompute a few quantities
  n       <- length(Y)
  distmat <- rdist(d)
  # compute confounder adjustment if not given
  # s        = Xgrid # assume X observed at a fine grid
  # alpha    = seq(0,30*pi,l=L)
  # if(is.null(zhat)) zhat = getzhat_gaussian_grid(rdist(d, s),m=alpha,rho=rhozhat,nu=nuzhat,X,rdist(s[1:2,])[2])
  nn_ordered <- find_ordered_nn(d,nn)
  Li       <- vecc(r,rhores,nures,d,nn_ordered) # SigmaRes = maternFun(distmat, c(r,rhores,nures,1/r-1))
  curdet   <- vecc_logdet(Li) # -determinant(SigmaRes)$mod
  
  for(iter in 1:iters){
    for(thinthin in 1:thin){
      
      # Update beta
      Yc       <- Y-zhat%*%betaL
      x.inv.x  <- vecc_qf_mat(Li,X,nn_ordered) # (t(cbind(1,X))%*%solve(SigmaRes,cbind(1,X)))
      x.inv.y  <- vecc_mult(X,Li,Yc,nn_ordered) # t(cbind(1,X))%*%solve(SigmaRes,Yc)
      VVV  <- solve(1/tau*x.inv.x + diag(0.01,1))
      MMM  <- 1/tau*x.inv.y
      betaX<- VVV%*%MMM+t(chol(VVV))%*%rnorm(1)
      # Update betazhat
      Yc     <- Y-X%*%betaX
      zhat.inv.zhat  <- vecc_qf_mat(Li,zhat,nn_ordered) #(t(zhat)%*%solve(SigmaRes,zhat))
      zhat.inv.y  <- vecc_mult(zhat,Li,Yc,nn_ordered) # (t(zhat)%*%solve(SigmaRes,Yc))
      VVV    <- solve(1/tau*zhat.inv.zhat + tauprecb*Q)
      MMM    <- 1/tau*zhat.inv.y
      M      <- VVV%*%MMM
      betaL  <- VVV%*%MMM+t(chol(VVV))%*%rnorm(L)
      
      # # setting the last one to zero. Do you need this?
      # condM  <- M[-L] + VVV[1:(L-1),L]/VVV[L,L]*(0-M[L])
      # condV  <- VVV[1:(L-1),1:(L-1)] + VVV[1:(L-1),L]%*%t(VVV[1:(L-1),L])/VVV[L,L]
      # betaL[1:(L-1)] <- condM + t(chol(condV))%*%rnorm(L-1)

      # update total variance and prior precision
      res  <- Y-X%*%betaX - zhat%*%betaL
      SSE  <- vecc_qf_vec(Li,res,nn_ordered) #t(res)%*%solve(SigmaRes,res)
      SSB  <- t(betaL)%*%Q%*%betaL
      tau  <- 1/rgamma(1,n/2+0.1,SSE/2+0.1)
      tauprecb <- rgamma(1,L/2+0.1,SSB/2+0.1)
      # VVV  <- 1/(1/tau*tauprecb*sum(Q) + 0.01)
      # MMM  <- 1/tau*tauprecb*sum(Q%*%betaL)
      # mub  <- rnorm(1,VVV*MMM,sqrt(VVV)) # Forgot what this is for
      
      block = F
      if(!block){
      # Update rhores, pior unif(0,1.5)
      MH     <- 0.25
      canrho <- exp(rnorm(1,log(rhores),MH))
      canLi  <- vecc(r,canrho,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
      candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
      canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
      
      canl <- 1/2*candet - 1/(2*tau)*canSSE + dunif(canrho, 0, 1.5, log = T)
      curl <- 1/2*curdet - 1/(2*tau)*SSE + dunif(rhores, 0, 1.5, log = T)
      R    <- canl-curl-log(rhores)+log(canrho)
      
      if(log(runif(1))<R){
        rhores   <- canrho
        Li       <- canLi # SigmaRes <- canSigmaRes
        curdet   <- candet
        SSE      <- canSSE# t(res)%*%solve(SigmaRes,res)
      }
      
      # Update r, pior unif(0,1)
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
          curdet <- candet
          SSE    <- canSSE
        }
      } 
      }else{
        #  need to check this block of code for block update
        # Update rhores
        MH     <- 0.25
        canrho <- exp(rnorm(1,log(rhores),MH))
        canr   <- pnorm(rnorm(1,qnorm(r),MH))
        canLi  <- vecc(canr,canrho,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
        candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
        canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
        
        canl <- 1/2*candet - 1/(2*tau)*canSSE + dunif(canrho, 0, 1.5, log = T)
        curl <- 1/2*curdet - 1/(2*tau)*SSE + dunif(rhores, 0, 1.5, log = T)
        R    <- canl-curl-log(rhores)+log(canrho)-dnorm(qnorm(r),log=TRUE)+dnorm(qnorm(canr),log=TRUE)
        
        if(log(runif(1))<R){
          rhores   <- canrho
          r        <- canr
          Li       <- canLi # SigmaRes <- canSigmaRes
          curdet   <- candet
          SSE      <- canSSE# t(res)%*%solve(SigmaRes,res)
        }
      }
    } # end thin
    
    if(iter%%update==0){
      plot(keep.betaX[,2],type="l")
      plot(keepers[,3],type="l")
    }
    
    
    dev[iter]         <- -2*(-n/2*log(tau) + 1/2*curdet - 1/(2*tau)*SSE)
    keep.betaX[iter,] <- betaX
    keep.betaL[iter,] <- betaL
    keepers[iter,]    <- c(tau,r,rhores,1/tauprecb)
  } # end iter
  
  #parameter means
  mn       <- colMeans(keepers[burn:iters,])
  betaX.mn <- colMeans(keep.betaX[burn:iters,])
  betaL.mn <- colMeans(keep.betaL[burn:iters,]) 
  res.mn <- Y-X%*%betaX.mn[2] - zhat%*%betaL.mn
  tau.mn <- mn[1]
  Li.mn  <- vecc(mn[2],mn[3],nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
  det.mn <- vecc_logdet(Li.mn) # -determinant(canSigmaRes)$mod 
  SSE.mn <- vecc_qf_vec(Li.mn,res.mn,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
  
  Dbar <- mean(dev[burn:iters])
  Dhat <- -2*(-n/2*log(tau.mn ) + 1/2*det.mn - 1/(2*tau.mn)*SSE.mn)
  pD   <- Dbar-Dhat
  DIC  <- Dbar + pD
  DIC  <- list(DIC=DIC,Dbar=Dbar,pD=pD)
  time = proc.time() - ti
  out <- list(betaX=keep.betaX,betaL=keep.betaL,dev=dev,keepers=keepers,DIC = DIC,time=time)
  return(out)}

