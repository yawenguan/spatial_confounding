# BMcausal_crho_vecc fit flexible bivariate matern using MCMC
# prameterization for smoothess:
# nuz = nux + 2C - C., nuxz = nux + C
# need to estimate C, C., rho, betaZ*sigmaZ, betax
# ---------------------
# constraints: (1) abs(rho) <= sqrt(nux*nuz)/nuxz (2) C. < u(0,nux+2c) to ensure nuz > 0
# - prior for spatial param (rho C C.) to impose proper constraints
# - p(rho|C C.)p(C.|C)p(C) = u(-sqrt(nux*nuz)/nuxz,sqrt(nux*nuz)/nuxz)*u(0,nux+2c)*halfcauchy(0,1e3)
# - normal(0,sd=10) for betaz*sigmaz and betax
# - invgamma(2,0.1) for tau (nugget variance)
# ---------------------
# the sign for cor and betaz*sigmaz aren't identifiable. 
# change prior for cor~u(0,sqrt(nux*nuz)/nuxz)

logit <- function(theta, a, b){log((theta-a)/(b-theta))}
logit.inv <- function(z, a, b){b-(b-a)/(1+exp(z))}

dinvgamma <- function(x, shape, scale=1, log=FALSE) {
  
  log.dens = shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) 
  - (scale/x)
  
  if(!log) {
    return(exp(log.dens))
  } else {
    return(log.dens)
  }
}
BMcausal_par_logl <- function(Y,X,d=NULL,L=1,nugget=TRUE,nn=30,invsig.x = NULL,nux,rhox,C0=0.5,
                          iters=10000,burn=1000,update=100,thin=1,dist_=rdist, parsimony = T,maxit=1e3){
  ####
  # parsimonious matern via MLE
  ####
  Y = matrix(Y,nc=1)
  if(is.null(ncol(X))) X = matrix(X,nc=1)
  distmat   <- dist_(d)
  sigx      <- maternFun(distmat,c(1,rhox,nux,0))
  invsig.x  <- solve(sigx,X)
  
  logit <- function(theta, a=0, b=1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a=0, b=1){b-(b-a)/(1+exp(z))}
  
  # initial param that are fixed
  rhoz <- rhoxz  <- rhox
  n         <- nrow(Y)
  distmat   <- dist_(d)
  nn_ordered<- find_ordered_nn(d,nn)
  
  logl.op  <- function(par){
    tau   = exp(par[1])
    r     = logit.inv(par[2])
    C     = exp(par[3])
    nuz    <- nux + 2*C
    nuxz   <- nux + C
    print(C)
    # form confounder
    zhat   <- maternFun(distmat,c(1,rhoxz,nuxz,0))%*%invsig.x # covparms: var, range,nu,nug/var
    xhat   <- cbind(X,zhat)
    # form residual var
    Li     <- vecc(r,rhoz,nuz,d,nn_ordered) 
    # get betahat profile beta out. 
    betahat<- solve(vecc_qf_mat(Li,xhat,nn_ordered), vecc_mult(xhat,Li,Y,nn_ordered))
    # get residuals
    # res    <- Y - X%*%betaX - zhat%*%betaL
    res    <- Y - xhat%*%betahat
    det    <- vecc_logdet(Li) # -determinant(SigmaRes)$mod
    SSE    <- vecc_qf_vec(Li,res,nn_ordered) #t(res)%*%solve(SigmaRes,res)
    logl   <- -n/2*log(tau) + 1/2*det - 1/(2*tau)*SSE
    -logl
  }
  
  # initial for optim
  C      <- C0
  tau    <- c(var(Y)/2) # total var
  r      <- ifelse(nugget,0.5,1) # var/totalvar
  
  par0 = c(log(tau),logit(r),log(C))
  # optimize it
  out = optim(par0,logl.op,control = list(trace = 6, maxit = maxit))
  
  # derive betahat for the estimated values
  par   = out$par   
  tau   = exp(par[1])
  r     = logit.inv(par[2])
  C     = exp(par[3])
  nuz    <- nux + 2*C
  nuxz   <- nux + C
  
  # form confounder
  zhat   <- maternFun(distmat,c(1,rhoxz,nuxz,0))%*%invsig.x # covparms: var, range,nu,nug/var
  xhat   <- cbind(X,zhat)
  # form residual var
  Li     <- vecc(r,rhoz,nuz,d,nn_ordered) 
  # get betahat
  betahat<- solve(vecc_qf_mat(Li,xhat,nn_ordered), vecc_mult(xhat,Li,Y,nn_ordered))
  betahatvar<- solve(vecc_qf_mat(Li,xhat,nn_ordered))
  
  return(list(out = out,beta = betahat,betavar=betahatvar))
}

BMcausal_crho_simple_logl <- function(Y,X,d=NULL,L=1,nugget=TRUE,nn=30,invsig.x = NULL,nux,rhox,C0=0.5,
                          iters=10000,burn=1000,update=100,thin=1,dist_=rdist, parsimony = F,maxit=1e3){
  # fit common range BM, constraint nuxz>nux, nuxz>nuz
  ## this is faster, but use marginal var of Z instead of the correct conditional var
 
  Y = matrix(Y,nc=1)
  if(is.null(ncol(X))) X = matrix(X,nc=1)
  
  logit <- function(theta, a=0, b=1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a=0, b=1){b-(b-a)/(1+exp(z))}
  
  # initial param that are fixed
  rhoz <- rhoxz  <- rhox
  n         <- nrow(Y)
  distmat   <- dist_(d)
  nn_ordered<- find_ordered_nn(d,nn)
  
  logl.op  <- function(par){
    tau   = exp(par[1])
    r     = logit.inv(par[2])
    C     = exp(par[3])
    C.    = exp(par[4])
    nuz   <- nux + 2*C - C. # nuxz + C-C.(nuxz can be <,= or > nuz)
    
    if(nuz>0){
    nuxz   <- nux + C
    print(C)
    # form confounder
    zhat   <- maternFun(distmat,c(1,rhoxz,nuxz,0))%*%invsig.x # covparms: var, range,nu,nug/var
    xhat   <- cbind(X,zhat)
    
    # this is faster, but use marginal var of Z instead of conditional var
    # form residual var
    Li     <- vecc(r,rhoz,nuz,d,nn_ordered)
    # get betahat profile beta out.
    betahat<- solve(vecc_qf_mat(Li,xhat,nn_ordered), vecc_mult(xhat,Li,Y,nn_ordered))
    
    # get residuals
    res    <- Y - xhat%*%betahat
    det    <- vecc_logdet(Li) # -determinant(SigmaRes)$mod
    SSE    <- vecc_qf_vec(Li,res,nn_ordered) #t(res)%*%solve(SigmaRes,res)
    
    logl   <- -n/2*log(tau) + 1/2*det - 1/(2*tau)*SSE
    -logl
    }else{
    Inf
    }
  }
  
  # initial for optim
  C      <- C0
  tau    <- c(var(Y)/2) # total var
  r      <- ifelse(nugget,0.5,1) # var/totalvar
  
  par0 = c(log(tau),logit(r),log(C), log(C+nux))
  # optimize it
  out = optim(par0,logl.op,control = list(trace = 6, maxit = maxit))
  
  # derive betahat for the estimated values
  par   = out$par   
  tau   = exp(par[1])
  r     = logit.inv(par[2])
  C     = exp(par[3])
  C.    = exp(par[4])
  nuz    <- C.#nux + 2*C - C.
  nuxz   <- nux + C
  
  # form confounder
  zhat   <- maternFun(distmat,c(1,rhoxz,nuxz,0))%*%invsig.x # covparms: var, range,nu,nug/var
  xhat   <- cbind(X,zhat)
  # form residual var
  Li     <- vecc(r,rhoz,nuz,d,nn_ordered) 
  # get betahat
  betahat<- solve(vecc_qf_mat(Li,xhat,nn_ordered), vecc_mult(xhat,Li,Y,nn_ordered))
  betahatvar<- solve(vecc_qf_mat(Li,xhat,nn_ordered))
  
  return(list(out = out,beta = betahat,betavar=betahatvar))
}

BMcausal_crho_logl <- function(Y,X,d=NULL,L=1,nugget=TRUE,nn=30,invsig.x = NULL,nux,rhox,C0=0.5,
                               iters=10000,burn=1000,update=100,thin=1,dist_=rdist, parsimony = F,maxit=1e3){
  ### slow, compute conditional var that depends on nuz, nuxz and rho
  # fit common range BM, MLE version of BMcausal_crho_vecc
  Y = matrix(Y,nc=1)
  if(is.null(ncol(X))) X = matrix(X,nc=1)
  
  logit <- function(theta, a=0, b=1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a=0, b=1){b-(b-a)/(1+exp(z))}
  
  # initial param that are fixed
  rhoz <- rhoxz  <- rhox
  n         <- nrow(Y)
  distmat   <- dist_(d)
  nn_ordered<- find_ordered_nn(d,nn)
  sigx      <- maternFun(distmat,c(1,rhox,nux,0))
  sigx.inv  = solve(sigx)
  
  logl.op  <- function(par){
    tau   = exp(par[1]) # this become the nug
    cor   = logit.inv(par[2]) # this become the correlation
    C     = exp(par[3])
    C.    = exp(par[4])
    betaz = par[5] # this is betaz*sigmaz
    
    nuz   <- nux + 2*C - C. # nuxz + C-C.(nuxz can be <,= or > nuz)

    if(nuz>0 & abs(cor)<=sqrt(nux*nuz)/nuxz){
      nuxz   <- nux + C
      print(C)
      # form confounder
      sigxz  <- maternFun(distmat, c(1,rhoxz,nuxz,0))
      zhat   <- sigxz%*%invsig.x # covparms: var, range,nu,nug/var
      xhat   <- cbind(X,zhat)
      
      # # this is slow, compute conditional var that depends on nuz, nuxz and rho
      # form residual var
      sigz  <- maternFun(distmat, c(1,rhoz,nuz,0))
      resvar<- betaz^2*(sigz - cor^2*sigxz%*%sigx.inv%*%sigxz) + diag(tau,n)
      
      resvar.inv <- solve(resvar)
      # get betahat profile beta out.
      betahat<- solve(t(X)%*%resvar.inv%*%X, t(X)%*%resvar.inv%*%(Y-zhat*cor*betaz))
      print(c(betahat,betaz))
      # get residuals
      res    <- Y - xhat%*%c(betahat,cor*betaz)
      det    <- -determinant(resvar)$mod
      SSE    <- t(res)%*%resvar.inv%*%res
      
      logl   <- 1/2*det - 1/2*SSE
      -logl
    }else{
      Inf
    }
  }
  
  # initial for optim
  C <- C.<- C0
  nuz    <- nux # nuxz + C-C.(nuxz can be <,= or > nuz)
  nuxz   <- nux + C
  tau    <- c(var(Y)/2)*0.1 # nug
  cor    = (sqrt(nux*nuz)/nuxz)*0.5
  betaz  = 1
  
  par0   = c(log(tau),logit(cor),log(C),log(C.),betaz)
  # optimize it
  out = optim(par0,logl.op,control = list(trace = 6, maxit = maxit))
  
  # derive betahat for the estimated values
  par   = out$par
  
  tau   = exp(par[1]) # this become the nug
  cor   = logit.inv(par[2]) # this become the correlation
  C     = exp(par[3])
  C.    = exp(par[4])
  betaz = par[5] # this is betaz*sigmaz
  
  nuz   <- nux + 2*C - C. # nuxz + C-C.(nuxz can be <,= or > nuz)
  nuxz  <- nux + C
  
  # form confounder
  sigxz  <- maternFun(distmat, c(1,rhoxz,nuxz,0))
  zhat   <- sigxz%*%invsig.x # covparms: var, range,nu,nug/var
  xhat   <- cbind(X,cor*zhat)

  # form residual var
  sigz  <- maternFun(distmat, c(1,rhoz,nuz,0))
  resvar<- betaz^2*(sigz - cor^2*sigxz%*%sigx.inv%*%sigxz) + diag(tau,n)
  resvar.inv <- solve(resvar)
  
  # get betahat 
  betahat<- solve(t(xhat)%*%resvar.inv%*%xhat, t(xhat)%*%resvar.inv%*%Y)
  betahatvar<- solve(t(xhat)%*%resvar.inv%*%xhat)
  
  # get betahatx
  betahat.uni<- solve(t(X)%*%resvar.inv%*%X, t(X)%*%resvar.inv%*%(Y-zhat*cor*betaz))
  betahatvar.uni<- solve(t(X)%*%resvar.inv%*%X)
  
  return(list(out = out,beta = betahat,betavar=betahatvar,beta.uni = betahat.uni,betavar.uni=betahatvar.uni))
}

BMcausal_crho_vecc<- function(Y,X,d=NULL,L=1,nugget=TRUE,nn=30,invsig.x = NULL,nux,rhox,theta=NULL,
                         iters=10000,burn=1000,update=100,thin=1,dist_=rdist,block=F){
  # can't use vecc to speed
  # fit flexible bivariate matern using Bayes
  library(fields)
  library(mvnfast)
  
  ti = proc.time()
  # initial param for MCMC
  if(is.null(theta)){
  C <- C.<- 0.5 # 0< C.< nux + 2c
  betax  <- coef(lm(Y~X-1))[1]
  nuz    <- nux + 2*C - C. # nuxz + C-C.(nuxz can be <,= or > nuz)
  nuxz   <- nux + C
  # tau    <- c(var(Y)/2)*0.1 # nug
  r      <- 0.5
  cor    <- (sqrt(nux*nuz)/nuxz)*0.5
  betaz  <- 3
  }else{
    betax  <- theta[1]
    nuz    <- theta[6] # nuxz + C-C.(nuxz can be <,= or > nuz)
    nuxz   <- theta[5]
    C      <- nuxz - nux
    C.     <- nux + 2*C - nuz
    r      <- theta[3]
    cor    <- theta[4]
    betaz  <- theta[2]
  }
  # initial param that are fixed
  rhoz <- rhoxz  <- rhox
  
  # place holder
  keep.betaX <- matrix(0,iters,1)
  keep.betaL <- matrix(0,iters,L)
  keepers    <- matrix(0,iters,5)
  dev        <- rep(0,iters)
  # a.C <- a.C.<- a.r <- a.tau <- a.betaz <- 0
  accep <- rep(0,6); MH <- rep(0.1,6)
  
  colnames(keepers) <- c("r","cor","nuxz","nuz","NA") # total var, ratio, resrange, xyrange
  # precompute a few quantities
  n         <- length(Y)
  distmat   <- dist_(d)
  sigx      <- maternFun(distmat,c(1,rhox,nux,0))
  invsig.x  <- solve(sigx,X)
  sigx.inv  <- solve(sigx)
  sigxz     <- maternFun(distmat, c(1,rhoxz,nuxz,0))
  zhat      <- sigxz%*%invsig.x # covparms: var, range,nu,nug/var
  zeros     <- rep(0,n)
  # form residual var
  sigz     <- maternFun(distmat, c(1,rhoz,nuz,0))
  res      <- Y - cbind(X,zhat)%*%c(betax,betaz*cor*sqrt(r))
  resvar   <- betaz*betaz*(r*sigz - r*cor*cor*sigxz%*%sigx.inv%*%sigxz + (1-r)*diag(n))
  curlogl  <- dmvn(c(res),zeros,resvar,log=TRUE) 
  
  for(iter in 1:iters){
    for(thinthin in 1:thin){
      
      # Update betaX 
      canbetax <- rnorm(1,betax,MH[6])
      canres  <- Y - cbind(X,zhat)%*%c(canbetax,betaz*cor*sqrt(r))
      
      canlogl  <- dmvn(c(canres),zeros,resvar,log=TRUE) 
      canl <- canlogl + dnorm(canbetax,0,10,log=T)
      curl <- curlogl + dnorm(betax,0,10,log=T)
      R    <- canl-curl
      
      if(log(runif(1))<R){
        betax    <- canbetax
        res      <- canres
        curlogl  <- canlogl
        accep[6] <- accep[6] +1
      }
      
      # below is gibbs update with normal prior
      # x.inv.x  <- t(X)%*%solve(resvar,X)
      # x.inv.y  <- t(X)%*%solve(resvar,(Y-zhat*betaz*cor*sqrt(r)))
      # VVV      <- solve(x.inv.x + 0.01)
      # MMM      <- x.inv.y
      # betax    <- VVV%*%MMM+t(chol(VVV))%*%rnorm(1)
      # # quantities depend on betax
      # res  <- Y - cbind(X,zhat)%*%c(betax,betaz*cor*sqrt(r))
      
      ### update betaz, exp(1) prior total variance####
      update.betaz = T
      if(update.betaz){
      canbeta <- exp(rnorm(1,log(betaz),MH[1]))
      canres  <- Y - cbind(X,zhat)%*%c(betax,canbeta*cor*sqrt(r))
      canresvar <- canbeta*canbeta*(r*sigz - r*cor*cor*sigxz%*%sigx.inv%*%sigxz + (1-r)*diag(n))
      
      canlogl  <- dmvn(c(canres),zeros,canresvar,log=TRUE) 
      canl <- canlogl + dgamma(canbeta,1,1,log=T)
      curl <- curlogl + dgamma(betaz,1,1,log=T)
      R    <- canl-curl-log(betaz)+log(canbeta) # +jacobian
      
      if(log(runif(1))<R){
        betaz    <- canbeta
        res      <- canres
        resvar   <- canresvar
        curlogl  <- canlogl
        accep[1] <- accep[1] +1
      }
      }
      ### update tau ####
      update.tau=F
      if(update.tau){ 
      # Update nug
      cantau   <- exp(rnorm(1,log(tau),MH[2]))
      canresvar<- betaz*betaz*(sigz - cor*cor*sigxz%*%sigx.inv%*%sigxz) + diag(cantau,n)
      
      canlogl <- dmvn(c(res),zeros,canresvar,log=TRUE) 
      canl    <- canlogl + dinvgamma(cantau,2,0.1,log=T)
      curl    <- curlogl + dinvgamma(tau,2,0.1,log=T)
      R       <- canl-curl-log(tau)+log(cantau)
      
      if(log(runif(1))<R){
        tau      <- cantau
        resvar   <- canresvar
        curlogl  <- canlogl
        accep[2] <- accep[2] +1
      }
      }
      ### update r ####
      update.r=T
      if(update.r){ 
        # Update nug
        canr     <- pnorm(rnorm(1,qnorm(r),MH[2]))
        canres   <- Y - cbind(X,zhat)%*%c(betax,betaz*cor*sqrt(canr))
        canresvar<- betaz*betaz*(canr*sigz - canr*cor*cor*sigxz%*%sigx.inv%*%sigxz + (1-canr)*diag(n))
        
        canlogl <- dmvn(c(canres),zeros,canresvar,log=TRUE) 
        canl    <- canlogl #unif prior cancel
        curl    <- curlogl 
        R       <- canl-curl-dnorm(qnorm(r),log=TRUE)+dnorm(qnorm(canr),log=TRUE)

        if(log(runif(1))<R){
          r        <- canr
          res      <- canres
          resvar   <- canresvar
          curlogl  <- canlogl
          accep[2] <- accep[2] +1
        }
      }
      
      ### update nuxz,nuz,cor ####
      if(!block){
      ### update C ####
      # both nuxz,nuz, zhat will change, half cauchy prior
      update.C = T
      if(update.C){
        # We have assumed rhores (rhoz) = rhox
        canC    <- exp(rnorm(1,log(C),MH[3]))
        cannuz  <- nux + 2*canC - C.
        cannuxz <- nux + canC
        if(cannuz>0){
          canbound <- sqrt(nux*cannuz)/cannuxz
          if(abs(cor)<=canbound){
            cansigxz<- maternFun(distmat, c(1,rhoxz,cannuxz,0))
            canzhat <- cansigxz%*%invsig.x # covparms: var, range,nu,nug/var
            canres  <- Y - cbind(X,canzhat)%*%c(betax,betaz*cor*sqrt(r))
            cansigz <- maternFun(distmat, c(1,rhoz,cannuz,0))
            canresvar<- betaz*betaz*(r*cansigz - r*cor*cor*cansigxz%*%sigx.inv%*%cansigxz + (1-r)*diag(n))
            
            canlogl <- dmvn(c(canres),zeros,canresvar,log=TRUE)
            canl    <- canlogl + dcauchy(canC,0,1e3,log=T)
          }else{canl <- -Inf}}else{canl <- -Inf}
        
        curl    <- curlogl + dcauchy(C,0,1e3,log=T)
        R       <- canl-curl-log(C)+log(canC)
        
        if(log(runif(1))<R){
          C        <- canC
          nuz      <- cannuz
          nuxz     <- cannuxz
          sigxz    <- cansigxz
          zhat     <- canzhat
          res      <- canres
          sigz     <- cansigz
          resvar   <- canresvar
          curlogl  <- canlogl
          accep[3] <- accep[3] +1
        }
      }
      
      # update C. #### 
      # nuz will change
      update.C.= T
      if(update.C.){
        # We have assumed rhores (rhoz) = rhox
        canC.    <- rtruncnorm(1,a=0,b=nux + 2*C,mean =C.,sd=MH[4]) #truncated normal(0, nux + 2c) 
        cannuz   <- nux + 2*C - canC.
        if(cannuz>0){
          canbound <- sqrt(nux*cannuz)/nuxz
          if(abs(cor)<=canbound){
            cansigz  <- maternFun(distmat, c(1,rhoz,cannuz,0))
            canresvar<- betaz*betaz*(r*cansigz - r*cor*cor*sigxz%*%sigx.inv%*%sigxz + (1-r)*diag(n))
            
            canlogl  <- dmvn(c(res),zeros,canresvar,log=TRUE)#unif prior cancel out
            canl     <- canlogl + log(dtruncnorm(C.,a=0,b=nux + 2*C,mean =canC.,sd=MH[4])) # asymmetric propsoal
          }else{canl <- -Inf}}else{canl <- -Inf}
        
        curl <- curlogl + log(dtruncnorm(canC.,a=0,b=nux + 2*C,mean =C.,sd=MH[4])) # asymmetric propsoal
        R    <- canl-curl

        if(log(runif(1))<R){
          C.       <- canC.
          nuz      <- cannuz
          sigz     <- cansigz
          resvar   <- canresvar
          curlogl  <- canlogl
          accep[4] <- accep[4] + 1
        }
      }
      # Update cor ####
      # uniform prior for r
      update.cor = T
      if(update.cor){
        bound    <- sqrt(nux*nuz)/nuxz
        cancor   <- rtruncnorm(1,a=-bound,b=bound,mean = cor,sd=MH[5]) #truncated normal(0, nux + 2c)
        canres   <- Y - cbind(X,zhat)%*%c(betax,betaz*cancor*sqrt(r)) 
        canresvar<- betaz*betaz*(r*sigz - r*cancor*cancor*sigxz%*%sigx.inv%*%sigxz + (1-r)*diag(n))
        
        canlogl <- dmvn(c(canres),zeros,canresvar,log=TRUE)#unif prior cancel out
        canl    <- canlogl + log(dtruncnorm(cor,a=-bound,b=bound,mean=cancor,sd=MH[5])) # asymmetric propsoal
        curl    <- curlogl + log(dtruncnorm(cancor,a=-bound,b=bound,mean=cor,sd=MH[5])) # asymmetric propsoal
        R       <- canl-curl
        
        if(log(runif(1))<R){
          cor     <- cancor
          res     <- canres
          resvar  <- canresvar
          curlogl <- canlogl
          accep[5]<- accep[5] +1
        }
      }
      }else{ 
        # block update C, C., Cor####
        # update C, both nuxz,nuz, zhat will change, half cauchy prior
          bound  <- sqrt(nux*nuz)/nuxz
          canC   <- exp(rnorm(1,log(C),MH[3]))
          canC.  <- rtruncnorm(1,a=0,b=nux + 2*canC,mean =C.,sd=MH[3]) #truncated normal(0, nux + 2c) 
          cannuz <- nux + 2*canC - canC.
          cannuxz<- nux + canC
          canbound  <- sqrt(nux*cannuz)/cannuxz
          cancor <- rtruncnorm(1,a=-canbound,b=canbound,mean = cor,sd=MH[3]/10) #truncated normal(-canbound, canbound)
          
          cansigxz<- maternFun(distmat, c(1,rhoxz,cannuxz,0))
          canzhat <- cansigxz%*%invsig.x # covparms: var, range,nu,nug/var
          canres  <- Y - cbind(X,canzhat)%*%c(betax,betaz*cancor*sqrt(r))
          cansigz <- maternFun(distmat, c(1,rhoz,cannuz,0))
          canresvar<- betaz*betaz*(r*cansigz - r*cancor*cancor*cansigxz%*%sigx.inv%*%cansigxz + (1-r)*diag(n))
          
          canlogl <- dmvn(c(canres),zeros,canresvar,log=TRUE)#unif prior cancel out
          canl    <- canlogl + dcauchy(canC,0,1e3,log=T) + 
                     dunif(canC.,0,nux + 2*canC,log=T) + dunif(cancor,-canbound,canbound,log=T)
          canl    <- canl + log(dtruncnorm(C.,a=0,b=nux+2*canC,mean=canC.,sd=MH[3]))+
                     log(dtruncnorm(cor,a=-canbound,b=canbound,mean=cancor,sd=MH[3]/10)) # asymmetric propsoal
          
          curl    <- curlogl + dcauchy(C,0,1e3,log=T) + 
                     dunif(C.,0,nux + 2*C,log=T) + dunif(cor,-bound,bound,log=T)
          curl    <- curl + log(dtruncnorm(canC.,a=0,b=nux + 2*canC,mean=C.,sd=MH[3]))+
                     log(dtruncnorm(cancor,a=-canbound,b=canbound,mean=cor,sd=MH[3]/10)) # asymmetric propsoal
          R       <- canl-curl-log(C)+log(canC) # +jacobian
          
          if(log(runif(1))<R){
            bound    <- canbound 
            cor      <- cancor
            C.       <- canC.
            C        <- canC
            nuz      <- cannuz
            nuxz     <- cannuxz
            sigxz    <- cansigxz
            zhat     <- canzhat
            res      <- canres
            sigz     <- cansigz
            resvar   <- canresvar
            curlogl  <- canlogl
            accep[3] <- accep[3] + 1
          }
        }
    } # end thin
    
    if(iter>burn &iter%%update==0){
      set.panel(3,2)
      plot(keep.betaX,type="l",main="betax")
      plot(keep.betaL,type="l",main="betaL")
      plot(keepers[,2],type="l",main="cor")
      plot(keepers[,1],type="l",main="r")
      plot(keepers[,3],type="l",main="nuxz")
      plot(keepers[,4],type="l",main="nuz")
      #keepers= c("r","cor","nuxz","nuz","NA")
    }
    
    # Adaptive tuning
    if(iter<burn & iter%%200==0){
      cat("accep. rate",accep/200,"\n")      
      for(j in 1:length(accep)){
      if(accep[j]>25){
      if(accep[j]/200 <0.2){MH[j] <- MH[j]*0.8}
      if(accep[j]/200 >0.5){MH[j] <- MH[j]*1.2}
      }
      accep[j]<-0
    }}
    
    if(iter>=burn & iter%%(iters*0.1)==0){
      cat("accep. betax",accep[6]/(iters-burn),"\n")
      cat("accep. betaz",accep[1]/(iters-burn),"\n")
      cat("accep. r",accep[2]/(iters-burn),"\n")
      cat("accep. C",accep[3]/(iters-burn),"\n")
      cat("accep. C.",accep[4]/(iters-burn),"\n")
      cat("accep. cor",accep[5]/(iters-burn),"\n")
    }
    
    dev[iter]         <- -2*curlogl
    keep.betaX[iter,] <- betax
    keep.betaL[iter,] <- betaz
    keepers[iter,]    <- c(r,cor,nuxz,nuz,rhoxz)# c("r","cor","nuxz","nuz","NA")
  } # end iter
  
  #parameter means
  mn       <- colMeans(keepers[burn:iters,])
  betaX.mn <- mean(keep.betaX[burn:iters,])
  betaL.mn <- mean(keep.betaL[burn:iters,]) 
  
  r.mn <- tau.mn   <- mn[1]; cor.mn <- mn[2];nuxz.mn<- mn[3];nuz.mn<-mn[4]
  sigxz.mn <- maternFun(distmat, c(1,rhoxz,nuxz.mn,0))
  zhat.mn  <- sigxz.mn%*%invsig.x # covparms: var, range,nu,nug/var
  res.mn   <- Y - cbind(X,zhat.mn)%*%c(betaX.mn,betaL.mn*cor.mn*sqrt(r.mn))
  sigz.mn  <- maternFun(distmat, c(1,rhoz,nuz.mn,0))
  resvar.mn<- betaL.mn*betaL.mn*(r.mn*sigz.mn - r.mn*cor.mn*cor.mn*sigxz.mn%*%sigx.inv%*%sigxz.mn + (1-r.mn)*diag(n))
  
  logl.mn  <- dmvn(c(res.mn),zeros,resvar.mn,log=TRUE)#unif prior cancel out
  
  Dbar <- mean(dev[burn:iters])
  Dhat <- -2*logl.mn
  pD   <- Dbar-Dhat
  DIC  <- Dbar + pD
  DIC  <- list(DIC=DIC,Dbar=Dbar,pD=pD)
  ti = proc.time()-ti
  time = ti[3]
  keepers = cbind(keep.betaL,keepers)
  colnames(keepers) <- c("var","r","cor","nuxz","nuz","NA")
  out <- list(betaX=keep.betaX,keepers=keepers,DIC = DIC,time = time)
  return(out)
}

BMcausal_vecc<- function(Y,X,d=NULL,L=1,nugget=TRUE,nn=30,invsig.x = NULL,nux,rhox,C0=0.5,
                    iters=10000,burn=1000,update=100,thin=1,dist_=rdist, parsimony = T){
  # mean = betax*X + cor*betaz*X
  # betaz = tot.var * sqrt(r/(1-cor*cor))
  # resid.var    <- [tot.var*tot.var]*(r*sigz +(1-r)*diag(n))
  library(fields)
  ti = proc.time()
  
  # fit parsimonious bivariate matern using Bayes
  # Bayes version of BMcausal_par_logl
  Y = matrix(Y,nc=1)
  if(is.null(ncol(X))) X = matrix(X,nc=1)
  
  # initial param for MCMC
  C      <- C0
  betax  <- 0
  nuz    <- nux + 2*C
  nuxz   <- nux + C
  tot.var<- c(var(Y)/2) # total var/ this the same as betaz in flex matern
  r      <- ifelse(nugget,0.5,1) # var/totalvar
  cor    <- sqrt(nux*nuz)/nuxz*0.7

  # initial param that are fixed
  rhoz <- rhoxz  <- rhox
  
  # place holder
  keep.betaX <- matrix(0,iters,1)
  keepers    <- matrix(0,iters,5)
  dev        <- rep(0,iters)
  accep      <- rep(0,5); MH <- rep(0.1,5)
  
  colnames(keepers) <- c("var","r","cor","nuxz","nuz") # ratio, cor, nuz
  # precompute a few quantities
  n         <- nrow(Y)
  distmat   <- dist_(d)
  sigx      <- maternFun(distmat,c(1,rhox,nux,0))
  invsig.x  <- solve(sigx,X)
  zhat      <- maternFun(distmat,c(1,rhoxz,nuxz,0))%*%invsig.x # covparms: var, range,nu,nug/var
  xhat      <- cbind(X,zhat)
  nn_ordered<- find_ordered_nn(d,nn)
  betaz     <- tot.var*sqrt(r/(1-cor*cor))
  res       <- Y - xhat%*%c(betax,betaz*cor)
  Li        <- vecc(r,rhoz,nuz,d,nn_ordered)
  SSE       <- vecc_qf_vec(Li,res,nn_ordered) #t(res)%*%solve(SigmaRes,res)
  curdet    <- vecc_logdet(Li) # -determinant(SigmaRes)$mod
  bound     <- sqrt(nux*nuz)/nuxz
  
  for(iter in 1:iters){
    for(thinthin in 1:thin){

      ### update sd=tot.var prior exp(1)####
      update.tot.var = T
      if(update.tot.var){
        cantot.var <- exp(rnorm(1,log(tot.var),MH[1]))
        canbetaz   <- cantot.var*sqrt(r/(1-cor*cor))
        canres  <- Y - cbind(X,zhat)%*%c(betax,canbetaz*cor)
        canSSE  <- vecc_qf_vec(Li,canres,nn_ordered) #t(res)%*%solve(SigmaRes,res)
        
        canl <- -n*log(cantot.var)- 1/(2*cantot.var*cantot.var)*canSSE  + dgamma(cantot.var,1,1,log=T)
        curl <- -n*log(tot.var)- 1/(2*tot.var*tot.var)*SSE  + dgamma(tot.var,1,1,log=T)
        R    <- canl-curl-log(tot.var)+log(cantot.var)
        
        if(log(runif(1))<R){
          tot.var  <- cantot.var
          betaz    <- canbetaz
          res      <- canres
          SSE      <- canSSE
          accep[1] <- accep[1] +1
        }
      } 
      # assmune common range
      # update rhoxz 
      update.rhoxz = F
      if(update.rhoxz){
        # we have assume rhoxz = rhox, no updating
        MH   <- 0.25
        cand0 <- log(rhoxz)
        candz <- rnorm(1,cand0,MH)
        canrhoxy <- logit.inv(candz,rhox,1) # (1*exp(candz)+rhox)/(exp(candz)+1)
        canzhat  <- maternFun(distmat,c(1,canrhoxy,nuxy,0))%*%invsig.x
        canres  <- Y-X%*%betaX - canzhat%*%betaL
        canSSE  <- vecc_qf_vec(Li,canres,nn_ordered)
        
        canl <- - 1/(2*tau)*canSSE
        curl <- - 1/(2*tau)*SSE
        R    <- (canl-log(rhoxy-rhox)-log(1-rhoxy)) -
          (curl-log(canrhoxy-rhox)-log(1-canrhoxy))
        
        if(log(runif(1))<R){
          cand0   <- candz
          rhoxy   <- canrhoxy
          zhat    <- canzhat
          res     <- canres
          SSE     <- canSSE# t(res)%*%solve(SigmaRes,res)
        }
      }
      
      update.rhoz = F
      if(update.rhoz){
        # We have assumed rhores (rhoz) = rhox
        MH   <- 0.25
        canrho <- exp(rnorm(1,log(rhores),MH))
        canLi  <- vecc(r,canrho,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
        candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
        canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
        curdet <- vecc_logdet(Li) # -determinant(SigmaRes)$mod
        
        canl <- 1/2*candet - 1/(2*tau)*canSSE
        curl <- 1/2*curdet - 1/(2*tau)*SSE
        R    <- canl-curl-log(rhores)+log(canrho)
        
        if(log(runif(1))<R){
          rhores   <- canrho
          Li       <- canLi # SigmaRes <- canSigmaRes
          SSE      <- canSSE# t(res)%*%solve(SigmaRes,res)
        }
      }
      
      # update C, both nuxz,nuz will change ####
      update.C = T
      if(update.C){
        # We have assumed rhores (rhoz) = rhox
        canC   <- exp(rnorm(1,log(C),MH[2]))
        cannuz <- nux + 2*canC
        cannuxz<- nux + canC
        canbound <- sqrt(nux*cannuz)/cannuxz
        
        if(abs(cor)<=canbound){
        canzhat <- maternFun(distmat,c(1,rhoxz,cannuxz,0))%*%invsig.x # covparms: var, range,nu,nug/var
        canres  <- Y - cbind(X,canzhat)%*%c(betax,betaz*cor)
        canLi   <- vecc(r,rhoz,cannuz,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,rho,cannures,1/r-1))
        canSSE  <- vecc_qf_vec(canLi,canres,nn_ordered) # t(res)%*%solve(SigmaRes,res)
        candet  <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
        
        canl <- 1/2*candet - 1/(2*tot.var*tot.var)*canSSE + dcauchy(canC,0,1e3,log=T)
        }else{canl <- -Inf}
        curl <- 1/2*curdet - 1/(2*tot.var*tot.var)*SSE + dcauchy(C,0,1e3,log=T)
        R    <- canl-curl-log(C)+log(canC)
        
        if(log(runif(1))<R){
          C        <- canC
          nuz      <- cannuz
          nuxz     <- cannuxz
          bound    <- canbound
          zhat     <- canzhat
          xhat     <- cbind(X,canzhat)
          res      <- canres
          Li       <- canLi # SigmaRes <- canSigmaRes
          SSE      <- canSSE# t(res)%*%solve(SigmaRes,res)
          curdet   <- candet
          accep[2] <- accep[2] +1
        }
      }
      
      # Update r, uniform prior for r ####
      update.r = T
      if(update.r){
        canr    <- pnorm(rnorm(1,qnorm(r),MH[3]))
        canbetaz<- tot.var*sqrt(canr/(1-cor*cor))
        canres  <- Y - cbind(X,zhat)%*%c(betax,canbetaz*cor)
        canLi   <- vecc(canr,rhoz,nuz,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,rho,cannures,1/r-1))
        canSSE  <- vecc_qf_vec(canLi,canres,nn_ordered) # t(res)%*%solve(SigmaRes,res)
        candet  <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
        
        canl <- 1/2*candet - 1/(2*tot.var*tot.var)*canSSE
        curl <- 1/2*curdet - 1/(2*tot.var*tot.var)*SSE
        R    <- canl-curl-dnorm(qnorm(r),log=TRUE)+dnorm(qnorm(canr),log=TRUE)
        
        if(log(runif(1))<R){
          r      <- canr
          betaz  <- canbetaz
          res    <- canres
          Li     <- canLi 
          SSE    <- canSSE
          curdet <- candet
          accep[3]  <- accep[3] +1
        }
      }
      
      # Update cor ####
      # uniform prior for cor
      update.cor = T
      if(update.cor){
        cancor   <- rtruncnorm(1,a=-bound,b=bound,mean = cor,sd=MH[4]) #truncated normal(0, nux + 2c)
        canbetaz <- tot.var*sqrt(r/(1-cancor*cancor))
        canres   <- Y - cbind(X,zhat)%*%c(betax,canbetaz*cancor) 
        canSSE   <- vecc_qf_vec(Li,canres,nn_ordered) # t(res)%*%solve(SigmaRes,res)
        
        canl    <- - 1/(2*tot.var*tot.var)*canSSE + 
                   log(dtruncnorm(cor,a=-bound,b=bound,mean=cancor,sd=MH[4])) # asymmetric propsoal
        curl    <- - 1/(2*tot.var*tot.var)*SSE + 
                   log(dtruncnorm(cancor,a=-bound,b=bound,mean=cor,sd=MH[4])) # asymmetric propsoal
        R       <- canl-curl
        
        if(log(runif(1))<R){
          cor     <- cancor
          betaz   <- canbetaz
          res     <- canres
          SSE     <- canSSE
          accep[4]<- accep[4] +1
        }
      }
      
      # Update betaX
      # canbetax<- rnorm(1,betax,MH[5])
      # canres  <- Y - cbind(X,zhat)%*%c(canbetax,betaz*cor)
      # canSSE  <- vecc_qf_vec(Li,canres,nn_ordered) # t(res)%*%solve(SigmaRes,res)
      # 
      # canl <- - 1/(2*tot.var*tot.var)*canSSE + dnorm(canbetax,0,10,log=T)
      # curl <- - 1/(2*tot.var*tot.var)*SSE + dnorm(betax,0,10,log=T)
      # R    <- canl-curl
      # 
      # if(log(runif(1))<R){
      #   betax    <- canbetax
      #   res      <- canres
      #   SSE      <- canSSE
      #   accep[5] <- accep[5] +1
      # }

      # udpate betax ####
      # below is gibbs update with normal prior
      Yc       <- Y - zhat*betaz*cor
      x.inv.x  <- vecc_qf_mat(Li,X,nn_ordered) # (t(cbind(1,X))%*%solve(SigmaRes,cbind(1,X)))
      x.inv.y  <- vecc_mult(X,Li,Yc,nn_ordered) # t(cbind(1,X))%*%solve(SigmaRes,Yc)
      VVV  <- solve(1/(tot.var*tot.var)*x.inv.x + 0.01)
      MMM  <- 1/(tot.var*tot.var)*x.inv.y
      betax <- VVV%*%MMM+t(chol(VVV))%*%rnorm(1)
      
      res   <- Y - cbind(X,zhat)%*%c(betax,betaz*cor)
      SSE   <- vecc_qf_vec(Li,res,nn_ordered)
    } # end thin
    
    if(iter>burn &iter%%update==0){
      set.panel(3,2)
      plot(keep.betaX,type="l",main="betax")
      plot(keepers[,1],type="l",main="tot.var(sd)")
      plot(keepers[,2],type="l",main="r")
      plot(keepers[,3],type="l",main="cor")
      plot(keepers[,4],type="l",main="nuxz")
      # c("var","r","cor","nuxz","nuz") 
    }
    
    # Adaptive tuning
    if(iter<burn & iter%%200==0){
      cat("accep. rate",accep/200,"\n")      
      for(j in 1:length(accep)){
        if(accep[j]>25){
          if(accep[j]/200 <0.2){MH[j] <- MH[j]*0.8}
          if(accep[j]/200 >0.5){MH[j] <- MH[j]*1.2}
        }
        accep[j]<-0
      }}
    
    if(iter>=burn & iter%%(iters*0.1)==0){
      cat("accep. betax",accep[5]/(iters-burn),"\n")
      cat("accep. tot.var",accep[1]/(iters-burn),"\n")
      cat("accep. C",accep[2]/(iters-burn),"\n")
      cat("accep. r",accep[3]/(iters-burn),"\n")
      cat("accep. cor",accep[4]/(iters-burn),"\n")
    }
    dev[iter]         <- -2*(1/2*curdet-n*log(tot.var)- 1/(2*tot.var*tot.var)*SSE)
    keep.betaX[iter,] <- betax
    keepers[iter,]    <- c(tot.var,r,cor,nuxz,nuz)
    # c("var","r","cor","nuxz","nuz") # ratio, cor, nuz
  } # end iter
  
  mn       <- colMeans(keepers[burn:iters,])
  betax.mn <- mean(keep.betaX[burn:iters,])
  tot.var.mn <- mn[1]; r.mn <- mn[2];cor.mn <- mn[3];nuxz.mn<- mn[3];nuz.mn<-mn[4]
  betaz.mn <- tot.var.mn*sqrt(r.mn/(1-cor.mn*cor.mn))
  zhat.mn <- maternFun(distmat,c(1,rhoxz,nuxz.mn,0))%*%invsig.x # covparms: var, range,nu,nug/var
  res.mn  <- Y - cbind(X,zhat.mn)%*%c(betax.mn,betaz.mn*cor.mn)
  Li.mn   <- vecc(r.mn,rhoz,nuz.mn,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,rho,cannures,1/r-1))
  SSE.mn  <- vecc_qf_vec(Li.mn,res.mn,nn_ordered) # t(res)%*%solve(SigmaRes,res)
  det.mn  <- vecc_logdet(Li.mn) # -determinant(canSigmaRes)$mod 
  logl.mn  <- -2*(1/2*det.mn-n*log(tot.var.mn)- 1/(2*tot.var.mn*tot.var.mn)*SSE.mn) 
  
  Dbar <- mean(dev[burn:iters])
  Dhat <- -2*logl.mn
  pD   <- Dbar-Dhat
  DIC  <- Dbar + pD
  DIC  <- list(DIC=DIC,Dbar=Dbar,pD=pD)
  ti = proc.time()-ti
  time = ti[3]
  
  pD   <- Dbar-Dhat
  DIC  <- Dbar + pD
  DIC  <- list(DIC=DIC,Dbar=Dbar,pD=pD)
  
  out <- list(betaX=keep.betaX,keepers=keepers,DIC = DIC,time =time)
  return(out)
  }