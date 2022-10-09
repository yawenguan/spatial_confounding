# Plots data on a square
square <- function(Y,main=""){
  library(viridis)
  library(fields)
  n <-sqrt(length(Y))
  image.plot(1:n,1:n,matrix(Y,n,n),col=viridis(100),
             main=main,xlab="",ylab="",cex.main=1.8,cex.axis=1.5,
             axis.args = list(cex.axis = 1.5))
}

make.prec.gp<- function(r,d){
  v <- r*d+1-r
  # r=signal/total var
  return(1/v)}

spec.cor <- function(X,Y,G,D,win=50){
  n <- length(X)
  X <- t(G)%*%X
  Y <- t(G)%*%Y
  out <- X
  for(i in 1:n){
    these  <- abs(1:n-i)<win
    out[i] <- lm(Y[these]~X[these])$coef[2]
  }
  plot(D,out,xlab="Precision",ylab="Corr")
  abline(0,0)
}

dev.gp <- function(theta,X,Y,d,L){
  r   <- pnorm(theta[1])
  w   <- make.prec.gp(r,d)
  if(L==1){
    beta <- sum(w*X*Y)/sum(w*X^2)
    R    <- Y-X*beta
  }             
  if(L>1){
    tXw  <- t(sweep(X,1,w,"*"))
    beta <- solve(tXw%*%X)%*%tXw%*%Y
    R    <- Y-X%*%beta
  }
  tau <- sum(w*R*R)/n
  dev <- -2*sum(dnorm(R,0,1/sqrt(1/tau*w),log=TRUE))  
  return(dev)}

BIC.gp <- function(Ystar,Xstar,L,D){
  B     <- mybs(D,L)
  Xstar <- sweep(B,1,Xstar,"*")  
  opt   <- optim(c(1,1),dev.gp,X=Xstar,Y=Ystar,d=D,L=L)
  bic   <- opt$val + log(length(D))*(L+3)
  bic   <- opt$val + 2*(L+3)
  return(bic)}


mybs <- function(x,L){
  n  <- length(x)
  x  <- x/max(x)
  if(L==1){B<-matrix(1,n,1)}
  if(L==2){B<-cbind(1-x,x)}
  if(L==3){B<-cbind(1-x,1-4*(x-0.5)^2,x)}
  if(L>=4){B<-bs(x,df=L,intercept=TRUE)}
  return(B)}


GP <- function(Y,X,G=NULL,d=NULL,A=NULL,L=1,nugget=TRUE,lam=0,
               iters=10000,burn=1000,update=1000,thin=1){
  library(splines)
  
  n     <- length(Y)
  Ystar <- as.vector(t(G)%*%Y)
  Xstar <- t(G)%*%X
  
  B     <- mybs(d,L)
  Xstar <- sweep(B,1,Xstar,"*")
  if(L>1){ ###??
    Q     <- rdist(1:L)==1
    Q     <- diag(rowSums(Q))-0.99*Q
    Q     <- diag(L)
  }
  
  tau   <- var(Ystar)/2
  tauprecb  <- 1
  mub   <- 0
  beta  <- rep(0,L)
  r     <- ifelse(nugget,0.5,1)
  w     <- make.prec.gp(r,d)
  
  if(L==1){
    X2   <- Xstar^2
    XY   <- Xstar*Ystar
  }
  
  keep.beta <- matrix(0,iters,L)
  keepers   <- matrix(0,iters,5)
  dev       <- rep(0,iters)
  colnames(keepers) <- c("sig","r","lambda","mub","sigb")
  
  for(iter in 1:iters){
    for(thinthin in 1:thin){
      
      # Update beta
      if(L==1){
        VVV  <- 1/(1/tau*sum(w*X2) + 0.01)
        MMM  <- 1/tau*sum(w*XY)
        beta <- rnorm(1,VVV*MMM,sqrt(VVV))
        res  <- Ystar-Xstar*beta
      }             
      if(L>1){ ### why tau*taub for beta prior
        tXw   <- t(sweep(Xstar,1,w,"*"))
        VVV  <- solve(1/tau*tXw%*%Xstar + 1/tau*tauprecb*Q)
        MMM  <- 1/tau*tXw%*%Ystar+mub*1/tau*tauprecb*rowSums(Q)
        beta <- VVV%*%MMM + t(chol(VVV))%*%rnorm(L)
        res  <- Ystar-Xstar%*%beta
      }
      
      # Update total variance
      if(L==1){
        SSE <- sum(w*res^2)
        tau <- 1/rgamma(1,n/2+0.1,SSE/2+0.1)
      }
      if(L>1){ 
        SSE  <- sum(w*res^2)
        SSB  <- t(beta-mub)%*%Q%*%(beta-mub)
        tau  <- 1/rgamma(1,(n+L)/2+0.1,(SSE+tauprecb*SSB)/2+0.1)
        tauprecb <- rgamma(1,(L)/2+0.1,1/tau*SSB/2+0.1)
        
        VVV  <- 1/(1/tau*tauprecb*sum(Q) + 0.01)
        MMM  <- 1/tau*tauprecb*sum(Q%*%beta)
        mub  <- rnorm(1,VVV*MMM,sqrt(VVV))
      }
      
      # # Update lambda
      #   MH   <- 0.25
      #   canl <- pnorm(rnorm(1,qnorm(lam),MH))
      #   canw <- make.prec(r,canl,d)
      #   R    <- sum(dnorm(res,0,1/sqrt(tau*canw),log=TRUE))-
      #           sum(dnorm(res,0,1/sqrt(tau*w),log=TRUE))+
      #           dnorm(qnorm(canl),log=TRUE)-
      #           dnorm(qnorm(lam),log=TRUE)
      #   if(log(runif(1))<R){
      #      lam <- canl
      #      w   <- canw
      #   }
      
      # Update r
      if(nugget){
        MH   <- 0.25
        canr <- pnorm(rnorm(1,qnorm(r),MH))
        canw <- make.prec.gp(canr,d)
        R    <- sum(dnorm(res,0,1/sqrt(1/tau*canw),log=TRUE))-
          sum(dnorm(res,0,1/sqrt(1/tau*w),log=TRUE))+
          dnorm(qnorm(canr),log=TRUE)-
          dnorm(qnorm(r),log=TRUE)
        if(log(runif(1))<R){
          r   <- canr
          w   <- canw
        }
      } 
    } # end thin
    
    if(iter%%update==0){
    }
    
    dev[iter]        <- -2*sum(dnorm(res,0,1/sqrt(1/tau*w),log=TRUE))
    keep.beta[iter,] <- beta
    keepers[iter,]   <- c(sqrt(tau),r,lam,mub,1/sqrt(tauprecb))
    
  } # end iter
  
  mn    <- colMeans(keepers[burn:iters,])
  sig   <- mn[1]/sqrt(make.prec.gp(mn[2],d))
  Dbar  <- mean(dev[burn:iters])
  if(L==1){
    b    <- mean(keep.beta[burn:iters])
    Dhat <- -2*sum(dnorm(Ystar,Xstar*b,sig,log=TRUE))
  }
  if(L>1){
    b    <- colMeans(keep.beta[burn:iters,])
    Dhat <- -2*sum(dnorm(Ystar,Xstar%*%b,sig,log=TRUE))
  }
  pD   <- Dbar-Dhat
  DIC  <- Dbar + pD
  DIC  <- list(DIC=DIC,Dbar=Dbar,pD=pD)
  
  out <- list(beta=keep.beta,keepers=keepers,dev=dev,DIC=DIC)
  return(out)}


# functions for GP
mlegp <- function(x,y,nu=1/2,nugget = T){
  library(fields)
  library(RandomFields)
  p     <- ncol(x)
  theta <- 0.2*max(x)
  sigma <- 1
  nug   <- 0.1
  
  fun2max <- function(logparam){
    theta <- exp(logparam[1])
    sigma <- exp(logparam[2])
    nug   <- exp(logparam[3])
    V     <- sigma*RFcovmatrix(RMwhittle(nu=nu,s=theta),x=x[,1],y=x[,2])
    diag(V) <- sigma + nug
    logl  <- -1/2*determinant(V)$modulus-1/2*t(y)%*%solve(V,y)
    -logl
  }
  logparamstart = log(c(theta,sigma,nug))
  fit = optim(logparamstart,fun2max)#,control = list(maxit=1e3,trace=10))
  
  # compute output to return
  param = exp(fit$par)
  names(param) <- c("range","variance","nugget")
  return(list(theta = param))
}

# estimate spectrum using Fuentes method (plug in zeros)
sqexp_kern <- function(r, nvec){
  v1 <- as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)), 
                    seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
  v2 <- as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                    seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
  v1arr <- v1[, rep(1, nvec[2])]
  v2arr <- t(v2[, rep(1, nvec[1])])
  kern <- exp(-(v1arr^2 + v2arr^2)/r^2)
  # kern <- exp(-(v1arr^2)/r^2)
  kern <- kern#/sum(kern)
  return(kern)
}

smooth_pgram <- function(pgram, kern, smoothlog = FALSE){
  # n <- prod(dim(pgram))
  n <- length(pgram)
  if (!smoothlog) {
    smpgram <- Re(1/n * fft(fft(pgram) * fft(kern), inverse = TRUE))
  }
  else if (smoothlog) {
    sumpgram <- sum(pgram)
    logpgram <- log(pgram)
    logsmpgram <- Re(1/n * fft(fft(logpgram) * fft(kern), 
                               inverse = TRUE))
    smpgram <- exp(logsmpgram)
    smpgram <- sumpgram * smpgram/sum(smpgram)
  }
  return(smpgram)
}

sim_matern<-function(nvec, parms){
  nvec0 <- nvec
  posdef <- FALSE
  maxcount <- 5
  counter <- 0
  while (!posdef & counter < maxcount) {
    dist1 <- matrix(c(0:(nvec0[1] - 1)), nvec0[1], nvec0[2])
    dist2 <- matrix(c(0:(nvec0[2] - 1)), nvec0[1], nvec0[2], 
                    byrow = TRUE)
    distarray <- sqrt(dist1^2 + dist2^2)
    covarray <- maternFun(distarray, parms) #covariance function
    covarray_wide <- cbind(covarray, covarray[1:nvec0[1], 
                                              nvec0[2]:2])
    covarray_embed <- rbind(covarray_wide, covarray_wide[nvec0[1]:2, 
                                                         1:(2 * nvec0[2] - 1)])
    spec_embed <- Re(fft(covarray_embed))
    if (min(spec_embed) <= 0) {
      counter <- counter + 1
      nvec0 <- 2 * nvec0
    }
    else {
      posdef <- TRUE
    }
  }
  if (counter == maxcount) {
    print("Did not converge, increase maxcount")
    return(NULL)
  }
  else {
    y <- uncond_sim(spec_embed)[1:nvec[1], 1:nvec[2]]
    return(y)
  }
}

gamfn <- function(omega,c1=0.02,c2=0.08,s=0.95) s*ifelse(omega<c1,0.5,0) + s*ifelse(omega<c2,0.5,0)

maternFun <- function(distmat, covparms){
  # covparms: var, range,nu,nug/var
  n1 <- dim(distmat)[1]
  n2 <- dim(distmat)[2]
  covmat <- matrix(0, n1, n2)
  scaledist <- distmat/covparms[2]
  covparms[3] <- min(max(covparms[3], 1e-12), 20)
  normcon <- 2^(covparms[3] - 1) * gamma(covparms[3])
  besspart <- besselK(scaledist, covparms[3])
  maternpart <- 1/normcon * (scaledist)^covparms[3] * besspart
  covmat[distmat != 0] <- covparms[1] * maternpart[distmat != 0]
  diag(covmat)  <- covparms[1] * (1 + covparms[4])
  covmat[distmat == 0] <- covparms[1] * (1 + covparms[4])
  return(covmat)
}


getzhat_laplacian_grid<-function(distmat,m,rho,nu=NULL,X,gridsize){
  # nu is lambda, the scale of laplace function
  library(geoR)
  L = length(m)
  zhat = sapply(1:L, function(l){
    temp = rho[l]^2/(rho[l]^2+4*pi^2*distmat^2)*cos(2*pi*m[l]*distmat)
    temp%*%X 
  })
  zhat = gridsize*zhat
  return(zhat)}

getzhat_gaussian_grid<-function(s1,s2,m,rho,nu=NULL,X,gridsize){
  # nu is lambda, the scale of laplace function
  library(geoR)
  L = length(m)
  if(length(rho)==1) rho=rep(rho,L)
  distmat = rdist(s1,s2)
  distmatsq <- distmat^2
  distmat1 =   abs(outer(s1[,1],s2[,1],"-"))
  distmat2 =   abs(outer(s1[,2],s2[,2],"-"))
  zhat = sapply(1:L, function(l){
    temp = exp(-1/2*distmatsq*rho[l])*cos(m[l]*(distmat1+distmat2))/(2*pi)^2
    temp%*%X 
  })
  zhat = gridsize*zhat/(2*pi)^2
  return(zhat)}

getzhat_gaussian_grid1D<-function(s1,s2,m,rho,nu=NULL,X,gridsize){
  # nu is lambda, the scale of laplace function
  library(geoR)
  L = length(m)
  if(length(rho)==1) rho=rep(rho,L)
  distmat = rdist(s1,s2)
  distmatsq <- distmat^2
  zhat = sapply(1:L, function(l){
    temp =  exp(-1/2*distmatsq*rho[l])*cos(m[l]*distmat)/sqrt(2*pi)
    temp%*%X 
  })
  zhat = gridsize*zhat/(2*pi)
  return(zhat)}

getzhat_gaussian_grid2D<-function(s1,s2,m,rho,nu=NULL,X,gridsize){
  m = as.matrix(expand.grid(m,m))
  # anisotropic
  library(geoR)
  L = nrow(m)
  distmat = rdist(s1,s2)
  distmatsq <- distmat^2
  distmat1 = abs(outer(s1[,1],s2[,1],"-"))
  distmat2 = abs(outer(s1[,2],s2[,2],"-"))
  zhat = sapply(1:L, function(l){
    temp = rho*exp(-1/2*distmatsq*rho)*cos(m[l,1]*distmat1+m[l,2]*distmat2)
    temp%*%X 
  })
  zhat = gridsize*zhat/(2*pi)^2
  return(zhat)}

# plot basis
gauss_mixture <- function(m,rho,ome=NULL){
  # rho is variance of norm
  L = length(m)
  if(is.null(ome)) ome = seq(-10*pi,10*pi,l=100)
  if(length(rho)==1) rho=rep(rho,L)
  mix = sapply(1:L, function(l){
    # 2*pi*rho[l]*(0.5*dnorm(ome, m[l], sd=sqrt(rho[l])) +
    # 0.5*dnorm(ome, -m[l], sd=sqrt(rho[l])))
    (0.5*dmvn(cbind(ome,ome), rep(m[l],2), sigma=diag(rho[l],2)) +
      0.5*dmvn(cbind(ome,ome),rep(-m[l],2), sigma=diag(rho[l],2)))
    # should be bivar normal density
    # sqrt(2*pi*rho[l])*dnorm(ome, m[l], sd=sqrt(rho[l])) # take only the positive
    
  })
  # matplot(ome,mix,type="l")
  mix
}
# plot basis
gauss_mixture_2D <- function(m,rho,ome=NULL){
  library(mvnfast)
  L = nrow(m)
  if(is.null(ome)) ome = seq(-10*pi,10*pi,l=100)
  if(length(rho)==1) rho=rep(rho,L)
  mix = sapply(1:L, function(l){
    0.5*dmvn(cbind(ome[,1],ome[,2]), c(m[l,1],m[l,2]), sigma=diag(rho[l],2)) +
       0.5*dmvn(cbind(ome[,1],ome[,2]),c(-m[l,1],-m[l,2]), sigma=diag(rho[l],2))
  })
  mix
}
# basis = gauss_mixture(rho_np,m,ome)

setup_bspline_grid<-function(s1,s2,X,dx = 0.5){
# dx = 0.5 or 1,should be divisible by del
  library(splines)
  library(geoR)
  deltax   = diff(sort(unique(s1[,1]))[1:2])
  deltay   = diff(sort(unique(s1[,2]))[1:2])
  delta    = c(deltax,deltay)
  max.del.x = floor(pi/max(delta))
  del.x = seq(0,max.del.x,by = dx)
  
  distmat = rdist(s1,s2)
  foo = lapply(del.x,function(xx) besselJ(distmat*xx, nu = 0)%*%X)
  foo = do.call(cbind,foo)
  return(list(del.x=del.x,foo = foo,delta = delta,dx = dx))
}

getzhat_bspline_grid<-function(s1,s2,L,X,del,setup=NULL){
  del.x = setup$del.x
  foo = setup$foo
  dx = setup$dx
  max.del.x = max(del.x)
  Kdelta = floor(max.del.x/del)
  
  deltax   = diff(sort(unique(s2[,1]))[1:2])
  deltay   = diff(sort(unique(s2[,2]))[1:2])
  delta.x    = c(deltax,deltay)
  
  # get bspline basis fn
  del.xx = c(del.x,seq(max.del.x+dx,(Kdelta+1)*del,by = dx))
  del.xx = c(rev(-seq(dx,3*del,by = dx)), del.xx)
  
  bsbasis = bs(del.xx,knots = seq(-3*del,(Kdelta+1)*del,by=del))
  # truncate to [0,wt]
  bsbasis = bsbasis[-c(which(del.xx<0),which(del.xx>max.del.x)),]
  del.xx  = del.xx[-c(which(del.xx<0),which(del.xx>max.del.x))]
  # keep only nonzero basis
  keep = colSums(bsbasis)!=0
  bsbasis = bsbasis[,keep]
  L = ncol(bsbasis)
  
  zhat = NULL
  for(l in 1:L){
    temp1 = 0
    for(xx in 1:length(del.xx)){
      temp1 = temp1 + 2*pi*del.xx[xx]*bsbasis[xx,l]*foo[,xx]
    }
    temp1 = temp1*diff(del.x)[1]*prod(delta.x)
    zhat = cbind(zhat,temp1)
  }
  zhat = zhat/(2*pi)^2
  return(list(zhat=zhat,del.xx=del.xx,bsbasis=bsbasis))}
 
getzhat_ncmatern_1<-function(distmat,m,rho,nu,X){
  # d = d(Z(s_0),X)
  library(geoR)
  convolX = geoR::matern(distmat,rho,nu)
  zhat = sapply(m, function(x){
    foo = cos(distmat*x)*convolX# eliment wide product
    foo%*%X
  })
  zhat = zhat/length(X)
  return(zhat)}

getzhat_ncmatern_grid<-function(distmat,m,rho,nu,X,grid){
  library(geoR)
  convolX = geoR::matern(distmat,rho,nu)
  zhat = sapply(m, function(x){
    foo = cos(2*pi*distmat*x)*convolX# eliment wide product
    # changed from cos(distmat*x) to cos(2*pi*distmat*x)
    foo%*%X 
  })
  zhat = grid*zhat
  return(zhat)}

getzhat_ncmatern_varyrho<-function(distmat,m,rho,nu,X,grid){
  # d = d(Z(s_0),X)
  library(geoR)
  zhat = sapply(1:length(m), function(x){
    convolX = geoR::matern(distmat,rho[x],nu)
    foo = cos(distmat*m[x])*convolX# eliment wide product
    foo%*%X 
  })
  zhat = grid*zhat
  return(zhat)}

gp.alpha<-function(dvec,m,rho,nu){
  library(geoR)
  convolX = geoR::matern(dvec,rho,nu)
  zhat = sapply(m, function(x){
    foo = cos(dvec*x)*convolX# eliment wide product
  })
  zhat = zhat*diff(dvec)[1]
  return(zhat)}

# # visusal of getzhat 
# d = sort(runif(100,0,1))
# plot(d,cos(d*m[1])*matern(d,0.289,1.5),ylim = c(-1,1))
# for(i in 1:100)
# lines(d,cos(d*m[i])*matern(d,0.289,1.5),col = rainbow(100)[i])

# noncentral muti t spec
mnct <- function(m,delta,sigma,df){
  library(mvtnorm)
  1/2*(dmvt(m,delta,sigma,df,log=F)+dmvt(m,-delta,sigma,df,log=F))
}
# #plots the multivar noncentral t
# m = seq(0,30*pi,l=100)
# M = as.matrix(expand.grid(m,m))
# rho_np = 0.89;nu_np=1.5
# MNCT = sapply(seq(0,30*pi,l=50),function(x) mnct(cbind(diag(matrix(M[,1],100)),diag(matrix(M[,2],100))),
#       delta=rep(x, 2),sigma=rho_np^2/(2*nu_np+1)*diag(2),df=2*nu_np+1))
# matplot(m,MNCT,type="l")
# plot(m,MNCT%*%colMeans(fitnp_grid$betaL),type="l")
# lapply(seq(-pi,pi,l=10),function(x) image.plot(matrix(mnct(M,delta=rep(x, 2),
# sigma=0.2^2/2*diag(2),df=2),length(m))))

mymatern <- function(h, nu,a){
  scaledist <- h/a
  normcon <- 2^(nu - 1) * gamma(nu)
  besspart <- besselK(scaledist, nu)
  maternpart <- 1/normcon * (scaledist)^nu * besspart
  return(maternpart)
}

uncond_sim_xy<-function(nvec,parms,gamfn,c1=pi/4,c2=pi/4,plt=T,D=2,indicator = F,gaus_kern=F){
  spec_embed_list <- vector("list",length=2)
  for(k in 1:2){
    nvec00 <-nvec0 <- nvec
    posdef <- FALSE
    maxcount <- 10
    counter <- 0
    
    while (!posdef & counter < maxcount){
      dist1 <- matrix(c(0:(nvec0[1] - 1)), nvec0[1], nvec0[2])
      dist2 <- matrix(c(0:(nvec0[2] - 1)), nvec0[1], nvec0[2], 
                      byrow = TRUE)
      distarray <- sqrt(dist1^2 + dist2^2)
      covarray  <- maternFun(distarray, parms[[k]]) #covariance function
      covarray_wide <- cbind(covarray,
                             covarray[1:nvec0[1],nvec0[2]:2])
      covarray_embed<- rbind(covarray_wide,
                             covarray_wide[nvec0[1]:2,1:(2*nvec0[2]-1)])
      spec_embed <- Re(fft(covarray_embed))
      if (min(spec_embed) <= 0) {
        counter <- counter + 1
        nvec0   <- 2*nvec0
      } else {
        posdef <- TRUE
      }
    }
    
    if (counter == maxcount) {
      print("Did not converge, increase maxcount")
      return(NULL)
    }
    spec_embed_list[[k]] <- spec_embed
  }
  if(prod(dim(spec_embed_list[[1]])) != prod(dim(spec_embed_list[[2]]))){
    cat("dim doesn't match",dim(spec_embed_list[[1]]),"&",dim(spec_embed_list[[2]]),"\n")
    covarray  <- maternFun(distarray, parms[[1]]) #covariance function
    covarray_wide <- cbind(covarray,
                           covarray[1:nvec0[1],nvec0[2]:2])
    covarray_embed<- rbind(covarray_wide,
                           covarray_wide[nvec0[1]:2,1:(2*nvec0[2]-1)])
    spec_embed <- Re(fft(covarray_embed))
    spec_embed_list[[1]] <- spec_embed
    cat("final dim",dim(spec_embed_list[[1]]),"&",dim(spec_embed_list[[2]]),"\n")
  }
  nvec <- dim(spec_embed_list[[1]])
  n <- prod(nvec) 
  v1 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)),
                         seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
  v2 <-  2*pi*as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                          seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
  v1arr <- v1[, rep(1, nvec[2])]
  v2arr <- t(v2[, rep(1, nvec[1])])
  omega <- sqrt(v1arr^2 + v2arr^2)
  # myspec = matern_spec(omega,parms[3],parms[2])
  # myspec/sum(myspec) ~ spec_embed/sum(spec_embed)
  if(indicator){gamma     <- gamfn(omega,c1=c1,c2=c2)}else{
    if(gaus_kern){
      gamma     <- dnorm(omega,sd=c1)/dnorm(0,sd=c1)
    }else{
      gamma     <- matern_spec(omega,nu=c1,s=c2)/max(matern_spec(omega,nu=c1,s=c2))}}
  
  x <- matrix(rnorm(n), nvec[1], nvec[2]) # simulate spectral process
  z <- matrix(rnorm(n), nvec[1], nvec[2]) # simulate spectral process
  
  fxy <- gamma*sqrt(spec_embed_list[[1]])*sqrt(spec_embed_list[[2]])
  resspec <- spec_embed_list[[2]]*(1-gamma^2)
  
  x1 <- sqrt(spec_embed_list[[1]]) * 1/sqrt(n) * fft(x)
  z1 <- sqrt(resspec) * 1/sqrt(n) * fft(z)
  z11 <- x1*gamma*sqrt(spec_embed_list[[2]])/sqrt(spec_embed_list[[1]])
  X <- (Re(1/sqrt(n) * fft(x1,inverse = TRUE)))[1:nvec00[1], 1:nvec00[2]]
  Z <- (Re(1/sqrt(n) * fft(z1,inverse = TRUE)))[1:nvec00[1], 1:nvec00[2]]
  Z1 <- (Re(1/sqrt(n) * fft(z11,inverse = TRUE)))[1:nvec00[1], 1:nvec00[2]]
  
  if(plt){set.panel(2,2); image.plot(X,main="X");
    image.plot(Z,main="Zres");image.plot(Z1,main="XZ")
    image.plot(Z1+Z,main="Z")}
  return(list(variable1=c(X),variable2=c(Z+Z1),variable3=c(Z1)))
}


# n1 = 100
# parms <- c(1,0.16*n1,1,0)
# ome = seq(0,2*pi,l=100)
# c1 = pi/3;c2 = pi/3
# plot(ome,gamfn(ome,c1=c1,c2=c2))
# plot(ome,matern_spec(ome,parms[3],parms[2]/n1,d=2),
#      col = gamfn(ome,c1=c1,c2=c2)*2+1)
# 
# set.seed(1)
# foo = uncond_sim_xy(c(n1,n1),parms,gamfn,c1=c1,c2=c2,plt = T)
# set.panel(2,2);par(mar=c(3,3,0,5))
# square(foo$variable1);square(foo$variable2)
# square(foo$variable3)
# square(foo$variable1-foo$variable3)
# save=c()
# for(i in 1:1000){
# foo = uncond_sim_xy(c(50,50),parms,gamfn,c1=100000,c2=100000)
# save = c(save,cor(foo$variable1,foo$variable2))
# }
# hist(save)

# # 1d case
# n = 100
# omega = seq(0,pi,l=n)
# s = seq(0,100,l=n)
# foo = matern_spec(omega,0.5,1/1.8,d=1)
# plot(omega,(foo),type="l")
# set.seed(0)
# std = rnorm(n)
# z = fft(std)/sqrt(n)
# 
# z = rnorm(n)+sqrt(-1+0i)*rnorm(n)
# plot(density(Re(z)),col="red")
# lines(density(Im(z)),col="blue")
# lines(density(std))
# 
# z = z*sqrt(foo)
# realy = Re(unlist(lapply(s,function(h) sum(exp(sqrt(-1+0i)*omega*h)*z))))
# plot(s,realy,type='l')
# foo3 =RFsimulate(RMwhittle(nu=0.5,s=1/1.8),x=s)
# plot(s,foo3@data$variable1,type="l")
# X = Re(fft(z,T))
# plot(s,X,col=2,type="l")
# 
# 
# plot((imat[1,1:n1]))
# plot(seq(0,2*pi,l=100),myimat,type="l")

grid.prep <- function(W, M, N, ext = 2) {
  cell.width <- diff(W$xrange)/M
  cell.height <- diff(W$yrange)/N
  
  mgrid <- seq(W$xrange[1], W$xrange[2], by = cell.width)
  ngrid <- seq(W$yrange[1], W$yrange[2], by = cell.height)
  mcens <- (mgrid + 0.5 * cell.width)[-(M + 1)]
  ncens <- (ngrid + 0.5 * cell.height)[-(N + 1)]
  
  if (ext <= 1) 
    mgrid.ext <- ngrid.ext <- mcens.ext <- ncens.ext <- M.ext <- N.ext <- NULL else {
      M.ext <- ext * M
      N.ext <- ext * N
      mgrid.ext <- seq(W$xrange[1], W$xrange[2] + (ext - 1) * diff(W$xrange), by = cell.width)
      ngrid.ext <- seq(W$yrange[1], W$yrange[2] + (ext - 1) * diff(W$yrange), by = cell.height)
      mcens.ext <- (mgrid.ext + 0.5 * cell.width)[-(M.ext + 1)]
      ncens.ext <- (ngrid.ext + 0.5 * cell.height)[-(N.ext + 1)]
    }
  
  return(list(M = M, N = N, mgrid = mgrid, ngrid = ngrid, mcens = mcens, ncens = ncens, 
              cell.width = cell.width, cell.height = cell.height, M.ext = M.ext, N.ext = N.ext, 
              mgrid.ext = mgrid.ext, ngrid.ext = ngrid.ext, mcens.ext = mcens.ext, ncens.ext = ncens.ext))
}
covariance.prep <- function(gp,variance,r,wrap,...){
  if(!wrap){
    cent <- expand.grid(gp$mcens,gp$ncens)
    mmat <- matrix(rep(cent[,1],gp$M*gp$N),gp$M*gp$N,gp$M*gp$N)
    nmat <- matrix(rep(cent[,2],gp$M*gp$N),gp$M*gp$N,gp$M*gp$N)
    D <- sqrt((mmat-t(mmat))^2+(nmat-t(nmat))^2)
    covmat <- variance*r(D,...)
  } else {
    Rx <- gp$M.ext*gp$cell.width
    Ry <- gp$N.ext*gp$cell.height
    m.abs.diff.row1 <- abs(gp$mcens.ext[1]-gp$mcens.ext)
    m.diff.row1 <- pmin(m.abs.diff.row1,Rx-m.abs.diff.row1)
    n.abs.diff.row1 <- abs(gp$ncens.ext[1]-gp$ncens.ext)
    n.diff.row1 <- pmin(n.abs.diff.row1,Ry-n.abs.diff.row1)
    cent.ext.row1 <- expand.grid(m.diff.row1,n.diff.row1)
    D.ext.row1 <- matrix(sqrt(cent.ext.row1[,1]^2+cent.ext.row1[,2]^2),gp$M.ext,gp$N.ext)
    C.tilde <- variance*r(D.ext.row1,...)
    return(D.ext.row1)
  }
}

uncond_sim_xy1<-function(nvec,parms,gamfn){
  m1 = nvec[1]
  m2 = nvec[2]
  n <- prod(nvec)
  # make omega
  w1 <- 2*pi*(seq(0,2*m1-1,1))/(2*m1)
  w2 <- 2*pi*(seq(0,2*m2-1,1))/(2*m2)
  w1 <- matrix(w1,2*m1,2*m2,byrow=FALSE)
  w2 <- matrix(w2,2*m1,2*m2,byrow=TRUE)
  omega  <- sqrt(w1^2+w2^2)
  gamma <- gamfn(omega,pi/4,pi/4)
  
  nvec_doub <- dim(omega)
  n_doub    <- prod(nvec_doub)
  x <- matrix(rnorm(n_doub), nvec_doub[1], nvec_doub[2]) # simulate spectral process
  z <- matrix(rnorm(n_doub), nvec_doub[1], nvec_doub[2]) # simulate spectral process
  spec <- matern_spec(omega,parms[3],s=parms[2])
  x1 <- sqrt(spec) * 1/sqrt(n_doub) * fft(x)
  z1 <- x1*gamma + sqrt(spec) * 1/sqrt(n_doub) * fft(z)
  
  X <- Re(1/sqrt(n_doub) * fft(x1,inverse = TRUE))
  Z <- Re(1/sqrt(n_doub) * fft(z1,inverse = TRUE))
  square(X);square(Z)
  return(list(variable1=c(X),variable2=c(Z)))
}

getspec <- function(X,m1,m2,K,BS=T){
  cmaq = matrix(X,m1,m2)
  
  # add taper
  taper_prop = 0.1
  taper1 <- spec.taper(x = rep(1,m1), p = taper_prop)
  taper2 <- spec.taper(x = rep(1,m2), p = taper_prop)
  taper <- outer(taper1,taper2)
  cmaq  <- cmaq*taper
  scaling <- sum(taper^2)
  
  #Fourier frequencies
  w1 <- 2*pi*(seq(0,m1-1,1))/m1
  w2 <- 2*pi*(seq(0,m2-1,1))/m2
  w1 <- matrix(w1,m1,m2,byrow=FALSE)
  w2 <- matrix(w2,m1,m2,byrow=TRUE)
  w  <- sqrt(w1^2+w2^2)
  
  #Accounting for aliasing
  wbar  <- sqrt(ifelse(w1>0,1,0)*(2*pi-w1)^2+
                  ifelse(w2>0,1,0)*(2*pi-w2)^2)
  delta <- ifelse(w<=wbar,w,wbar)
  #fft
  z <- fft(cmaq,inverse=TRUE)
  
  Xtilde <- matrix(NA,nrow=m1*m2,ncol=K)
  #Create the constructed covariates using bernstein
  if(!BS){
    for(k in 1:K){
      W           <- dbinom(k-1, K-1,delta/(2*pi))
      x           <- fft(z*W)
      
      xtilde <- Re(x)/scaling
      Xtilde[,k] = xtilde
    }
  }else{
    #  use bspline
    foo = bs(delta/(2*pi),K,intercept = T)
    for(k in 1:K){
      W           <- matrix(foo[,k],m1)
      x           <- fft(z*W)
      
      xtilde <- Re(x)/scaling
      Xtilde[,k] = xtilde
    }
  }
  return(Xtilde)
}

getspec1 <- function(X,m1,m2,K,BS=T){
  
  cmaq = matrix(X,m1,m2)
  nvec = c(m1,m2)
  
  # add taper
  taper_prop = 0.1
  taper1 <- spec.taper(x = rep(1,m1), p = taper_prop)
  taper2 <- spec.taper(x = rep(1,m2), p = taper_prop)
  taper <- outer(taper1,taper2)
  cmaq  <- cmaq*taper
  scaling <- sum(taper^2)
  #Fourier frequencies
  w1 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)), 
                    seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
  w2 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                    seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
  w1 <- w1[, rep(1, nvec[2])]
  w2 <- t(w2[, rep(1, nvec[1])])
  w  <- sqrt(w1^2+w2^2)
  
  #Accounting for aliasing
  wbar  <- sqrt(ifelse(w1>0,1,0)*(2*pi-w1)^2+
                  ifelse(w2>0,1,0)*(2*pi-w2)^2)
  delta <- ifelse(w<=wbar,w,wbar)
  #fft
  
  z <- fft(cmaq,inverse=TRUE)
  
  Xtilde <- matrix(NA,nrow=m1*m2,ncol=K)
  #Create the constructed covariates using bernstein
  if(!BS){
    for(k in 1:K){
      W           <- dbinom(k-1, K-1,delta/max(delta))
      x           <- fft(z*W)
      
      xtilde <- Re(x)/scaling
      Xtilde[,k] = xtilde
    }
  }else{
    #  use bspline
    foo = bs(delta,K,intercept = T)
    for(k in 1:K){
      W           <- matrix(foo[,k],m1)
      x           <- fft(z*W)
      
      xtilde <- Re(x)/scaling
      Xtilde[,k] = xtilde
    }
  }
  return(Xtilde)
}

getspec_lap <- function(X,m1,m2,K){
  
  cmaq = matrix(X,m1,m2)
  nvec = c(m1,m2)
  
  # add taper
  taper_prop = 0.1
  taper1 <- spec.taper(x = rep(1,m1), p = taper_prop)
  taper2 <- spec.taper(x = rep(1,m2), p = taper_prop)
  taper <- outer(taper1,taper2)
  cmaq  <- cmaq*taper
  scaling <- sum(taper^2)
  #Fourier frequencies
  w1 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)), 
                         seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
  w2 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                         seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
  w1 <- w1[, rep(1, nvec[2])]
  w2 <- t(w2[, rep(1, nvec[1])])
  delta <- w  <- sqrt(w1^2+w2^2)
  
  z <- fft(cmaq,inverse=TRUE)
  
  Xtilde <- matrix(NA,nrow=m1*m2,ncol=K)
  #  use lapalce
  
  m = seq(0,2*pi,l=K)
  rho_np = rep(2,l=K)
  
  lap_mixture = sapply(1:K, function(l){
    0.5*laplace(rho_np[l],m[l],delta) +
      0.5*laplace(rho_np[l],m[l],-delta)
  })
  
  for(k in 1:K){
    W           <- matrix(lap_mixture[,k],m1)
    x           <- fft(z*W)
    
    xtilde <- Re(x)/scaling
    Xtilde[,k] = xtilde
  }
    
  return(Xtilde)
}

smoothfre <- function(X,m1,m2,K,BS=T){
  X=X-mean(X)
  cmaq = matrix(X,m1,m2)
  nvec = c(m1,m2)
  # add taper
  taper_prop = 0.1
  taper1 <- spec.taper(x = rep(1,m1), p = taper_prop)
  taper2 <- spec.taper(x = rep(1,m2), p = taper_prop)
  taper <- outer(taper1,taper2)
  cmaq  <- cmaq*taper
  scaling <- sum(taper^2)
  #Fourier frequencies
  w1 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)), 
                         seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
  w2 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                         seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
  w1 <- w1[, rep(1, nvec[2])]
  w2 <- t(w2[, rep(1, nvec[1])])
  w  <- sqrt(w1^2+w2^2)
  
  #Accounting for aliasing
  wbar  <- sqrt(ifelse(w1>0,1,0)*(2*pi-w1)^2+
                  ifelse(w2>0,1,0)*(2*pi-w2)^2)
  delta <- ifelse(w<=wbar,w,wbar)
  #fft
  
  z <- fft(cmaq,inverse=TRUE)
  
  #Create the constructed covariates using bernstein
  if(!BS){
    for(k in 1:K){
      W           <- dbinom(k-1, K-1,delta/max(delta))
      Xtilde[,k] = W
    }
  }else{
    #  use bspline
    Xtilde = bs(delta,K,intercept = T)
  }
  return(list(freq=c(z),basis=Xtilde))
}

getspecxy1 <- function(X,Y,m1,m2,K,BS=T,DS=F,Kern=T){
  cmaqX = matrix(X,m1,m2)
  cmaqY = matrix(Y,m1,m2)
  nvec = c(m1,m2)
  #Fourier frequencies
  w1 <- as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)), 
                    seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
  w2 <- as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                    seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
  w1 <- w1[, rep(1, nvec[2])]
  w2 <- t(w2[, rep(1, nvec[1])])
  w  <- sqrt(w1^2+w2^2)
  
  #Accounting for aliasing
  delta <- w 
  #fft
  x <- fft(cmaqX,inverse=TRUE)
  y <- fft(cmaqY,inverse=TRUE)
  
  Xtilde <- matrix(NA,nrow=m1*m2,ncol=K)
  #  use bspline
  foo = bs(delta,K,intercept = T)
  for(k in 1:K){
    W           <- matrix(foo[,k],m1)
    fo          <- fft(sqrt(y*x)*W)
    xtilde <- Re(fo)/(m1*m2)
    Xtilde[,k] = xtilde
  }
  return(Xtilde)
}

getspec2 <- function(X,m1,m2,K,BS=T,DS=F,Kern=F){
  library(splines)
  cmaq = matrix(X,m1,m2)
  #Fourier frequencies
  w1 <- 2*pi*(seq(0,m1-1,1))/m1
  w2 <- 2*pi*(seq(0,m2-1,1))/m2
  w1 <- matrix(w1,m1,m2,byrow=FALSE)
  w2 <- matrix(w2,m1,m2,byrow=TRUE)
  w  <- sqrt(w1^2+w2^2)
  
  #Accounting for aliasing
  wbar  <- sqrt(ifelse(w1>0,1,0)*(2*pi-w1)^2+
                  ifelse(w2>0,1,0)*(2*pi-w2)^2)
  delta <- ifelse(w<=wbar,w,wbar)
  #fft
  z <- fft(cmaq,inverse=TRUE)
  Xtilde <- matrix(NA,nrow=m1*m2,ncol=K)
  #Create the constructed covariates using bernstein
  if(!BS){
    for(k in 1:K){
      W           <- dbinom(k-1, K-1,delta/(2*pi))
      x           <- fft(z*W)
      
      xtilde <- Re(x)/(m1*m2)
      Xtilde[,k] = xtilde
    }
  }else{
    #  use bspline
    foo = bs(delta/(2*pi),K,intercept = T)
    for(k in 1:K){
      W           <- matrix(foo[,k],m1)
      x           <- fft(z*W)
      
      xtilde <- Re(x)/(m1*m2)
      Xtilde[,k] = xtilde
    }
  }
  if(DS){
    freqseg = seq(0,2*pi,by=2*pi/K)
    for(k in 1:(length(freqseg)-1)){
      W           <- dbinom(k-1, K-1,delta/(2*pi))
      W[which(delta < freqseg[k] | delta >= freqseg[k+1])]=0
      W[which(delta >= freqseg[k] & delta < freqseg[k+1])]=1
      x           <- fft(z*W)
      xtilde <- Re(x)/(m1*m2)
      Xtilde[,k] = xtilde
    }
    if(Kern){
      freqseg = seq(0,2*pi,by=2*pi/K)
      for(k in 1:(length(freqseg)-1)){
        W           <- dbinom(k-1, K-1,delta/(2*pi))
        W[which(delta < freqseg[k] | delta >= freqseg[k+1])]=0
        W[which(delta >= freqseg[k] & delta < freqseg[k+1])]=1
        x           <- fft(z*W)
        xtilde <- Re(x)/(m1*m2)
        Xtilde[,k] = xtilde
      }
    }
    
  }
  return(Xtilde)
}

getspecxy <- function(X,Y,m1,m2,K,BS=T,DS=F,Kern=T){
  cmaqX = matrix(X,m1,m2)
  cmaqY = matrix(Y,m1,m2)
  #Fourier frequencies
  w1 <- 2*pi*(seq(0,m1-1,1))/m1
  w2 <- 2*pi*(seq(0,m2-1,1))/m2
  w1 <- matrix(w1,m1,m2,byrow=FALSE)
  w2 <- matrix(w2,m1,m2,byrow=TRUE)
  w  <- sqrt(w1^2+w2^2)
  
  #Accounting for aliasing
  wbar  <- sqrt(ifelse(w1>0,1,0)*(2*pi-w1)^2+
                  ifelse(w2>0,1,0)*(2*pi-w2)^2)
  delta <- ifelse(w<=wbar,w,wbar)
  #fft
  x <- fft(cmaqX,inverse=TRUE)
  y <- fft(cmaqY,inverse=TRUE)
  
  Xtilde <- matrix(NA,nrow=m1*m2,ncol=K)
  #  use bspline
  foo = bs(delta/(2*pi),K,intercept = T)
  for(k in 1:K){
    W           <- matrix(foo[,k],m1)
    fo          <- fft(sqrt(y*x)*W)
    xtilde <- Re(fo)/(m1*m2)
    Xtilde[,k] = xtilde
  }
  return(Xtilde)
}

get_alasedspec<-function(nvec, parms){
  nvec0 <- nvec
  posdef <- FALSE
  maxcount <- 5
  counter <- 0
  while (!posdef & counter < maxcount) {
    dist1 <- matrix(c(0:(nvec0[1] - 1)), nvec0[1], nvec0[2])
    dist2 <- matrix(c(0:(nvec0[2] - 1)), nvec0[1], nvec0[2], 
                    byrow = TRUE)
    distarray <- sqrt(dist1^2 + dist2^2)
    covarray <- maternFun(distarray, parms) #covariance function
    covarray_wide <- cbind(covarray, covarray[1:nvec0[1], 
                                              nvec0[2]:2])
    covarray_embed <- rbind(covarray_wide, covarray_wide[nvec0[1]:2, 
                                                         1:(2 * nvec0[2] - 1)])
    spec_embed <- Re(fft(covarray_embed))
    if (min(spec_embed) <= 0) {
      counter <- counter + 1
      nvec0 <- 2 * nvec0
    }
    else {
      posdef <- TRUE
    }
  }
  if (counter == maxcount) {
    print("Did not converge, increase maxcount")
    return(NULL)
  }
  else {
    y <- uncond_sim(spec_embed)[1:nvec[1], 1:nvec[2]]
    return(y)
  }
}

getspec.cor <- function(X,Z,m1,m2){
  nvec = c(m1,m2)
  # add taper to boundary values
  taper_prop = 0.1
  taper1 <- spec.taper(x = rep(1,m1), p = taper_prop)
  taper2 <- spec.taper(x = rep(1,m2), p = taper_prop)
  taper <- outer(taper1,taper2)
  scaling <- sum(taper^2)
  #Fourier frequencies
  w1 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[1] - 1)/2)), 
                         seq(floor((nvec[1] - 1)/2), 1, by = -1)))/nvec[1]
  w2 <- 2*pi*as.matrix(c(seq(0, ceiling((nvec[2] - 1)/2)),
                         seq(floor((nvec[2] - 1)/2), 1, by = -1)))/nvec[2]
  w1 <- w1[, rep(1, nvec[2])]
  w2 <- t(w2[, rep(1, nvec[1])])
  w  <- sqrt(w1^2+w2^2)
  delta = w
  
  # #Accounting for aliasing # not sure how to use this for the new code
  # wbar  <- sqrt(ifelse(w1>0,1,0)*(pi-w1)^2+
  #                 ifelse(w2>0,1,0)*(pi-w2)^2)
  # delta <- ifelse(w<=wbar,w,wbar)
  
  #fft
  
  Rexspec <- matrix(NA,nr=nrow(X),nc=ncol(X))
  Rezspec <- matrix(NA,nr=nrow(Z),nc=ncol(Z))
  for(ind in 1:ncol(X)){
  Xmat = matrix(X[,ind],m1,m2)
  Zmat = matrix(Z[,ind],m1,m2)
  Xmat  <- Xmat*taper
  Zmat  <- Zmat*taper
  
  xspec <- fft(Xmat,inverse=TRUE)
  zspec <- fft(Zmat,inverse=TRUE)
  
  Rexspec[,ind] = c(Re(xspec))
  Rezspec[,ind] = c(Re(zspec))
  }
  v = seq(min(delta),max(delta),l=100)
  i <- findInterval(delta,v)
  
  corV <-c()
  for(ii in 1:length(v)){
    id = which(i == ii)
    corV = c(corV, cor(c(Rexspec[id,]),c(Rezspec[id,])))
  }
  
  MarX <-c() 
  for(ii in 1:length(v)){
    id = which(i == ii)
    MarX = c(MarX, var(c(Rexspec[id,])))
  }
  
  plot(v,corV,type="l")
  marX = matern_spec(v,nu1,s1)
  dat = cbind(c(delta),c(Re(xspec)),c(Re(zspec)),i)
  colnames(dat) <- c("delta","x","z","i")
  # foo=by(dat,dat[,"i"], function(ii) cor(ii[,"x"],ii[,"z"]))
  return(Xtilde)
}


lik = function(param,s1,c1,RinvX=NULL,nu1){
  nu12=nu1
  # rho    = 0.3*exp(param[3])/(exp(param[3])+1)
  # s2      =  (1*exp(param[3])+s1)/(exp(param[3])+1)
  s12      = (1*exp(param[4])+s1)/(exp(param[4])+1)
  tau2   = exp(param[2])
  sigmaz2= 0#exp(param[3])
  # cat("rho",rho,"c",s12/s1,"tau2",tau2,"\n")
  # cat("sigmaz2",sigmaz2,"c",s12/s1,"tau2",tau2,"\n")
  
  if(is.null(RinvX)){
    R = RFcovmatrix(RMwhittle(nu=nu1,s=s1),x=s[,1],y=s[,2]) #exp(-D/rhox)
    RinvX = solve(R,X)
  }
  # 2^(1-nu)*gamma(nu)^{-1}*(h/s)^nu*besselK(h/s,nu)
  SXZ = RFcovmatrix(RMwhittle(nu=nu12,s=s12),x=s[,1],y=s[,2])
  # SZ = RFcovmatrix(RMwhittle(nu=nu1,s=s2),x=s[,1],y=s[,2])
  # w = sigmaz2*SZ - rho^2*sigmaz2*t(SXZ)%*%solve(R,SXZ)
  # diag(w) = diag(w)+tau2
  w = diag(sigmaz2+tau2,n)
  # Zhat = rho*sigmaz2/c1*SXZ%*%solve(R,X)
  Zhat = SXZ%*%RinvX
  # Yhat = solve(R,Y)
  
  XM = cbind(X,Zhat)
  # foo=0
  # foo = tryCatch(bhat <- solve(t(XM)%*%solve(w,XM),t(XM)%*%solve(w,Y)),
  #                  error = function(e) {
  #                    print("___failed__________");
  #                    print(param);
  #                    return(param);
  #                  })
  # 
  # print(bhat)
  
  bhat0 = solve(t(XM)%*%XM,t(XM)%*%Y)
  # print(bhat0)
  e = Y-XM%*%bhat0
  nlogl = n*log(diag(w)[1]) + t(e)%*%e/diag(w)[1]
  
  
  # e = Y-XM%*%bhat
  # if(length(foo)>2) {foo <<- param ;break}
  # nlogl = determinant(w)$mod + t(e)%*%solve(w,e)
  return(nlogl)
}	

lik_w = function(param,s1,c1,totalvar,s2,nu1,nu2){
  # r: sigmaz2/sigmaz2+tau2
  r      = (1*exp(param[2]))/(exp(param[2])+1)
  rho    = 0#exp(param[3])/(exp(param[3])+1)
  # s2      =  (1*exp(param[3])+s1)/(exp(param[3])+1)
  s12    = (1*exp(param[4])+max(c(s1,s2)))/(exp(param[4])+1)
  sigmaz2= r*totalvar
  
  R = RFcovmatrix(RMwhittle(nu=nu1,s=s1),x=s[,1],y=s[,2]) #exp(-D/rhox)
  SXZ = RFcovmatrix(RMwhittle(nu=nu12,s=s12),x=s[,1],y=s[,2])
  SZ = RFcovmatrix(RMwhittle(nu=nu1,s=s2),x=s[,1],y=s[,2])
  w = sigmaz2*SZ - rho^2*sigmaz2*t(SXZ)%*%solve(R,SXZ)
  diag(w) = totalvar
  Zhat = sigmaz2/c1*SXZ%*%solve(R,X)
  # Zhat = SXZ%*%solve(R,X)
  # Yhat = solve(R,Y)
  XM = cbind(X,Zhat)
  # foo=0
  # foo = tryCatch(bhat <- solve(t(XM)%*%solve(w,XM),t(XM)%*%solve(w,Y)),
  #                  error = function(e) {
  #                    print("___failed__________");
  #                    print(param);
  #                    return(param);
  #                  })
  # 
  # print(bhat)
  
  # bhat0 = solve(t(XM)%*%XM,t(XM)%*%Y)
  # e = Y-XM%*%bhat0
  # nlogl = n*log(diag(w)[1]) + t(e)%*%e/diag(w)[1]
  
  bhat <- solve(t(XM)%*%solve(w,XM),t(XM)%*%solve(w,Y))
  e = Y-XM%*%bhat
  # if(length(foo)>2) {foo <<- param ;break}
  nlogl = determinant(w)$mod + t(e)%*%solve(w,e)
  cat("sigmaz2",sigmaz2,"rho",rho,"s12",s12,"\n")
  cat("bhat",bhat,"\n")
  return(nlogl)
}	

lik_np = function(param,s1,c1){
  tau2   = exp(param[1])
  R = RFcovmatrix(RMwhittle(nu=nu1,s=s1),x=s[,1],y=s[,2]) #exp(-D/rhox)
  SXZ2 = RFcovmatrix(RMwhittle(nu=1,s=2*s1),x=s[,1],y=s[,2])
  SXZ3 = RFcovmatrix(RMwhittle(nu=1,s=4*s1),x=s[,1],y=s[,2])
  SXZ5 = RFcovmatrix(RMwhittle(nu=5,s=2*s1),x=s[,1],y=s[,2])
  SXZ6 = RFcovmatrix(RMwhittle(nu=5,s=4*s1),x=s[,1],y=s[,2])
  
  w = diag(tau2,n)
  RinvX = solve(R,X)
  Zhat = cbind(
    SXZ2%*%RinvX,
    SXZ3%*%RinvX,
    SXZ5%*%RinvX,
    SXZ6%*%RinvX)
  
  XM = cbind(X,Zhat)
  bhat0 = solve(t(XM)%*%XM,t(XM)%*%Y)
  print(bhat0)
  e = Y-XM%*%bhat0
  nlogl = n*log(diag(w)[1]) + t(e)%*%e/diag(w)[1]
  
  return(nlogl)
}	

bhatfn = function(param,s1,c1,nu1,RinvX=NULL){
  nu12 = nu1
  s12      = (1*exp(param[4])+s1)/(exp(param[4])+1)
  tau2   = exp(param[2])
  sigmaz2= 0#exp(param[3])
  
  if(is.null(RinvX)){
    R = RFcovmatrix(RMwhittle(nu=nu1,s=s1),x=s[,1],y=s[,2]) #exp(-D/rhox)
    RinvX = solve(R,X)
  }
  
  SXZ = RFcovmatrix(RMwhittle(nu=nu12,s=s12),x=s[,1],y=s[,2])
  w = diag(sigmaz2+tau2,n)
  Zhat = SXZ%*%RinvX
  XM = cbind(X,Zhat)
  bhat0 = solve(t(XM)%*%XM,t(XM)%*%Y)
  bhatsd = sqrt(diag(tau2*solve(t(XM)%*%XM)))
  cat("s12",s12,"tau2",tau2,"\n")
  return(list(bhat = bhat0,bsd = bhatsd))
}

bhatfn_w = function(param,s1,c1,totalvar,s2){
  nu1=nu12=1.5
  r      = (1*exp(param[2]))/(exp(param[2])+1)
  rho    = 0#exp(param[3])/(exp(param[3])+1)
  s12    = (1*exp(param[4])+max(c(s1,s2)))/(exp(param[4])+1)
  sigmaz2= r*totalvar
  tau2 = (1-r)*totalvar
  
  R = RFcovmatrix(RMwhittle(nu=nu1,s=s1),x=s[,1],y=s[,2]) #exp(-D/rhox)
  SXZ = RFcovmatrix(RMwhittle(nu=nu12,s=s12),x=s[,1],y=s[,2])
  SZ = RFcovmatrix(RMwhittle(nu=nu1,s=s2),x=s[,1],y=s[,2])
  w = sigmaz2*SZ - rho^2*sigmaz2*t(SXZ)%*%solve(R,SXZ)
  diag(w) = totalvar
  Zhat =  sigmaz2/c1*SXZ%*%solve(R,X)
  XM = cbind(X,Zhat)
  
  bhat <- solve(t(XM)%*%solve(w,XM),t(XM)%*%solve(w,Y))
  bhatsd = sqrt(diag(tau2*solve(t(XM)%*%XM)))
  cat("s12",s12,"s2",s2,"sigmaz2",sigmaz2,"\n")
  return(list(bhat = bhat,bsd = bhatsd))
}

library(GpGp)

vecc <- function(r,rho,nu,S,nn){
  Li  <- vecchia_Linv(c(r,rho,nu,1/r-1),
                      "matern_isotropic",
                      S,nn)
  return(Li)}

# log(det(inverse_cov))
vecc_logdet <- function(Li){
  2*sum(log(Li[,1]))
}

# t(y)%*%inverse_cov%*%y for vector y
vecc_qf_vec <- function(Li,y,nn){
  sum(Linv_mult(Li,y,nn)^2)
}

# t(y)%*%inverse_cov%*%y for matrix y
vecc_qf_mat <- function(Li,y,nn){
  Ly <- NULL
  for(j in 1:ncol(y)){
    Ly <- cbind(Ly,Linv_mult(Li,y[,j],nn))
  }
  out <- t(Ly)%*%Ly
  return(out)}


# t(x)%*%inverse_cov%*%y
vecc_mult <- function(x,Li,y,nn){
  if(is.vector(x)){
    Lx <- Linv_mult(Li,x,nn)
  }
  if(is.matrix(x)){
    Lx <- NULL
    for(j in 1:ncol(x)){
      Lx <- cbind(Lx,Linv_mult(Li,x[,j],nn))
    }
  }
  if(is.vector(y)){
    Ly <- Linv_mult(Li,y,nn)
  }
  if(is.matrix(y)){
    Ly <- NULL
    for(j in 1:ncol(y)){
      Ly <- cbind(Ly,Linv_mult(Li,y[,j],nn))
    }
  }
  out <- t(Lx)%*%Ly
  return(out)}