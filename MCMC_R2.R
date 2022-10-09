
# Form of the precision in the spectral domain
make.prec<- function(r,lam,omega){
   v <- 1/(1-lam+lam*omega)
   v <- length(v)*v/sum(v)
   v <- r*v+(1-r)
return(1/v)}

# Function to create the b-spline basis
mybs <- function(x,L){
  library(splines)
  n  <- length(x)
  x  <- x/max(x)
  if(L==1){B<-matrix(1,n,1)}
  if(L==2){B<-cbind(1-x,x)}
  if(L==3){B<-cbind(1-x,1-4*(x-0.5)^2,x)}
  if(L>=4){B<-bs(x,df=L,intercept=TRUE)}
return(B)}

# Plots data on a square
square <- function(Y,main=""){
   library(viridis)
   library(fields)
   n <-sqrt(length(Y))
   image.plot(1:n,1:n,matrix(Y,n,n),col=viridis(100),
              main=main,xlab="",ylab="",cex.main=1.8,cex.axis=1.5,
              axis.args = list(cex.axis = 1.5))
}

# Performs a Metropolis-Hastings step for a standard deviation
# with PCP/expenential prior

update.sig.pcp <- function(sig,Y,scale=-log(0.01)*0.31/0.50,A=NULL,MH=0.1,n=NULL,threshold=0.001){
    Y    <- as.vector(Y)
    if(is.null(n)){n  <- length(Y)}
    if(is.null(A)){Q  <- sum(Y^2)}
    if(!is.null(A)){Q <- as.vector(t(Y)%*%A%*%Y)}

    for(reps in 1:25){
     can  <- rgamma(1,sig/MH,1/MH)
     R    <- (-n*log(can) - 0.5*Q/can^2) -
             (-n*log(sig) - 0.5*Q/sig^2) +
             dexp(can/scale,1,log=TRUE)-
             dexp(sig/scale,1,log=TRUE)+
             dgamma(sig,can/MH,1/MH,log=TRUE)-
             dgamma(can,sig/MH,1/MH,log=TRUE)
     if(can>threshold & !is.na(R)){
       sig  <- ifelse(log(runif(1))<R,can,sig)
     }
    }
return(sig)}

update.sig.R2 <- function(sig,Y,scale=1,A=NULL,MH=0.1,n=NULL,threshold=0.001){
    Y    <- as.vector(Y)
    if(is.null(n)){n  <- length(Y)}
    if(is.null(A)){Q  <- sum(Y^2)}
    if(!is.null(A)){Q <- as.vector(t(Y)%*%A%*%Y)}

    sig2 <- sig^2
    for(reps in 1:25){
     can  <- rgamma(1,sig2/MH,1/MH)
     R    <- (-0.5*n*log(can) - 0.5*Q/can) -
             (-0.5*n*log(sig) - 0.5*Q/sig) -
             2*log(can+1)+2*log(sig2+1)+
             dgamma(sig2,can/MH,1/MH,log=TRUE)-
             dgamma(can,sig2/MH,1/MH,log=TRUE)
     if(can>threshold & !is.na(R)){
       sig2  <- ifelse(log(runif(1))<R,can,sig2)
     }
    }
return(sqrt(sig2))}


make.Q <- function(L,linear=TRUE){
 Q      <- as.matrix(dist(1:L))==1
 Q      <- diag(rowSums(Q))-Q
 if(linear){
   Q[1,L] <- 1/(L-1)
   Q[L,1] <- 1/(L-1)
   Q[1,1] <- 1-1/(L-1)
   Q[L,L] <- 1-1/(L-1) 
 }
return(Q)}

gen_inv <- function(Q,thresh=1/10^10){
 E    <- eigen(Q)
 D    <- ifelse(E$val>thresh,1/E$val,0)
 Qinv <- E$ve%*%diag(D)%*%t(E$ve)
return(Qinv)}



###################################################################
#
# This is the main MCMC function
# 
# Model
#   Ystar[i]  ~ N(Xstar[i]*beta[i],prec = tauY*make.prec(r,lam,d[i]))
#   beta     <- B%*%b
#   b         ~ Normal(0,prec = tauY Q)
#   tauY,tauB ~ Gamma(eps,eps)
#   r,lam     ~ Unif(0,1)
#
# Inputs
#   L      := number of columns in B
#   iters  := number MCMC iterations
#   burn   := number of MCMC iterations to be discarded as burn-in
#   thin   := Degree of thinning
#
###################################################################

CAR_AR <- function(Ystar,Xstar,Mstar,d,L=10,
                   eps=0.1,prior="R2",
                   iters=10000,burn=1000,update=iters,thin=1){
 
  tick <- proc.time()[3]

  Ystar <- as.vector(Ystar)
  Xstar <- as.vector(Xstar)
  d     <- as.vector(d)
  n     <- length(Ystar)

  # Create basis functions

  if(L==1){
    XB <- cbind(Mstar,Xstar)
    Q  <- 0*diag(2)
  }
  if(L>1){
   B  <- mybs(d,L)
   XB <- sweep(B,1,Xstar,"*")
   Q  <- make.Q(L)
   Q  <- sum((XB%*%gen_inv(Q))*XB)*Q/n
   XB <- cbind(Mstar,XB)
   Q  <- cbind(0,Q)
   Q  <- rbind(0,Q)
  }

  # Initial values

  b     <- rep(0,L+1)
  lam   <- 0.9
  r     <- 0.9
  w     <- make.prec(r,lam,d)
  sigY  <- sd(Ystar*sqrt(w))
  tauY  <- 1/sigY^1
  sigB  <- 0.1
  tauB  <- 1/sigB^1

  # Store results

  keep.b    <- matrix(0,iters/thin,L+1)
  keepers   <- matrix(0,iters/thin,5)
  dev       <- rep(0,iters/thin)
  colnames(keepers) <- c("sigY","sigB","r","lambda","beta0")
  
  # Start MCMC
  for(iter in 1:(iters/thin)){
    for(thinthin in 1:thin){

      # Update b/beta
        tXBw <- t(sweep(XB,1,w,"*"))
        VVV  <- solve(tauY*tXBw%*%XB + tauY*tauB*Q)
        MMM  <- tauY*tXBw%*%Ystar
        b    <- VVV%*%MMM + t(chol(VVV))%*%rnorm(L+1)
        res  <- Ystar-XB%*%b

      # Update precisions
        bQb  <- as.vector(t(b)%*%Q%*%b)
        tauY <- rgamma(1,n/2+(L-2)/2+eps,
                         sum(res*res*w)/2 + tauB*bQb/2 + eps)
        sigY <- 1/sqrt(tauY)
        if(L>1){
          if(prior=="PCP"){
            sigB  <- update.sig.pcp(sigB,b,A=tauY*Q,n=L-2)
            tauB  <- 1/sigB^2
          }
          if(prior=="R2"){
            sigB  <- update.sig.R2(sigB,b,A=tauY*Q,n=L-2)
            tauB  <- 1/sigB^2
          }
          if(prior=="Gamma"){
            tauB <- rgamma(1,(L-2)/2+eps,tauY*bQb/2+eps)
            sigB <- 1/sqrt(tauB)
          }
        }

      # Update lambda
        canl <- pnorm(rnorm(1,qnorm(lam),0.25))
        canw <- make.prec(r,canl,d)
        R    <- sum(dnorm(res,0,1/sqrt(tauY*canw),log=TRUE))-
                sum(dnorm(res,0,1/sqrt(tauY*w),log=TRUE))+
                dnorm(qnorm(canl),log=TRUE)-
                dnorm(qnorm(lam),log=TRUE)
        if(log(runif(1))<R){
           lam <- canl
           w   <- canw
        }

      # Update r
         canr <- pnorm(rnorm(1,qnorm(r),0.25))
         canw <- make.prec(canr,lam,d)
         R    <- sum(dnorm(res,0,1/sqrt(tauY*canw),log=TRUE))-
                 sum(dnorm(res,0,1/sqrt(tauY*w),log=TRUE))+
                 dnorm(qnorm(canr),log=TRUE)-
                 dnorm(qnorm(r),log=TRUE)
         if(log(runif(1))<R){
            r   <- canr
            w   <- canw
         }
     } # end thin

     if(iter%%update==0){if(L>1){
        plot(b,main=iter)
     }}

     dev[iter]      <- -2*sum(dnorm(res,0,1/sqrt(tauY*w),log=TRUE))
     keep.b[iter,]  <- b
     keepers[iter,] <- c(sigY,sigB,r,lam,b[1])

  } # end iter
  
  #burn mcmc
  keep.b <- keep.b[-(1:(burn/thin)),]
  keepers<- keepers[-(1:(burn/thin)),]
  dev    <- dev[-(1:(burn/thin))]
  
  mnY   <- XB%*%colMeans(keep.b)
  mn    <- colMeans(keepers)
  sig   <- mn[1]/sqrt(make.prec(mn[3],mn[4],d))
  Dbar  <- mean(dev)
  Dhat  <- -2*sum(dnorm(Ystar,mnY,sig,log=TRUE))
  pD    <- Dbar-Dhat
  DIC   <- Dbar + pD
  DIC   <- list(DIC=DIC,Dbar=Dbar,pD=pD)

  tock <- proc.time()[3]

  out <- list(b=keep.b[,-1],keepers=keepers,dev=dev,DIC=DIC,seconds=tock-tick)
return(out)}

