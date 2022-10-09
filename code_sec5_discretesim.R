###################################################################################################
# example code to perform discrete-space simulation study
# X ~ CAR(0,sigx^2,lam), we use Leroux CAR here 
# Z|X ~ CAR(bXZ*W*X,sigz^2,lam)
# Y|X,Z ~ N(bX*X + bZ*Z, sigma^2)
# W is kernel smoothing matrix with bandwidth phi=bw
#
# Notes:
# - To save time, this script simulates two data sets for each scenario in Table 1.
#   For each simulated data set, we compare four methods:
#     the standard, parametric, semiparametric-PCP and semiparametric-R2
#   This takes about 50 minutes to complete
# 
# - The actual simulation study for Table 1 is performed on a high performance computer 
#   where each scenario has 500 simulated data sets and takes about 48 hours
###################################################################################################

library(fields)
library(splines)
library(mvtnorm)
library(eCAR)
source("MCMC_R2.R")


for(design in 1:5){
print(design)

# nsims <- 500 # use 500 replicates to produce Table 1 result
nsims <- 2 # number of replicates
n1    <- 40  # n1xn1 grid
n     <- n1^2 # n locations
lam   <- 0.95  
bX    <- 0.5
bZ    <- 0.5
sigma <- 0.25
sigx  <- sqrt(1.7)
sigz  <- sqrt(1.0)
iters <- 25000
burn  <- 5000
thin  <- 20
Ls    <- c(1,5,10,20,30,40) # try different number of basis for semi-parametric model

if(design==1){bXZ <- 0; bw=1}
if(design==2){bXZ <- 1; bw=1}
if(design==3){bXZ <- 2; bw=1}
if(design==4){bXZ <- 1; bw=2}
if(design==5){bXZ <- 2; bw=2}


# Define rook neighborhood structure 
s     <- expand.grid(1:n1,1:n1)
A     <- rdist(s)==1 # adjacency
R     <- diag(rowSums(A))-A # precision of ICAR
E     <- eigen(R) # eigen component
G     <- E$vec
D     <- E$val
rm(E,R)

# Kernel smoothing matrix
W     <- exp(-(rdist(s)/bw)^2)
W     <- sweep(W,1,rowSums(W),"/")


# To store results
b_mn      <- matrix(0,nsims,4)
b_sd <- b_mn
LDIC      <- matrix(0,nsims,2)
beta.hat1 <- beta.hat2 <- beta.hat3 <- matrix(0,n,nsims)

boot_rmse_se <- function(x,iters=100){
     y <- matrix(0,iters,ncol(x))
     for(i in 1:iters){
       samp  <- sample(1:nrow(x),replace=TRUE)
       y[i,] <- sqrt(colMeans(x[samp,]^2))
     }
return(apply(y,2,sd))}


for(sim in 1:nsims){
  print(paste(sim,"of",nsims))
  set.seed(sim*0820+design)

  w    <- 1/sqrt(1-lam+lam*D)
  X    <- sigx*G%*%rnorm(n,0,w) # simulate X from Leroux CAR model
  Z    <- bXZ*W%*%X + sigz*G%*%rnorm(n,0,w)
  Y    <- bX*X + bZ*Z + rnorm(n,0,sigma)
  
  # assume same adjacency structure and project all variables 
  # into the spectral domain using graph Fourier transform
  Xstar <- t(G)%*%X
  Ystar <- t(G)%*%Y
  Zstar <- t(G)%*%Z
  Mstar <- t(G)%*%rep(1,n)

  
  # if(sim==1){ ## plot the simulated data
  #    pdf(paste0("sim",design,".pdf"))
  #    square(Z,main="Z")
  #    square(X,main="X")
  #    square(Y,main="Y")
  #    dev.off()
  # }

  #### Standard Leroux CAR Base model ####
  fit0 <- CAR_AR(Ystar,Xstar,Mstar,d=D,L=1,iters=iters,burn=burn,thin = thin)

  #### Spectral Semiparametric ####
  # PCP model # try different L and use the one with smallest DIC
  L1 <- Ls
  if(length(Ls)>1){
   dic <- Ls
   for(ddd in 1:length(Ls)){
     dic[ddd] <- CAR_AR(Ystar,Xstar,Mstar,d=D,L=Ls[ddd],prior="PCP",iters=5000,burn=1000)$DIC$DIC
   }
   print(dic-dic[1])
   L1          <- Ls[which.min(dic)]
  }
  LDIC[sim,1] <- L1
  if(L1==1){
    fit1            <- fit0
    fit1$b          <- matrix(fit0$b,nc=1)
    beta.hat1[,sim] <- mean(fit1$b)
  }
  if(L1>1){
    fit1            <- CAR_AR(Ystar,Xstar,Mstar,d=D,L=L1,prior="PCP",iters=iters,burn=burn,thin=thin)
    beta.hat1[,sim] <- mybs(D,L=L1)%*%colMeans(fit1$b)
  }
  
  #### Spectral Semiparametric ####
  # R2 model # try different L and use the one with smallest DIC
  L2 <- Ls
  if(length(Ls)>1){
    dic <- Ls
    for(ddd in 1:length(Ls)){
       dic[ddd] <- CAR_AR(Ystar,Xstar,Mstar,d=D,L=Ls[ddd],prior="R2",iters=5000,burn=1000)$DIC$DIC
    }
    print(dic-dic[1])
    L2          <- Ls[which.min(dic)]
  }
  LDIC[sim,2] <- L2
  if(L2==1){
    fit2            <- fit0
    fit2$b          <- matrix(fit0$b,nc = 1)
    beta.hat2[,sim] <- mean(fit1$b)
  }
  if(L2>1){
    fit2            <- CAR_AR(Ystar,Xstar,Mstar,d=D,L=L2,prior="R2",iters=iters,burn=burn, thin=thin)
    beta.hat2[,sim] <- mybs(D,L=L2)%*%colMeans(fit2$b)
  }
  
  #### Spectral Parametric ####
  fit3 <- par.eCAR.Leroux(y=Y, x=X, W=A, model="Gaussian", draws=iters, burn=burn, thin=thin, verbose=FALSE)
  
  # organize the output to get beta(omega) estimate
  poster <- fit3$posterior_draws
  numdraws = length(poster$beta)
  beta.func <- matrix(NA, nrow=numdraws, ncol=length(D))
  
  for(kk in 1:numdraws){
    beta.func[kk,] <- poster$beta[kk] + 
      poster$rho[kk]*sqrt(poster$sig2z[kk]/poster$sig2x[kk]) * 
      sqrt((1-poster$lamx[kk] + poster$lamx[kk]*rev(D))/
             (1-poster$lamz[kk] + poster$lamz[kk]*rev(D)))
  }
  beta.hat3[,sim] <- colMeans(beta.func)
  
  boxplot(fit1$b,outline=FALSE)
  abline(bX,0,col=3,lwd=2)
  abline(quantile(fit0$b,0.5),0,col=2)
  abline(quantile(fit0$b,0.05),0,col=2,lty=2)
  abline(quantile(fit0$b,0.95),0,col=2,lty=2)

  b_mn[sim,1] <- mean(fit0$b)
  b_mn[sim,2] <- mean(fit1$b[,L1])
  b_mn[sim,3] <- mean(fit2$b[,L2])
  b_mn[sim,4] <- mean(beta.func[,n])
  
  b_sd[sim,1] <- sd(fit0$b)
  b_sd[sim,2] <- sd(fit1$b[,L1])
  b_sd[sim,3] <- sd(fit2$b[,L2])
  b_sd[sim,4] <- sd(beta.func[,n])
}

BIAS <- colMeans(b_mn-bX)
SD   <- colMeans(b_sd)
RMSE <- sqrt(colMeans((b_mn-bX)^2))
COV  <- colMeans(abs(b_mn-bX)<2*b_sd)
POW  <- colMeans(abs(b_mn)>2*b_sd)

BIAS_se <- apply(b_mn-bX,2,sd)/sqrt(nsims)
RMSE_se <- boot_rmse_se(b_mn-bX)
SD_se   <- apply(b_sd,2,sd)/sqrt(nsims)
COV_se  <- apply(abs(b_mn-bX)<2*b_sd,2,sd)/sqrt(nsims)
POW_se  <- apply(abs(b_mn)>2*b_sd,2,sd)/sqrt(nsims)

out     <- cbind(RMSE,RMSE_se,BIAS,BIAS_se,SD,SD_se,COV,COV_se)

save.image(paste0("DiscreteSim/sim",design,"_nsim",nsims,".RData"))

tex <- function(x,dig){
   x <- format(round(x,dig),nsmall=dig)

   cat(paste0(x[1],
          " (",
          x[2],
          ") & ",
          x[3],
          " (",
          x[4],
          ") & ",
          x[5],
          " (",
          x[6],
          ") & ",
          x[7],
          " (",
          x[8],
          ") \\"),"\n")
}

cat("design",design,"\n")
for(i in 1:4){
  tex(100*out[i,],1)
}

}
 
