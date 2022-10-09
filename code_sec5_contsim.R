###################################################################################################
# example code to perform continuous-space simulation study
# X ~ GP_Matern(0,sigx^2,smoothness=nu1,range=s1), we use Leroux CAR here 
# Z|X ~ GP_Matern(bXZ*W*X,smoothness=nu1,range=s1)
# Y|X,Z ~ N(bX*X + bZ*Z, sigma^2)
# W is kernel smoothing matrix with bandwidth phi=bw
#
# Notes:
# - To save time, this script simulates 1 data set for one scenario in supplementary material Table 1.
#   For each simulated data set, we compare four methods:
#     the standard Matern, flexible Matern, parsimonious Matern and semiparametric
#   This takes about 15min to complete
# 
# - The actual simulation study for supplementary material Table 1 is performed on a high performance computer 
#   where each scenario has 100 simulated data sets. 
#   each data set takes about 5 hours 
###################################################################################################
library(fields)
library(RandomFields)
library(truncnorm)
library(mvnfast)
library(splines)

source("MCMC_GP.R")  #utility functions for the continuous case
source("SP_vecc.R")  #standard Matern spatial model (naive)
source("BMcausal_vecc.R") #common range Flexible Bivariate Matern and common range Parsimonious Bivariate Matern
source("SemiPcausal_vecc.R") #semiparametric model 

case = 5 # cases 5-9 are presented in the table
nsim = 1 # nsim = 100; number of replicates. 
iters <- 2e2 # small mcmc iteration to make simulation run faster. Longer MCMC (iters=2e3) were ran for manuscript results

for(sim in 1:nsim){
n1    <- 100 # grid size; simulate on a large grid, then down sampled to 23x23 
n     <- n1^2
bX    <- 1
bZ    <- 1
sigma <- 0.25
design = 2

# simulate data on a grid
if(case==5){bXZ <- 1; bw=1/15;convol  = T}
if(case==6){bXZ <- 2; bw=1/15;convol  = T}
if(case==7){bXZ <- 1; bw=2/15;convol  = T}
if(case==8){bXZ <- 2; bw=2/15;convol  = T}
if(case==9){bXZ <- 0; bw=1/15;convol  = T}
s   <- as.matrix(expand.grid(seq(0,1,l=n1),seq(0,1,l=n1)))

# truealpha <- function(omega,bw,bXZ) bXZ*1/(2*pi)*exp(-omega^2/(4/bw^2)) 
# omega = seq(0,20*pi,l=100)
# plot(omega,truealpha(omega,bw,bXZ),type="l"); abline(v=pi/(1/15))

RFoptions(install="no")
RFoptions(seed=sim*0820+50+case)
set.seed(seed=sim*0820+50+case)

# simulate from exponential with convolution
nu1=0.5; s1=0.1
model <- RMbiwm(nu=c(nu1,10,nu1), c=c(1,0,1), s=c(s1,s1,s1))
foo   <- RFsimulate(model,x=seq(0,1,l=n1),y=seq(0,1,l=n1))
X     <- scale(foo$variable1)
delta <- foo$variable2
W     <- exp(-(rdist(s)/bw)^2)
W     <- 1/(pi*bw^2)*W*(1/n1^2)
EZX   <- bXZ*W%*%X 
Z     <- EZX + delta
Y     <- bX*X + bZ*Z + rnorm(n,0,sigma)


# sample subgrided data
if(design ==2) rnd.idx <- matrix(1:n,n1)[seq(6,(n1-5),by=4),seq(6,(n1-5),by=4)]
# smaller
if(design ==1) rnd.idx <- matrix(1:n,n1)[seq(6,(n1-5),by=6),seq(6,(n1-5),by=6)]

Z.s <- Z[rnd.idx];X.s <- X[rnd.idx];Y.s <- Y[rnd.idx];s.s <-s[rnd.idx,]

beta_sd <-beta_mn <- matrix(0,100,30)
beta_mn[sim,1] <- coef(lm(Y.s~X.s))[2]
beta_mn[sim,2] <- coef(lm(Y.s~X.s+Z.s))[2]
beta_sd[sim,1] <- summary(lm(Y.s~X.s))$coef[2,2]
beta_sd[sim,2] <- summary(lm(Y.s~X.s+Z.s))$coef[2,2]

### standard Matern spatial model (naive) ####
naive = SP_vecc(matrix(Y.s,nc=1),matrix(X.s,nc=1),d=s.s,L=1,iters=iters,nn=20,nures=1)
beta_mn[sim,3] <- mean(naive$betaX[(0.5*iters):iters])
beta_sd[sim,3] <- sd(naive$betaX[(0.5*iters):iters])


### common range Flexible Bivariate Matern ####
BMPar.Bayes = BMcausal_crho_vecc(Y.s,X.s,d=s.s,nux=nu1,rhox=s1,iters=iters,burn=0.5*iters,update =  iters*0.2)
print(BMPar.Bayes$time)
beta_mn[sim,4] <- mean(BMPar.Bayes$betaX[(0.5*iters):iters])
beta_sd[sim,4] <- sd(BMPar.Bayes$betaX[(0.5*iters):iters])

### common range Parsimonious Bivariate Matern ####  
BMPar.parsi.Bayes = BMcausal_vecc(Y.s,X.s,d=s.s,nux=nu1,rhox=s1,iters=iters,burn=0.5*iters,update =  iters*0.2)
print(BMPar.parsi.Bayes$time)
beta_mn[sim,5] <- mean(BMPar.parsi.Bayes$betaX[(0.5*iters):iters])
beta_sd[sim,5] <- sd(BMPar.parsi.Bayes$betaX[(0.5*iters):iters])

savebetaX = NULL
alldic = NULL

### semiparametric model ####
s.gri.idx  <- matrix(1:n,n1)[seq(1,(n1),by= 2),seq(1,(n1),by= 2)]
s.expand = s[s.gri.idx,]
X.expand = X[s.gri.idx]
setup = setup_bspline_grid(s.s,s.expand,X.expand,dx = 1)

for(del in seq(1,40,by=5)){
  zhat = getzhat_bspline_grid(s.s, s.expand,L,X=X.expand,del=del,setup = setup) #calculate confounder adjustment zhat
  fitsp_bs = SemiPcausal_vecc(Y.s,X.s,d=s.s,L=ncol(zhat$zhat),
                              zhat=zhat$zhat,iters=iters,burn =0.5*iters,nn=20)
  alldic = c(alldic, fitsp_bs$DIC$DIC)
  savebetaX = cbind(savebetaX, fitsp_bs$betaX[,1])
  
  assign(paste0("fitsp_bs_",del),fitsp_bs)
  assign(paste0("zhat_",del),zhat)
}

file  = paste0("ContSim/testrun_RMsen_par_parsi_bs_case",case,"_design",design,"_id",sim,".RData")
save(naive,BMPar.Bayes,BMPar.parsi.Bayes,X,Y,Z,s,rnd.idx,alldic,
     fitsp_bs_1,fitsp_bs_6,fitsp_bs_11,fitsp_bs_16,fitsp_bs_21,fitsp_bs_26,
     fitsp_bs_31,fitsp_bs_36,
     zhat_1,zhat_6,zhat_11,zhat_16,zhat_21,zhat_26,zhat_31,zhat_36,
     file=file)

}

# plot MCMC samples from beta_x posterior distribution
set.panel(1,1)
del = seq(1,40,by=5)[which.min(alldic)]
fitsp_bs=get(paste0("fitsp_bs_",del))
allbetax <- cbind(naive$betaX,BMPar.Bayes$betaX,BMPar.parsi.Bayes$betaX,fitsp_bs$betaX[,1])
allbetax <-  allbetax[(0.25*iters+1):iters,]
colnames(allbetax) <- c("standard","flexible","parsimonious","semi-par")
boxplot(allbetax,main="beta_x posterior")
abline(h=bX)

