###################################################################################################
# R file code_sec5_discretesim.R performs discrete-case simulation study
# Outputs of simulation study are saved under folder DiscreteSim
# 
# After simulation study is completed:
# Run this R script to make manuscript Table 1 & Fi2 & Fi3, and supp material Figure 1
# 
# Notes:
# - To save time, code_sec5_discretesim.R simulates two data sets for each scenario in Table 1.
# - The actual simulation study for Table 1 is performed on a high performance computer 
#   where each scenario has 500 simulated data sets and takes about 48 hours.
#   The actual simulation study outputs are saved in folder DiscreteSim
#   
#     For each simulated data set, we compare four methods:
#     the standard Leroux CAR model ("standard"), parsimonious bivariate CAR ("parametric"), 
#     semi-parametric with penalized complexity prior ("semip-PCP") and semi-parametric R2 ("semip-R2")
######################################################################################################################################################################################################

# Manuscript Table 1. Discrete-space simulation study results ####
tex <- function(names,x,dig){
  x <- format(round(x,dig),nsmall=dig)
  
  cat(paste0(names,x[1],
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
table_row <- c("standard  ","semip(PCP) ","semip(R2)  ","parametric ")

sink(paste0("Table_Figures/MainTable1.txt"))
print(c("RMSE","RMSE_se","BIAS","BIAS_se","SD","SD_se","COV","COV_se"))

for(design in 1:5){
  cat("design",design,"\n")
  # the R object 'out' contains cbind(RMSE,RMSE_se,BIAS,BIAS_se,SD,SD_se,COV,COV_se)
  load(paste0("DiscreteSim/sim",design,"_nsim500.RData"))
  
  for(i in 1:4){
    tex(table_row[i],100*out[i,],1)
  }
}
sink()

library(fields)
# Plots data on a square
square <- function(Y,main=""){
  library(viridis)
  library(fields)
  n <-sqrt(length(Y))
  image.plot(1:n,1:n,matrix(Y,n,n),col=viridis(100),
             main=main,xlab="",ylab="",cex.main=1.8,cex.axis=1.5,
             axis.args = list(cex.axis = 1.5))
}

# Supp material Figure 1. Simulated data in discrete space ####
pdf("Table_Figures/SuppMaterial_Fig1.pdf",height = 3,width = 9.5)
set.panel(1,3)
par(mar=c(3,4,3,5))
for(design in 2:5){
  n1    <- 40  # n1xn1 grid
  n     <- n1^2 # n locations
  lam   <- 0.95  
  bX    <- 0.5
  bZ    <- 0.5
  sigma <- 0.25
  sigx  <- sqrt(1.7)
  sigz  <- sqrt(1.0)
  
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
  
  sim = 1
  set.seed(sim*0820+design)
  
  w    <- 1/sqrt(1-lam+lam*D)
  X    <- sigx*G%*%rnorm(n,0,w) # simulate X from Leroux CAR model
  Z    <- bXZ*W%*%X + sigz*G%*%rnorm(n,0,w)
  Y    <- bX*X + bZ*Z + rnorm(n,0,sigma)
  
  if(sim==1){ ## plot the simulated data
     square(X,main="X")
     square(Y,main="Y")
     square(Z,main="Z")
  }
}
dev.off()



# Figure 3 shows beta(omega) for the discrete simulation study ####
pdf("Table_Figures/Main_Fig3.pdf")
par(cex.lab=1.5,cex.main=2,cex.axis=1.5,mar=c(5, 4+0.5, 4, 2) + 0.1)
ylim=c(0.25,1.3)
for(design in 1:5){
  load(paste0("DiscreteSim/sim",design,"_nsim500.RData"))
  
  if(design==1){title <- expression("(a)" ~ beta[XZ]==0)}
  if(design==2){title <- expression("(b)" ~ phi==1 ~ "and" ~ beta[XZ]==1)}
  if(design==3){title <- expression("(c)" ~ phi==1 ~ "and" ~ beta[XZ]==2)}
  if(design==4){title <- expression("(d)" ~ phi==2 ~ "and" ~ beta[XZ]==1)}
  if(design==5){title <- expression("(e)" ~ phi==2 ~ "and" ~ beta[XZ]==2)}
  
  #beta.hat1 contains beta(omega) estimates from semi-par (PCP)
  Q1 <- apply(t(beta.hat1),2,quantile,c(0.025,0.500,0.975))
  #beta.hat3 contains beta(omega) estimates from spectral-par
  Q3 <- apply(t(beta.hat3),2,quantile,c(0.025,0.500,0.975))
  #b_mn[,1] contains beta estimates from standard Leroux CAR
  q <- quantile(b_mn[,1],c(0.025,0.500,0.975))
  
  plot(D,Q1[2,],lwd=2,col=3,ylim=ylim,type="l",
       xlab=expression("Eigenvalue," ~ omega[k]),ylab=expression(beta[k]),
       main=title)
  lines(D,Q1[1,],lty=2,col=3)
  lines(D,Q1[3,],lty=2,col=3)
  # use rev(D) for plotting, since beta(omega) in the parametric model used increasing order of eigen values (sec4.1)
  lines(rev(D),Q3[1,],lty=2,col=4)
  lines(rev(D),Q3[2,],lwd=2,col=4)
  lines(rev(D),Q3[3,],lty=2,col=4)
  abline(q[2],0,lwd=2,col=2)
  abline(q[1],0,lty=2,col=2)
  abline(q[3],0,lty=2,col=2)
  abline(0.5,0,lwd=2,col=1)
  
  if(design==1){
    legend("topright",c("True","Standard","Semiparametric","Parametric"),col=1:4,lwd=2,cex=1.5,bty="n")
  }
}
dev.off()




# Figure 2 Correlations in the spectral domain for the simulation study ####
library(fields)
nsims <- 500 # use 500 replicates to produce Table 1 result
n1    <- 40  # n1xn1 grid
n     <- n1^2 # n locations
lam   <- 0.95  
bX    <- 0.5
bZ    <- 0.5
sigma <- 0.25
sigx  <- sqrt(1.7)
sigz  <- sqrt(1.0)

# Define rook neighborhood structure 
s     <- expand.grid(1:n1,1:n1)
A     <- rdist(s)==1 # adjacency
R     <- diag(rowSums(A))-A # precision of ICAR
E     <- eigen(R) # eigen component
G     <- E$vec
D     <- E$val
rm(E,R)

corest <- c()
corsp <- c()
for(design in 2:5){
  if(design==1){bXZ <- 0; bw=1}
  if(design==2){bXZ <- 1; bw=1}
  if(design==3){bXZ <- 2; bw=1}
  if(design==4){bXZ <- 1; bw=2}
  if(design==5){bXZ <- 2; bw=2}
  
  # Kernel smoothing matrix
  W     <- exp(-(rdist(s)/bw)^2)
  W     <- sweep(W,1,rowSums(W),"/")
  
  outXstar <- c()
  outZstar <- c()
  outcor <- c()
  for(sim in 1:nsims){
    set.seed(sim*0820+design)
    
    w    <- 1/sqrt(1-lam+lam*D)
    X    <- sigx*G%*%rnorm(n,0,w) # simulate X from Leroux CAR model
    Z    <- bXZ*W%*%X + sigz*G%*%rnorm(n,0,w)
    Y    <- bX*X + bZ*Z + rnorm(n,0,sigma)
    
    # assume same adjacency structure and project all variables 
    # into the spectral domain using graph Fourier transform
    Xstar <- t(G)%*%X
    Zstar <- t(G)%*%Z
    outXstar <- cbind(outXstar,Xstar)
    outZstar <- cbind(outZstar,Zstar)
    outcor  <- c(outcor,cor(X,Z))
  }
  corest_design <- unlist(lapply(1:n,function(i) cor(outXstar[i,],outZstar[i,])))
  corest <- cbind(corest,corest_design)
  corsp  <- cbind(corsp,outcor)
}

pdf("Table_Figures/Main_Fig2.pdf",width = 5, height=4) 
par(mar=c(4,4,1,1))
matplot(D, lty = 1,corest,type="l",xlab=expression("Eigenvalue," ~ omega[k]),ylab="Correlation",ylim = c(-0.2,1))
legend("topright",c(expression(phi==1~ " " ~ beta[XZ]==1~ " "),
         expression(phi==1 ~ " "~ beta[XZ]==2~ "   "),
         expression(phi==2 ~ " "~ beta[XZ]==1~ "   "),
         expression(phi==2 ~ " "~ beta[XZ]==2~ "   ")),col=1:4,lty=1,bty = "n")
dev.off()
print(round(colMeans(corsp),2))
