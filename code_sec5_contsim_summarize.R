###################################################################################################
# R file code_sec5_contsim.R performs continuous-case simulation study
# Outputs of simulation study are saved under folder ContSim
# 
# After simulation study is completed:
# Run this R script to make manuscript Figure 1 and supp material Table 1
# 
# Notes:
# - To save time, code_sec5_contsim.R simulates 1 data set for one scenario in supp material Table 1.
# - The actual simulation study for supp material Table 1 is performed on a high performance computer 
#   where each scenario has 100 simulated data sets. Simulation study outputs are saved in ContSim
#   
#     For each simulated data set, we compare four methods:
#     the standard Matern, flexible Matern, parsimonious Matern and semiparametric
# 
###################################################################################################

# supp material Table 1. Continuous-space simulation study results ####
out   <- NULL
nsims <- 100
bX    <- 1
for(case in c(9,5:8)){
  beta_sdsave = beta_msave = beta_lcisave = beta_ucisave = matrix(NA,nr=4,nc = 100)
  for(sim in c(1:nsims)){
    # store results
    file  = paste0("ContSim/RMsen_par_parsi_bs_case",case,"_design2_id",sim,".RData")
    if(file.exists(file)){
      load(file=file)
      del = seq(1,40,by=5)[which.min(alldic)]
      fitsp_bs=get(paste0("fitsp_bs_",del))
      
      iters = nrow(naive$betaX)
      allbetax <- cbind(naive$betaX,BMPar.Bayes$betaX,BMPar.parsi.Bayes$betaX,fitsp_bs$betaX[,1])
      allbetax <-  allbetax[(0.25*iters+1):iters,]
      beta_sdsave[,sim] <- apply(allbetax,2,sd)
      beta_msave[,sim]  <- colMeans(allbetax)
      beta_lcisave[,sim]<- apply(allbetax,2, quantile,p = 0.025)
      beta_ucisave[,sim]<- apply(allbetax,2, quantile,p = 0.975)
      rm(alldic,naive,BMPar.Bayes,BMPar.parsi.Bayes,fitsp_bs)
    }
  }
  BIAS <- rowMeans(beta_msave-bX,na.rm=T)
  SD   <- rowMeans(beta_sdsave,na.rm=T)
  RMSE <- sqrt(rowMeans((beta_msave-bX)^2,na.rm=T))
  COV  <- rowMeans(beta_lcisave<bX & bX<beta_ucisave,na.rm=T)
  
  boot_rmse_se <- function(x,iters=100){
    y <- matrix(0,iters,ncol(x))
    for(i in 1:iters){
      samp  <- sample(1:nrow(x),replace=TRUE)
      y[i,] <- sqrt(colMeans(x[samp,]^2,na.rm=T))
    }
    return(apply(y,2,sd))}
  
  BIAS_se <- apply(beta_msave-bX,1,sd,na.rm=T)/sqrt(nsims)
  RMSE_se <- boot_rmse_se(t(beta_msave-bX))
  SD_se   <- apply(beta_sdsave,1,sd,na.rm=T)/sqrt(nsims)
  COV_se  <- apply(beta_lcisave<bX & bX<beta_ucisave,1,sd,na.rm=T)/sqrt(nsims)
  
  out_case     <- cbind(RMSE,RMSE_se,BIAS,BIAS_se,SD,SD_se,COV,COV_se)
  out <-rbind(out,out_case)
}
# summarize results
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

table_row <- rep(c("stand ","flexi ","parsi ","semip "),5)
sink(paste0("Table_Figures/SuppMaterialTable1.txt"))
print(c("RMSE","RMSE_se","BIAS","BIAS_se","SD","SD_se","COV","COV_se"))
for(i in 1:nrow(out)){
  tex(table_row[i],100*out[i,],1)
}
sink()


# Figure 1 Example confounder adjustment for the bivariate Matern ####
library(fields)
source("MCMC_GP.R")

set.seed(0)
d         <- as.matrix(expand.grid(1:50,1:50))
distmat   <- rdist(d)
sigmax <- 1; rangex <- 1; nux <- 1
sigx      <- maternFun(distmat=distmat,covparms=c(sigmax,rangex,nux,0)) #covparms: var,range,nu,nug/var
X         <- t(chol(sigx))%*%rnorm(nrow(sigx))
invsig.x  <- solve(sigx,X)
sigx.inv  <- solve(sigx)

pdf("Table_Figures/Main_Fig1.pdf",width = 3,height= 3)
par(mar=c(1,1,2,1))
for(C in c(1,3,5)){
  sigxz     <- maternFun(distmat, c(sigmax,rangex,C*nux,0))
  zhat      <- sigxz%*%invsig.x # covparms: var, range,nu,nug/var
  image(1:50,1:50,matrix(scale(zhat),50),xaxt='n',yaxt='n',
        main= substitute(paste(nu[XZ]," = ",nn,nu[X]), list(nn=C)),
        xlab="",ylab="",col=viridis(64),)
}
dev.off()