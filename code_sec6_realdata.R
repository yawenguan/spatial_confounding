#############################################################################################
# this script reproduces the content in sec 6 in the paper: 
# fig 4 and 5, section 6.1 (Scottish lip cancer); fig 6 and 7 section 6.2 (covid-19 USA)
# ## note: results in the paper are run longer; draws=10000, burn=0.5*draws
#############################################################################################
draws = 1e3
# required libraries to fit model
library(INLA)
library(eCAR) 
# required libraries for plots
library(usmap)
library(ggplot2)
library(viridis)
library(ggthemes)
library(broom) # "tidy" function
library(dplyr) # "left_join" function
library(maptools)

## Lip cancer Scotland data analysis####
# sec 6.1 in the paper
# Four methods are compared:
# Spectral Parametric
# Spectral Semiparametric
# Standard non-spatial
# Standard spatial
##

data(lipcancer)
W <- lipcancer$neighborhood.Matrix
y <- lipcancer$data$observed
x <- lipcancer$data$pcaff
E <- lipcancer$data$expected

#### Spectral Semiparametric ####
# fit the semipar with L=10,20,30,40,50 and select via DIC
fit.semi.10 <- semipar.eCAR.Leroux(y=y, x=x, W=W, E=E, C=NULL,
                                   pcprior.sd = c(0.1,1),
                                   model="Poisson",
                                   L=10,
                                   method = "spectral",
                                   verbose=T)

fit.semi.20 <- semipar.eCAR.Leroux(y=y, x=x, W=W, E=E, C=NULL,
                                   pcprior.sd = c(0.1,1),
                                   model="Poisson",
                                   L=20,
                                   method = "spectral",
                                   verbose=T)

fit.semi.30 <- semipar.eCAR.Leroux(y=y, x=x, W=W, E=E, C=NULL,
                                   pcprior.sd = c(0.1,1),
                                   model="Poisson",
                                   L=30,
                                   method = "spectral",
                                   verbose=T)

fit.semi.40 <- semipar.eCAR.Leroux(y=y, x=x, W=W, E=E, C=NULL,
                                   pcprior.sd = c(0.1,1),
                                   model="Poisson",
                                   L=40,
                                   method = "spectral",
                                   verbose=T)

fit.semi.50 <- semipar.eCAR.Leroux(y=y, x=x, W=W, E=E, C=NULL,
                                   pcprior.sd = c(0.1,1),
                                   model="Poisson",
                                   L=50,
                                   method = "spectral",
                                   verbose=T)

fit.semi.10$DIC   # best
fit.semi.20$DIC
fit.semi.30$DIC
fit.semi.40$DIC
fit.semi.50$DIC

#### Spectral Parametric ####
set.seed(101)
fit.par <- par.eCAR.Leroux(y=y, x=x, W=W, E=E, C=NULL, model="Poisson",
                           draws=10000, burn=5000, thin=1, verbose=FALSE,
                           joint_prior_lamx_lamz=FALSE, m=0, s2=4)

#### Standard non-spatial ####
# Poisson regression
r.standard.pois <- inla(y ~   x,
                        data=list(y=y, x=x),
                        family="poisson", E=E,
                        control.predictor = list( compute=TRUE),
                        control.compute = list(config=TRUE, dic=TRUE),
                        verbose=F)

# calculate posterior mean and 95% intervals
mean.beta.standard1 <- inla.emarginal(function(x) exp(x),
                                      r.standard.pois$marginals.fixed$x)
q975.beta.standard1 <- inla.qmarginal(0.975,
                                      inla.tmarginal(function(x) exp(x),
                                                     r.standard.pois$marginals.fixed$x))
q025.beta.standard1 <- inla.qmarginal(0.025,
                                      inla.tmarginal(function(x) exp(x),
                                                     r.standard.pois$marginals.fixed$x))

#### Standard spatial ####
# Poisson regression with residuals modeled as CAR Leroux
M <- rowSums(W)
R <- diag(M)-W
r.standard.leroux <- inla(y ~   x +
                            f(id.z, model = "generic1",
                              Cmatrix = Diagonal(x=1, n=56)-inla.as.sparse(R),
                              hyper = list(prec=list(
                                prior="pc.prec",
                                param=c(1/0.31, 0.01)))),
                          data=list(y=y, x=x, id.z=1:56),
                          family="poisson", E=E,
                          control.predictor = list( compute=TRUE),
                          control.compute = list(config=TRUE, dic=TRUE),
                          verbose=F)

# calculate posterior mean and 95% intervals
mean.beta.standard2 <- inla.emarginal(function(x) exp(x),
                                      r.standard.leroux$marginals.fixed$x)
q975.beta.standard2 <- inla.qmarginal(0.975,
                                      inla.tmarginal(function(x) exp(x),
                                                     r.standard.leroux$marginals.fixed$x))
q025.beta.standard2 <- inla.qmarginal(0.025,
                                      inla.tmarginal(function(x) exp(x),
                                                     r.standard.leroux$marginals.fixed$x))

#### Compare the four methods ####
# Fig. 5: Effect of percent of workforce in AFF on lip cancer in Scotland: ####
pdf("Table_Figures/Main_Fig5.pdf")
# fig 5 in the paper
par(mar=c(5.5, 5.5, 4, 2) + 0.1)
# Spectral Semiparametric fit; L=10 (best DIC)
plot(fit.semi.10$beta_omega[,"omega"], fit.semi.10$beta_omega[,"beta.mn"], type="l",
     ylab=expression(exp(beta[k])),
     xlab=expression(paste("", omega[k])),
     ylim=c(0.85,1.3),
     lwd=2, cex.axis=1.5, cex.lab=2, col=3)
lines(fit.semi.10$beta_omega[,"omega"], fit.semi.10$beta_omega[,"beta.q025"], lty=2, col=3)
lines(fit.semi.10$beta_omega[,"omega"], fit.semi.10$beta_omega[,"beta.q975"], lty=2, col=3)

# Spectral Parametric fit using par_eCAR
lines(fit.par$beta_omega[,4], fit.par$beta_omega[,1], lty=1, col=1, lwd=2)
lines(fit.par$beta_omega[,4], fit.par$beta_omega[,2], lty=2, col=1)
lines(fit.par$beta_omega[,4], fit.par$beta_omega[,3], lty=2, col=1)

# Standard non-spatial fit
abline(h=mean.beta.standard1,  col=2, lwd=2)
abline(h=q975.beta.standard1, lty=2, col=2)
abline(h=q025.beta.standard1, lty=2, col=2)

# Standard Spatial fit; use CAR Leroux for random effects
abline(h=mean.beta.standard2,  col=4, lwd=2)
abline(h=q975.beta.standard2, lty=2, col=4)
abline(h=q025.beta.standard2, lty=2, col=4)

# exp(beta)=1; no effect
abline(h=1, lty=3, lwd=3, col="orange")

legend("topright",
       legend=c("Spectral Parametric",
                "Spectral Semiparametric",
                "Standard non-spatial",
                "Standard spatial"),
       col=c(1,3,2,4),
       lty=1, lwd=c(2,2),
       bty="n", cex = 1.7)
dev.off()

# Fig. 4: in the paper; Maps for Scotland lip cancer: ####
# standard mortality ratio (left) and workforce in AFF (right)
pdf("Table_Figures/Main_Fig4.pdf")
library(ggplot2)
library(broom) # "tidy" function
library(dplyr) # "left_join" function
library(viridis)
library(maptools)

load(file="lip.Rdata")
source("combine.data.shapefile.R") # "combine.data.shapefile" function
lipdata$name <- (lipdbf$dbf[,1])
lipdbf$dbf <- lipdbf$dbf[ ,c(2,1)]
data.combined <- combine.data.shapefile(data=lipdata, shp=lipshp, dbf=lipdbf)

sp.df <- tidy(data.combined, region = "name")
data.combined$id <- data.combined@data$name
sp.df <- left_join(sp.df,data.combined@data,by = "id")
sp.df$RR <- sp.df$observed/sp.df$expected

# fig 4 left panel
ggplot(data = sp.df, aes(x = long, y = lat, group = group, fill=RR)) +
  geom_polygon(color="black") +
  theme_bw() +
  labs(fill = "Standardized Mortality Ratio") +
  scale_fill_viridis(discrete = FALSE) +
  theme(legend.position = "bottom",
        legend.title=element_text(size=24),  # font size of the legend
        legend.text=element_text(size=14),
        axis.title.x=element_blank(),  # remove axis, title, ticks
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),  # remove axis, title, ticks
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.width= unit(2,"line"),
        legend.key.height=unit(2,"line"))

# fig 4 right panel 
ggplot(data = sp.df, aes(x = long, y = lat, group = group, fill=pcaff)) +
  geom_polygon(color="black") +
  theme_bw() +
  labs(fill = "% Workforce in AFF") +
  scale_fill_viridis(discrete = FALSE) +
  theme(legend.position = "bottom",
        legend.title=element_text(size=24),  # font size of the legend
        legend.text=element_text(size=14),
        axis.title.x=element_blank(),  # remove axis, title, ticks
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),  # remove axis, title, ticks
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.width= unit(2,"line"),
        legend.key.height=unit(2,"line"))
dev.off()




##COVID US data analysis#######
# sec 6.2 in the paper
##
# load data 'county_neighborhood_matrix.txt' and 'covid_data2.txt'
adj.matrix    <- as.matrix(read.table(file="county_neighborhood_matrix.txt"))
covid <- read.table(file="covid_data2.txt", header = TRUE)

# rescale covariates
covid$education <- scale(covid$education)
covid$poverty <- scale(covid$poverty)
covid$pct_owner_occ  <-  scale(covid$pct_owner_occ)
covid$pct_blk <-  scale(covid$pct_blk)
covid$hispanic <-  scale(covid$hispanic)
covid$prime_pecent <- scale(covid$prime_pecent)
covid$mid_pecent <- scale(covid$mid_pecent)
covid$older_pecent <- scale(covid$older_pecent)
covid$obese <- scale(covid$obese)
covid$smoke <- scale(covid$smoke)
covid$popdensity <- scale(covid$popdensity)
covid$medhouseholdincome <- scale(covid$medhouseholdincome)
covid$medianhousevalue <- scale(covid$medianhousevalue)
covid$mean_summer_temp <- scale(covid$mean_summer_temp)
covid$mean_winter_temp <- scale(covid$mean_winter_temp)
covid$mean_summer_rm <- scale(covid$mean_summer_rm)
covid$mean_winter_rm <- scale(covid$mean_winter_rm)
covid$beds <- scale(covid$beds)
covid$date_since_social <- scale(covid$date_since_social)
covid$date_since <- scale(covid$date_since)
pop <- covid$population
pop.scaled <- scale(covid$population, center = FALSE)

X.cov <- cbind(covid$popdensity, covid$poverty, covid$medianhousevalue,
               covid$medhouseholdincome,covid$pct_owner_occ,covid$education,
               covid$pct_blk,covid$hispanic,covid$older_pecent,covid$prime_pecent,covid$mid_pecent,
               covid$date_since_social,
               covid$date_since,covid$beds,covid$obese,covid$smoke,
               covid$mean_summer_temp,covid$mean_winter_temp,covid$mean_summer_rm,covid$mean_winter_rm)




## Analysis on the reduced dataset ####
# only counties that reported at least 10 cases are considered (1977 counties)
y.na <- covid$Deaths
ind <- covid$Confirmed>=10
y.na[!ind] <- NA
pop.scaled.na <- pop.scaled
pop.scaled.na[!ind] <- NA

### Spectral Semiparametric ####
### The commented code below try different L values to choose an 'optimal' L; 
### To save time, YOU CAN SKIP THIS, as we have chosen L=8 according to DIC 
# L0 <- 6 #try a few different L values
# step.L <- 2
# max.numL <- 10
# res <- vector("list", length = max.numL)
# 
# for(j in 1:max.numL){
#   if (j==1) {
#     L <- L0
#     stored.L <- L
#   } else {
#     stored.L <- L
#     L <- L + step.L
#   }
#   
#   res[[j]] <- semipar.eCAR.Leroux(y=y.na,
#                                   x=covid$mean_pm25,
#                                   W=adj.matrix,
#                                   E=pop.scaled,
#                                   C = X.cov,
#                                   names.covariates = c(
#                                     "p_dens","poverty","house_val","house_inc",
#                                     "owner","education","blk","hisp",
#                                     "old", "prime","mid",
#                                     "d_since_social", "d_since",
#                                     "bed","obesity","smoke","s_temp",
#                                     "w_temp","s_rm","w_rm"),
#                                   pcprior.sd = c(0.1,1), s2=10,
#                                   model="Negative Binomial",
#                                   L=L,
#                                   method= "spectral",
#                                   verbose=F)
#   
#   
#   print(paste("iter: ", j," ; L=", L," : DIC=", round(res[[j]]$DIC,3), sep=""))
#   if (j == 1) {
#     stored.DIC <- res[[j]]$DIC
#   } else {
#     if (stored.DIC < res[[j]]$DIC) {
#       print(paste("optimal L=", stored.L," : DIC=",round(stored.DIC,3),sep=""))
#       opt.L <- stored.L
#       break
#     }
#     stored.DIC <- res[[j]]$DIC
#   }
# }

### fit Spectral Semiparametric with L=8 ####
L <- 8 # chosen according to DIC
nb.spectral.semipar.reduced <- semipar.eCAR.Leroux(y=y.na,
                                                   x=covid$mean_pm25,
                                                   W=adj.matrix,
                                                   E=pop.scaled,
                                                   C = X.cov,
                                                   names.covariates = c(
                                                     "p_dens","poverty","house_val","house_inc",
                                                     "owner","education","blk","hisp",
                                                     "old", "prime","mid",
                                                     "d_since_social", "d_since",
                                                     "bed","obesity","smoke","s_temp",
                                                     "w_temp","s_rm","w_rm"),
                                                   pcprior.sd = c(0.1,1), s2=10,
                                                   model="Negative Binomial", L=L,
                                                   method= "spectral",
                                                   verbose=T)

### Standard Spatial or 'naive' ####
nb.naive.reduced <- semipar.eCAR.Leroux(y=y.na,
                                        x=covid$mean_pm25,
                                        W=adj.matrix,
                                        E=pop.scaled,
                                        C = X.cov,
                                        names.covariates = c(
                                          "p_dens","poverty","house_val","house_inc",
                                          "owner","education","blk","hisp",
                                          "old", "prime","mid",
                                          "d_since_social", "d_since",
                                          "bed","obesity","smoke","s_temp",
                                          "w_temp","s_rm","w_rm"),
                                        pcprior.sd = c(0.1,1), s2=10,
                                        model="Negative Binomial",
                                        method= "naive",
                                        verbose=T)

### Spectral Parametric ####
set.seed(101)
nb.spectral.par.reduced <- par.eCAR.Leroux(y=covid$Deaths[ind], 
                                           x=covid$mean_pm25[ind], 
                                           W=adj.matrix[ind,ind], 
                                           E=pop.scaled.na[ind], 
                                           C=X.cov[ind,], 
                                           model="Negative Binomial",
                                           draws=draws, burn=0.5*draws, thin=1, verbose=TRUE,
                                           joint_prior_lamx_lamz=FALSE, m=0, s2=4)
## note: results in the paper are run longer; draws=10000, burn=5000

# compute the fit for Parametric 
# note: doing it manually using the output of par.eCAR.Leroux() cause in this case we have NA in the data vector; 
# par.eCAR.Leroux() not yet optimized for situations where we have missing values
emp.hpd <- function(x, conf = 0.95) {
  conf <- min(conf, 1 - conf)
  n <- length(x)
  nn <- round(n * conf)
  x <- sort(x)
  xx <- x[(n - nn + 1):n] - x[1:nn]
  m <- min(xx)
  nnn <- which(xx == m)[1]
  return(c(x[nnn], x[n - nn + nnn]))
}
nout <- 500
Dseq <- seq(0, 15, length = 1000)
c.beta <- matrix(NA, nrow = nout, ncol = length(Dseq))
out <- nb.spectral.par.reduced$posterior_draws
for (i in 1:nout) {
  c.beta[i, ] <- out$beta[i] + 
    out$gamma[i] * sqrt((1 - out$lamx[i] + out$lamx[i] * Dseq)/
                          (1 - out$lamz[i] + out$lamz[i] * Dseq))
}
beta.mn.spectral.par.reduced <- matrix(apply(exp(c.beta), 2, function(x) mean(x)), 
                                       nrow = 1000, byrow = TRUE)
beta.q025.spectral.par.reduced <- matrix(apply(exp(c.beta), 2, function(x) emp.hpd(x))[1, ], nrow = 1000, byrow = TRUE)
beta.q975.spectral.par.reduced <- matrix(apply(exp(c.beta), 2, function(x) emp.hpd(x))[2, ], nrow = 1000, byrow = TRUE)


##  Fig. 7 panel (b)  Analysis on the reduced dataset ####
pdf("Table_Figures/Main_Fig7b.pdf")
plot(nb.spectral.semipar.reduced$beta_omega[,"omega"],
     nb.spectral.semipar.reduced$beta_omega[,"beta.mn"], type="l", main = "(b)",
     ylab=expression(exp(beta[k])),
     xlab=expression(paste("", omega[k])), 
     ylim=c(0.85,1.7),  lwd=2, cex.axis=1.5, cex.lab=2.5, cex.main=2, col=3)
lines(nb.spectral.semipar.reduced$beta_omega[,"omega"],
      nb.spectral.semipar.reduced$beta_omega[,"beta.q025"], col=3, lty=2, lwd=1)
lines(nb.spectral.semipar.reduced$beta_omega[,"omega"],
      nb.spectral.semipar.reduced$beta_omega[,"beta.q975"], col=3, lty=2, lwd=1)

lines(Dseq,
      beta.mn.spectral.par.reduced, col=1, lwd=2)
lines(Dseq,
      beta.q025.spectral.par.reduced, col=1, lwd=1, lty=2)
lines(Dseq,
      beta.q975.spectral.par.reduced, col=1, lwd=1, lty=2)

lines(nb.naive.reduced$beta_omega[,"omega"],
      nb.naive.reduced$beta_omega[,"beta.mn"], col=2, lwd=2)
lines(nb.naive.reduced$beta_omega[,"omega"],
      nb.naive.reduced$beta_omega[,"beta.q025"], col=2, lty=2, lwd=1)
lines(nb.naive.reduced$beta_omega[,"omega"],
      nb.naive.reduced$beta_omega[,"beta.q975"], col=2, lty=2, lwd=1)
# exp(beta)=1; no effect
abline(h=1, lty=3, lwd=3, col="orange")

legend("topright",
       legend=c("Spectral Semiparametric",
                "Spectral Parametric",
                "Standard Spatial"),
       col=c(3,1,2),
       lty=1, lwd=c(2,2,2),
       bty="n", cex = 1.7)
dev.off()

## Analysis on the entire dataset ####
# all counties are considered (3109 counties)

### Spectral Semiparametric ####
### The commented code below try different L values to choose an 'optimal' L; 
### To save time, YOU CAN SKIP THIS, as we have chosen L = 4 according to DIC 
# L0 <- 4 # try a few different L values
# step.L <- 2
# max.numL <- 10
# res <- vector("list", length = max.numL)
# 
# for(j in 1:max.numL){
#   if (j==1) {
#     L <- L0
#     stored.L <- L
#   } else {
#     stored.L <- L
#     L <- L + step.L
#   }
#   
#   res[[j]] <- semipar.eCAR.Leroux(y=covid$Deaths,
#                                   x=covid$mean_pm25,
#                                   W=adj.matrix,
#                                   E=pop.scaled,
#                                   C = X.cov,
#                                   names.covariates = c(
#                                     "p_dens","poverty","house_val","house_inc",
#                                     "owner","education","blk","hisp",
#                                     "old", "prime","mid",
#                                     "d_since_social", "d_since",
#                                     "bed","obesity","smoke","s_temp",
#                                     "w_temp","s_rm","w_rm"),
#                                   pcprior.sd = c(0.1,1), s2=10,
#                                   model="Negative Binomial",
#                                   L=L,
#                                   method= "spectral",
#                                   verbose=F)
#   
#   
#   print(paste("iter: ", j," ; L=", L," : DIC=", round(res[[j]]$DIC,3), sep=""))
#   if (j == 1) {
#     stored.DIC <- res[[j]]$DIC
#   } else {
#     if (stored.DIC < res[[j]]$DIC) {
#       print(paste("optimal L=", stored.L," : DIC=",round(stored.DIC,3),sep=""))
#       opt.L <- stored.L
#       break
#     }
#     stored.DIC <- res[[j]]$DIC
#   }
# }

### fit Spectral Semiparametric with L=4 ####
L <- 4 # chosen according to DIC
nb.spectral.semipar <- semipar.eCAR.Leroux(y=covid$Deaths,
                                           x=covid$mean_pm25,
                                           W=adj.matrix,
                                           E=pop.scaled,
                                           C = X.cov,
                                           names.covariates = c(
                                             "p_dens","poverty","house_val","house_inc",
                                             "owner","education","blk","hisp",
                                             "old", "prime","mid",
                                             "d_since_social", "d_since",
                                             "bed","obesity","smoke","s_temp",
                                             "w_temp","s_rm","w_rm"),
                                           pcprior.sd = c(0.1,1), s2=10,
                                           model="Negative Binomial", L=L,
                                           method= "spectral",
                                           verbose=T)

### Spatial Standard or 'naive' ####
nb.naive <- semipar.eCAR.Leroux(y=covid$Deaths,
                                x=covid$mean_pm25,
                                W=adj.matrix,
                                E=pop.scaled,
                                C = X.cov,
                                names.covariates = c(
                                  "p_dens","poverty","house_val","house_inc",
                                  "owner","education","blk","hisp",
                                  "old", "prime","mid",
                                  "d_since_social", "d_since",
                                  "bed","obesity","smoke","s_temp",
                                  "w_temp","s_rm","w_rm"),
                                pcprior.sd = c(0.1,1), s2=10,
                                model="Negative Binomial",
                                method= "naive",
                                verbose=T)

### Spectral Parametric ####
set.seed(101)
nb.spectral.par <- par.eCAR.Leroux(y=covid$Deaths,
                                   x=covid$mean_pm25,
                                   W=adj.matrix,
                                   E=pop.scaled,
                                   C=X.cov,
                                   model="Negative Binomial",
                                   draws=draws, burn=0.5*draws, thin=1,
                                   verbose=TRUE,
                                   joint_prior_lamx_lamz=FALSE, m=0, s2=4)
## note: results in the paper are run longer; draws=10000, burn=5000

## Fig 7 panel (a) Analysis on the entire dataset; ####
pdf("Table_Figures/Main_Fig7a.pdf")
plot(nb.spectral.semipar$beta_omega[,"omega"],
     nb.spectral.semipar$beta_omega[,"beta.mn"], type="l", main = "(a)",
     ylab=expression(exp(beta[k])),
     xlab=expression(paste("", omega[k])), 
     ylim=c(0.85,1.7),  lwd=2, cex.axis=1.5, cex.lab=2.5, cex.main=2, col=3)
lines(nb.spectral.semipar$beta_omega[,"omega"],
      nb.spectral.semipar$beta_omega[,"beta.q025"], col=3, lty=2, lwd=1)
lines(nb.spectral.semipar$beta_omega[,"omega"],
      nb.spectral.semipar$beta_omega[,"beta.q975"], col=3, lty=2, lwd=1)

lines(nb.spectral.par$beta_omega[,"omega"],
      nb.spectral.par$beta_omega[,1], col=1, lwd=2)
lines(nb.spectral.par$beta_omega[,"omega"],
      nb.spectral.par$beta_omega[,2], col=1, lwd=1, lty=2)
lines(nb.spectral.par$beta_omega[,"omega"],
      nb.spectral.par$beta_omega[,3], col=1, lwd=1, lty=2)

lines(nb.naive$beta_omega[,"omega"],
      nb.naive$beta_omega[,"beta.mn"], col=2, lwd=2)
lines(nb.naive$beta_omega[,"omega"],
      nb.naive$beta_omega[,"beta.q025"], col=2, lty=2, lwd=1)
lines(nb.naive$beta_omega[,"omega"],
      nb.naive$beta_omega[,"beta.q975"], col=2, lty=2, lwd=1)
# exp(beta)=1; no effect
abline(h=1, lty=3, lwd=3, col="orange")

legend("topright",
       legend=c("Spectral Semiparametric",
                "Spectral Parametric",
                "Standard Spatial"),
       col=c(3,1,2),
       lty=1, lwd=c(2,2,2),
       bty="n", cex = 1.7)
dev.off()

## Fig 6 PM2.5 exposure and COVID-19 mortality by US county:####
pdf("Table_Figures/Main_Fig6.pdf")
covid$Log_mortality <- log(covid$Deaths/covid$population)
# covid$Log_mortality[!ind] <- NA
library(usmap)
library(ggplot2)
library(viridis)
plot_usmap(regions="county",
           data=covid,
           values="mean_pm25") + scale_fill_viridis()
library(ggthemes)

# fig 6 left panel
plot_usmap(regions="county",
           data=covid,
           values="mean_pm25") + 
  theme_map() +
  theme(legend.position="bottom")+
  scale_fill_viridis()+
  labs(fill = "mean PM2.5")

# fig 6 right panel
plot_usmap(regions="county",
           data=covid,
           values="Log_mortality") + 
  theme_map() +
  theme(legend.position="bottom")+
  scale_fill_viridis()+
  labs(fill = "Log(Deaths/pop)")
dev.off()
