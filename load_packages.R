# package required for running the code
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages("RandomFieldsUtils")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/RandomFields/RandomFields_3.3.14.tar.gz"
install.packages(packageurl, repos=NULL)

load_packages <- c("GpGp","ggplot2", "broom", "dplyr", "viridis","usmap","fields","RandomFieldsUtils",
                   "splines","mvtnorm","eCAR","geoR","ggthemes","truncnorm","mvnfast","maptools","sp")                                        # Specify your packages
not_installed <- load_packages[!(load_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)   
