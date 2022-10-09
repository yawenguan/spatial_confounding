# Code for paper A spectral adjustment for spatial confounding

This repository contains R-scripts and other files that were employed to carry out simulation studies and data analysis described in the paper "A Spectral Adjustment for Spatial Confounding" 

## First install all required R packages 
	source("load_packages.R")
###  To produce all tables and figures in the manuscript and supplementary materials. 
	source("code_sec5_discretesim_summarize.R")
	source("code_sec5_contsim_summarize.R")
	source("code_sec6_realdata.R")
###  To run discrete and continuous simulation study.  
	source("code_sec5_discretesim.R")
	source("code_sec5_contsim.R")
###  To run data analysis in Section 6. 
	source("code_sec6_realdata.R")

## A brief description of the files follows
Simulation Study in Section 5 of the paper
  
  1. code_sec5_discretesim.R - executes the discrete-case simulation study (~50 min)
	(1). To save time, this script simulates two data sets for each scenario in Table 1.
	(2). For each simulated data set, four methods are compared:
	     	** the standard Leroux, parametric, semiparametric-PCP and semiparametric-R2
	(3). The actual simulation study (500 replicates) is performed on a high performance computer (~48hrs per scenario)
	(4). outputs are saved in folder 'DiscreteSim'
  2. code_sec5_discretesim_summarize.R - organizes output generated from 'code_sec5_discretesim.R'
	* This R script creates Table 1, Fig2 and Fig3 in the manuscript 
	* and Fig1 in the supplementary material 
  3. code_sec5_contsim.R - executes the continuous-case simulation study (~15 min)
	* To save time, this script simulates 1 data set for one scenario in supplementary material Table 1.
	* For each simulated data set, four methods are compared:
	     	** the standard Matern, flexible Matern, parsimonious Matern and semiparametric
		** shorter MCMC are drawn to save time. 
	* The actual simulation study (100 replicates) is performed on a high performance computer (~5hrs per dataset)
	* outputs are saved in folder 'ContSim'
  4. code_sec5_contsim_summarize.R - organizes output generated from 'code_sec5_contsim.R'
	* This R script creates Fig 1 in the manuscript 
	* and Table 1 in the supplementary material 

Data analysis in Section 6 of the paper
  1. code_sec6_realdata.R - R-script that reproduces all output in Section 6
	* This R script creates Figures 4 - 7 in the manuscript 
  2. count_neighborhood_matrix.txt - file that contains neighborhood information of covid example
  3. covid_data2.txt - data file that contains covid data
  4. lip.Rdata - workspace containing lip cancer data

Folder Table_Figures contains table and figures generated from the above
Folder DiscreteSim contains discrete case simulation study outputs
Folder ContSim contains continuous case simulation study output

Other R files
  1. Utilities.R - R-script that contains code that builds the neighborhood matrices.
  2. combine.data.shapefile.R - organize shape file to plot Scotland lip cancer data.
  3. BMcausal_vecc.R - main function to fit parametric model
  4. SemiPcausal_vecc.R - main function to fit semiparametric model
  5. SP_vecc.R - main function to fit standard matern model
  6. MCMC_R2.R and MCMC_GP.R contain utility functions for discrete and continuous methods respectively. 
	
All codes have been tested in Windows 11pro and Linux Ubuntu 20.04 with below information
-----Linux session information-----
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
-----Windows session information-----
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22000)

Matrix products: default
