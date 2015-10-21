# Software for Survival Analysis in Patients with Locally Advanced Cervical Cancer

## Introduction
The enclosed code performs analyses presented in the paper ‘Predicting Disease-Free Survival In Locally Advanced Cervical Cancer: A Prospective DCE-MRI Study’, which is currently in review. The code runs survival analysis on a cohort of patients with locally advanced cervix cancer recruited as part of a prospective clinical trial performed at the Christie Hospital NHS Foundation Trust, Manchester, UK. The analysis compares the prognostic value of standard clinicopathologic variables to those obtained from dynamic contrast enhanced imaging data (Tofts, extended Tofts, two-compartment exchange model parameters). The null hypothesis of no difference in predictive accuracy between a model trained using clinicopathologic variables alone and a model trained using the top 6 clinicopathologic and tracer kinetic variables is tested. 

## Other software requirements 

The code is written in the statistical programming language R. The code requires R version 3.1 or higher, which is freely available for Windows and OS X from https://www.r-project.org. There are several ways to run R code, however I recommend use of the RStudio environment, which is also freely available online. The remainder of this README assumes use of RStudio. The code has been tested on both Windows and OS X systems. 

The survival analysis software requires the following R packages from CRAN, which you may need to install yourself (e.g. use RStudio’s Packages tab):

* randomForestSRC
* survival
* rms
* prodlim
* boot
* survcomp
* MASS

The software also requires a package not available on CRAN: **survcomp**. This is available from Bioconductor by typing the following in the R console:

source("http://bioconductor.org/biocLite.R")

biocLite(“survcomp”)

## Installation

To install the software, clone the repository to a suitable location. 

## Running the code

Open the R script ’SurvivalAnalysis.R’ in RStudio. Click the ‘Source’ button in the top right of the main console. The code will take around 10 minutes to run on a standard desktop computer. 

## What the code does

Running ’SurvivalAnalysis.R’ will does the following, in order:

1. Performs Kaplan-Meier analysis on microvascular variables obtained from dynamic contrast enhanced MRI data. For definitions of variables see the paper. 
2. Performs Kaplan-Meier analysis on clinicopathologic variables. 
3. Fits a multivariate random survival forest model to all clinicopathologic and DCE-MRI variables in a bootstrapping framework. Computes point estimates and 95% confidence intervals on median prognostic importance (VIMP). 
4. Defines a null model as one containing only clinicopathologic variables.
5. Selects the top 6 predictors (highest VIMP) for an alternative model. 
6. Performs a leave-one-out analysis on null and alternative models
7. Tests the null hypothesis that the predictive accuracy of the alternative model was no better than the null model. 
8. Makes plots of Figures 3,4,5,6,7 presented in the paper and write pdfs of these into a folder called ‘Figures’. 

## Problems/Bugs/Queries

If you download the software and have trouble making it work, or have a query regarding the paper, 
