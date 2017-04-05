# Software for Survival Analysis in Patients with Locally Advanced Cervical Cancer

## Introduction
The enclosed code runs survival analysis on a cohort of cervical cancer patients recruited as part of a prospective clinical trial at the Christie Hospital NHS Foundation Trust, Manchester, UK. The analysis compares the prognostic value of standard clinicopathologic variables to those obtained from dynamic contrast enhanced (DCE) MRI data (two-compartment exchange model parameters, 2CXM). Univariate and multivariate survival analyses are performed to infer the most prognostic parameters for disease-free survival. The null hypothesis of no difference in predictive accuracy between a model trained using clinicopathologic variables alone and a model trained using the top 6 clinicopathologic and 2CXM variables is tested. Results in press in British Journal of Cancer. 

## Other software requirements 

The code is written in the statistical programming language R. The code requires R version 3.1 or higher, which is freely available for Windows and OS X from https://www.r-project.org. There are several ways to run R code, however I recommend use of the RStudio environment, which is also freely available online. The remainder of this README assumes use of RStudio. The code has been tested on both Windows and OS X systems. 

The survival analysis software requires the following R packages from CRAN, which you may need to install yourself (e.g. use RStudio’s Packages tab):

* randomForestSRC
* survival
* rms
* prodlim
* boot
* MASS
* datasets
* grDevices
* colorspace 
* ggplot2

The software also requires a package not available on CRAN: **survcomp**. This is available from Bioconductor by typing the following in the R console:

source("http://bioconductor.org/biocLite.R")

biocLite(“survcomp”)

## Installation

To install the software, clone the repository to a suitable location. 

## Running the code

Open the R script ’SurvivalAnalysis.R’ in RStudio. Click the ‘Source’ button in the top right of the main console. The code will take around 10 minutes to run on a standard desktop computer. 

## What the code does

Running ’SurvivalAnalysis.R’ will does the following, in order:

1. Generates Kaplan-Meier curves for all 2CXM clinicopathologic variables.  
2. Fits a multivariate random survival forest model to all clinicopathologic and 2CXM variables in a bootstrapping framework. Computes point estimates and 95% confidence intervals on median variable importance (VIMP). 
3. Defines a null model as one containing only clinicopathologic variables.
4. Selects the top 6 predictors (highest VIMP) for an alternative model. 
5. Performs a leave-one-patient-out analysis on null and alternative models, computing the accuracy of models on the left out patient
6. Tests the null hypothesis of no difference in predictive accuracy between null and alternative models.
7. Make a plot showing DFS risk predictions for each variable in the alternative model, adjusted for all other variables in the model.  
## Problems/Bugs/Queries

If you download the software and have trouble making it work, or have a query regarding the paper, please email: ben.dickie@postgrad.manchester.ac.uk
