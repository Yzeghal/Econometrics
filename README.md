# Econometrics
This repository contains codes to run a regression adapted to jointly exchangeable and dissociated arrays. 
Everything from heteroscedastic robust OLS to JEDA 2SLS was recoded "by hand" with matrix formulas.
Most of the code is vectorised,, but some function remain built on nested for loops, but a set of 10  000 observation does not make the difference tangible. Wald test was taken form aod library, but I checked it performed the exact matrix calculation (see sandbow script)

There are 2 main scripts :
Jointly Exchangeable Dissociated Array OLS.R

and

Sandbox_tester.R

You will find all the functions required to run a JEDA-2SLS in the first one.
The second one contains many tests of the functions, but also illustration.
I mainly used it to generate 2d JEDAS with endogeneity, control variables and heteroscedasticity, using Aldous(1981) representation.
Basic Regression.R contains code that are not directly useful, but work well and are used as benchmarks in the sandbox.
It is part of our StatApp Project

The different libraries used are : 
library(matrixcalc) #for matrix inversion and others
library(aod) #contains functions such that wald test we use in 2SLS. Homemade version is in Sandbox_tester
library(AER) #for a 2SLS built-in estimator (much better than homemade one in case of NAs)
library (R.utils) #to wrap array if we use built-in methods instead of homemade ones

## How to load it 
To use OLS and IV_LS on your own data you will need to use the functions that are in Jointly Exchangeable Dissociated Array.R
To load it, type :
source("Script Location") #with the path to Jointly Exchangeable Dissociated Array.R in "script location"

## How to use IV_LS :
All data inputs should be in n*n*d arrays (with n in common). the first 2 dimensions contain the index of observations, and vectors valued regressors are along the 3rd one.
inputs :
G : Control variables, containing only 1 in the first index along the 3rd dimension
(G[,,1] contains 1 or NAs)
Z : Instruments, (without including G)
X : Endogenous variables
Y : Regressed variable. 

Hyp : (default at 0) vector of hyopothesis to test (must be the same length as the estimated beta_2SLS.

built_in_reg : (default at TRUE) TRUE to use R lm() and AER methods, FALSE for homemade method. They are equivalent if no NAs are in the data, but might be significantly different if NAs are too numerous.

Outputs : 
list('SLS'=list('coefs'=values, 'var'=asvar, 'F'=F_test,"R2"= R2),
      'FLS'=list('R2'=reg1$R2,'F'=reg1$F)))
A nested list containing informations on both regressions. 
The Beta vector contains the control coefficients (starting with 1) and X coefficients, in the same order as given.
F statistics are built with the JEDA estimated asymptotic variance in the 1st and 2nd stage.

Unless you are sure not to have NAs out of your diagonal set (observations with the same index), please use all parameters as set to default. 
## How to use OLS :
All data inputs should be in n*n*d arrays (with n in common). the first 2 dimensions contain the index of observations, and vectors valued regressors are along the 3rd one.

inputs :
X : Explicative variables, only 1 in the first index along the 3rd dimension

Y : Regressed variable.

model : "JEDA"    -> (default method)will compute the asymptotic variance using our method
        "BASIC"   -> will compute heteroscedastic asymptotic variance
        "NO ASVAR"-> will compute nothing but the Beta_hat OLS estimator
        
Hyp : (default at 0) vector of hyopothesis to test (must be the same length as the estimated beta_hat.

built_in_reg : TRUE to use R lm() method, FALSE for homemade method. They are equivalent if no NAs are in the data, but might be significantly different if NAs are too many.

Outputs : 
list('SLS'=list('coefs'=values, 'var'=asvar, 'F'=F_test,"R2"= R2),
      'FLS'=list('R2'=reg1$R2,'F'=reg1$F)))
A nested list containing informations on both regressions.
F statistics are built with the JEDA estimated asymptotic variance in the 1st and 2nd stage.

Have fun !

Y.

