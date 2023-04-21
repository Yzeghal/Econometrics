## Econometrics
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

It is part of our StatApp Project

The different libraries used are : 
library(matrixcalc) #for matrix inversion and others
library(MASS) #to generate multivariate data
library(lmtest) #to perform test with specified indications with coeftest
library(sandwich) # to compute heteroscedasticityvariance-covariance.
Use vcovHC(reg, type = "HC0") to have variance covariance matrix calculated with the sandwich formula.
library(aod) #contains functions such that wald test, whiwh we use in 2SLS.
source("C:/Users/tayoy/Documents/GitHub/Econometrics/Basic Regression.R", local = b <- new.env())

Except matrixcalc, none of them is required, since they only serve to build benchmarks and compare results.
