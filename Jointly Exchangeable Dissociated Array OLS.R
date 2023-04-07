# Script for MCO in the case of jointly echangeable and dissociated arrays.
library(matrixcalc) #for matrix inversion and others
library(MASS)#to generate multivariate data
library(lmtest) #to perform test with specified indications with coeftest
library(sandwich) # to compute heteroscedasticity robust tests : vcovHC(reg, type = "HC0")


#generating data
Beta_0 = matrix(c(1,2,4)) #constant coefficient included
n=10
Ui=runif(n,0,1)
Uj=runif(n,0,1)
Uij=matrix(runif(n^2,0,1),ncol=n)


f<-function(Ui,Uj,Uij){
  #function f
  Xij=matrix(c(1,log(Ui),Ui+Uj)) #constant added here
  # print(Xij)
  Yij=t(Beta_0)%*%matrix(Xij)+Uij/10-0.05
  # print(Yij)
  return (array(c(Xij,Yij),c(1,length(Beta_0)+1)))
}
M=array(0,c(n,n,length(Beta_0)+1))

for(i in 1:n){
  for(j in 1:n){
    M[i,j,]=f(Ui[i],Uj[j],Uij[i,j])
  }
}
Xij=M[,,1:length(Beta_0)]
Yij=M[,,length(Beta_0)+1]

# M is a n*n*length(Beta_0)+1 array : 
# Its 2 first dimensions are lines and columns.
# The last dimension is :
# (1, X1,...,Xk,Yk) with length 1+length(Beta_0) because Beta_0 includes
# a constant coef

source(file = "C:/Users/tayoy/Documents/GitHub/Econometrics/Basic Regression.R")
#to recycle function used in Basic Regression.R

x1=as.vector(Xij[,,2])
x2=as.vector(Xij[,,3])
x=matrix(c(x1,x2),ncol=2)
y=matrix(as.vector(Yij[,]))
Beta_hat=coefs_barbarian(x,y)
Beta_hat

