library(matrixcalc) #for matrix inversion and others
library(MASS)#to generate multivariate data
library(lmtest) #to perform test with specified indications with coeftest
library(sandwich) # to compute heteroscedasticity robust tests : vcovHC(reg, type = "HC0")
#library(lmtest)
#library(sandwich)
generator<-function(dim1=2, dim2=2, range=100){
  #to generate random matrices and test functions
  l=dim1*dim2
  v=sample(-100:100,size=l,replace=TRUE)
  v=matrix(v,dim1,dim2)
}
  
#generating data
#set.seed(0)
eps=rnorm(1000, mean=0, sd=0.2)
l=2
d=diag(l)
d2=cbind(d,matrix(rep(0,l)))[,l:(l+1)]
d<-d+d2
sigma<-(d+t(d))/2

beta= matrix(c(2,4))
x=mvrnorm(1000,c(100,50),sigma[1:l,1:l])
y=x%*%beta+eps



coefs<-function(X,Y){
  nx=dim(X)[1]
  ny=dim(y)[1]
  if (nx!=ny){stop('not the same number of observations')}
  M=t(X)%*%X
  X=cbind(rep(1,nx),X)
  if (det(M)==0){stop('Matrix X\'X non inversible !')}
  matrix.inverse(t(X)%*%X)%*%(t(X)%*%Y)
}

estimate<-function(X,Beta_hat){
  nx=dim(X)[1]
  X=cbind(matrix(rep(1),nx),X)
  Y_hat=X%*%Beta_hat
}

asymptotic_variance<-function(X,Y){
  nx=dim(X)[1]
  beta_hat=coefs(X,Y)
  eps_hat=Y-estimate(X,beta_hat)
  d_eps=diag(as.vector(eps_hat))
  print(dim(d_eps))
  X=cbind(rep(1,nx),X)
  X_eps = d_eps%*%X
  print(dim(X_eps))
  middle_term = t(X_eps)%*%(X_eps)*nx
  print(dim(middle_term))
  sandwich_term= matrix.inverse(t(X)%*%X)
  print(dim(sandwich_term))
  var = sandwich_term %*% middle_term %*% sandwich_term
}

standard_errors<-function (var,nx){
  s=matrix(sqrt(diag(var)/nx))
}

student_t<-function(Beta_hat,var,nx,hyp = 0){
  sd=sqrt(diag(var)/nx)
  t=(Beta_hat-hyp)/sd
}
#homemade model fitting
Beta_hat = coefs(x,y)
asvar=asymptotic_variance(x,y)
se=standard_errors(asvar, 1000)

#built-in model fitting
reg=lm(y~x)
reg=lm(y~x)
tests=coeftest(reg, vcov = vcovHC(reg, type = "HC0")) #performs tests on the model reg with the variance matrix calculated with vcovHC

#comparaison
print(Beta_hat-tests[,1]) #difference between built-in and homemade estimators
print(se-tests[,2]) #difference between built-in and homemade std error

#difference are around 1e-10 with Hc0 and 1e-4 with Hc3 

