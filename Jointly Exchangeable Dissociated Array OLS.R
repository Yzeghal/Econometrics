# Script for MCO in the case of jointly echangeable and dissociated arrays.
library(matrixcalc) #for matrix inversion and others
library(MASS)#to generate multivariate data
library(lmtest) #to perform test with specified indications with coeftest
library(sandwich) # to compute heteroscedasticity robust tests : vcovHC(reg, type = "HC0")


#generating data
Beta_0 = matrix(c(1,2,4)) #constant coefficient included
n=100
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

#Dispatching NAs to check robustness to missing values.

dispatch_na<-function(Arr,ind=2){
  r=runif(dim(Arr)[1]*dim(Arr)[2],0,1)>0.9
  #print(r)
  frame=Arr[,,ind]
  frame[r]<-NA
  #print(frame)
  Arr[,,ind]=frame
  Arr
}

M=dispatch_na(M, ind=2)
M=dispatch_na(M, ind=3)
M=dispatch_na(M, ind=4)
# for (i in 1:n){
#   M[i,i,]=NA
# }

Xij=M[,,1:length(Beta_0)]
Yij=M[,,length(Beta_0)+1]

# M is a n*n*length(Beta_0)+1 array : 
# Its 2 first dimensions are lines and columns.
# The last dimension is :
# (1, X1,...,Xk,Yk) with length 1+length(Beta_0) because Beta_0 includes
# a constant coef

#transforming data to give them to lm()
x1=as.vector(Xij[,,2])
x2=as.vector(Xij[,,3])
x=matrix(c(x1,x2),ncol=2)
y=matrix(as.vector(Yij[,]))



#Computations 
spot_na<-function(A){
  any(is.na(A))
}
map_na<-function(Arr){
  apply(Arr,c(1,2),spot_na) #applies spot_na to all dimensions having dimensions 1 and 2 fixed.
}
J<-function(X){
  if (dim(X)[1]!=dim(X)[2]){
    stop("Matrix X is not square !")
  }
  not_na=!map_na(X)
  nb_obs=sum(not_na)
  M=matrix(0,dim(X)[3],dim(X)[3])
  
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[1]){
      if (not_na[i,j]){
        M=M+matrix(X[i,j,])%*% t(matrix(X[i,j,]))
      }
      
    }
  }
  M/nb_obs
}
Beta2<-function(X,Y){
  if (any(dim(Y)!=dim(X)[1:2])){
    messageX = paste("2 first dimensions of X", dim(X)[1], dim(X)[2])
    messageY = paste("\n 2 first dimensions of Y", dim(Y)[1], dim(Y)[2])
    stop(paste(messageX,messageY))
  }
  if (dim(X)[1]!=dim(X)[2]){
    stop("Matrix X is not square !")
  }
  not_naX=!map_na(X)
  not_naY=!map_na(Y)
  viable = not_naX&not_naY
  #print(viable)
  if (any (not_naX!=not_naY)){
    messageX = paste("\nnb of values in X:",sum(not_naX))
    messageY = paste("\nnb of values in Y:",sum(not_naY))
    warning(paste(messageX,messageY,"\nMissing values in X and Y are not at the same places !"))
  }
  J_hat=J(X)
  #print(J_hat)
  second_term = matrix(0,dim(X)[3],1)
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[1]){
      if (viable[i,j]){
        second_term=second_term+X[i,j,]*Y[i,j]
      }
    }
  }
  print(paste(sum(viable),'viable terms'))
  second_term=second_term/sum(viable)
  #print(second_term)
  matrix.inverse(J_hat)%*%second_term
}
Beta<-function(X,Y){ #Does the match with lm results (which simply drops all observations) and seems more precise...
  if (any(dim(Y)!=dim(X)[1:2])){
    messageX = paste("2 first dimensions of X", dim(X)[1], dim(X)[2])
    messageY = paste("\n 2 first dimensions of Y", dim(Y)[1], dim(Y)[2])
    stop(paste(messageX,messageY))
  }
  if (dim(X)[1]!=dim(X)[2]){
    stop("Matrix X is not square !")
  }
  not_naX=!map_na(X)
  not_naY=!map_na(Y)
  viable = not_naX&not_naY
  #print(viable)
  if (any (not_naX!=not_naY)){
    messageX = paste("\nnb of values in X:",sum(not_naX))
    messageY = paste("\nnb of values in Y:",sum(not_naY))
    message(paste(messageX,messageY,"\nMissing values in X and Y are not at the same places !"))
  }
  for (i in 1:dim(X)[3]){
    X[,,i]=ifelse(viable ,X[,,i],NA)
  }
  #print(X)
  J_hat=J(X)
  #print(J_hat)
  second_term = matrix(0,dim(X)[3],1)
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[1]){
      if (viable[i,j]){
        second_term=second_term+X[i,j,]*Y[i,j]
      }
    }
  }
  print(paste(sum(viable),'viable terms'))
  second_term=second_term/sum(viable)
  #print(second_term)
  matrix.inverse(J_hat)%*%second_term
}

#Difference with built in
l<-lm(y~x)
Beta_hat=Beta(Xij,Yij)



