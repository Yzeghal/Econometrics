source('C:/Users/tayoy/Documents/GitHub/Econometrics/Jointly Exchangeable Dissociated Array OLS.R')
source("C:/Users/tayoy/Documents/GitHub/Econometrics/Basic Regression.R", local = b <- new.env())
#If you want to compare with built-in functions :
library(MASS)#to generate multivariate data
library(lmtest) #to perform test with specified indications with coeftest
library(sandwich) # to compute heteroscedasticityvariance-covariance.
# Use vcovHC(reg, type = "HC0") to have variance covariance matrix calculated with the sandwich formula.
library(AER) #for a 2SLS built-in functions benchmark 

#Case of OLS -----

dispatch_na<-function(Arr,ind=2){
  #Dispatches NAs to check robustness to missing values.
  r=runif(dim(Arr)[1]*dim(Arr)[2],0,1)>0.9
  #print(r)
  frame=Arr[,,ind]
  frame[r]<-NA
  #print(frame)
  Arr[,,ind]=frame
  Arr
}
diag_na<-function(Arr){
  #Puts NAs in diagonal so that those value won't count
  d=dim(Arr)[1]
  ind=(1:d-1)*(d+1)+1
  m=matrix(Arr[,,2])
  m[ind]<-NA
  Arr[,,2]=m
  return(Arr)
}

Beta_0 = matrix(c(1,2,4))
n=50
Ui=runif(n,0,1)
Uj=runif(n,0,1)
Uij=matrix(runif(n^2,0,1),ncol=n)

f<-function(UiUjUij){
  #function f
  Ui=UiUjUij[1]
  Uj=UiUjUij[2]
  Uij=UiUjUij[3]
  c=10
  e=0.1
  Xij=matrix(c(1,Ui,Uj)) #constant added here
  # print(Xij)
  epsilon=-(Ui+Uj)*c/2+Ui*Uj*c+c/4-e/2+e*Uij
  Yij=t(Beta_0)%*%matrix(Xij)+epsilon
  # print(Yij)
  return (array(c(Xij,Yij),c(1,length(Beta_0)+1)))
}
N = array(0,c(n,n,3))
N[,,1] = matrix(rep(Ui,n),c(n,n), byrow=FALSE) 
N[,,2] = matrix(rep(Uj,n),c(n,n), byrow=TRUE)
N[,,3] = Uij
M = apply(N,MARGIN = c(1,2),FUN = f)
M=aperm(M,c(2,3,1))

# data generation
#----
M=array(0,c(n,n,3))
M[,,1] = matrix(rep(Ui,n),c(n,n), byrow=FALSE) 
M[,,2] = matrix(rep(Uj,n),c(n,n), byrow=TRUE)
M[,,3] = Uij
A = apply(M,MARGIN = c(1,2),FUN = g2)
A=aperm(A,c(2,3,1))
#A=dispatch_na(A)
#A=diag_na(A)
Y = A[,,1]
G = A[,,2:3]
Z = A[,,4:5]
X = A[,,6:7]
#----
e1=as.vector(A[,,8])
e2=as.vector(A[,,9])
e3=as.vector(A[,,10])
x1=as.vector(X[,,1])
x2=as.vector(X[,,2])
z1=as.vector(Z[,,1])
z2=as.vector(Z[,,2])
g1=as.vector(G[,,2])
cor(e3,z1)
cor(e3,z2)
reg1=FS_OLS(G,Z,X) #regression of X on G and Z
reg1
B = reg1$Betas
D=predIV(G,Z,X,first_Betas=B)
G=array(G,dim=c(dim(G),1))
GD=array(c(G,D),dim=c(dim(G)[1:2],dim(D)[3]+dim(G)[3]))
GX=array(c(G,X),dim=c(dim(G)[1:2],dim(X)[3]+dim(G)[3])) #to compare with X

reg2=OLS(GD,Y,model="JEDA") #regression of Y on G and D
noIV<-OLS(GX,Y,model="BASIC")

reg2
noIV



# M is a n*n*length(Beta_0)+1 array :
# Its 2 first dimensions are lines and columns.
# The last dimension is :
# (1, X1,...,Xk,Yk) with length 1+length(Beta_0) because Beta_0 includes
# a constant coef

# M=dispatch_na(M, ind=2)
# M=dispatch_na(M, ind=3)
# M=dispatch_na(M, ind=4)

# M=diag_na(M)
Xij=M[,,1:length(Beta_0)]
Yij=M[,,length(Beta_0)+1]
#transforming data to give them to lm()
x1=as.vector(Xij[,,2])
x2=as.vector(Xij[,,3])
x=matrix(c(x1,x2),ncol=2)
y=matrix(as.vector(Yij[,]))

#Difference with built in
reg<-lm(y~x)
Beta_hat=reg$coefficients
tests=coeftest(reg, vcov = vcovHC(reg, type = "HC0")) #performs tests on the model reg with the variance matrix calculated with vcovHC


JEDA<-OLS(Xij,Yij, model="JEDA")
LS<-OLS(Xij,Yij, model="BASIC")

JEDA
LS
tests
F_=wald.test(vcovHC(reg, type="HC0"), b=Beta_hat,df=3, L=diag(3))$result$Ftest[1]
F_



#Case of 2SLS no control variable----
Beta_0 = matrix(c(1,2,3)) #constant coefficient included :(cst, X1,X2)
Beta_1 = matrix(c(10,5,7)) #(cst,Z1, Z2) #boost cst coef to max mean(x) and the endogeneity error (on cst coef mainly)
Beta_2 = matrix(c(20,4,8)) #(cst,Z1, Z2)
n=100
Ui=runif(n,0,1)
Uj=runif(n,0,1)
Uij=matrix(runif(n^2,0,1),ncol=n)

g<-function(UiUjUij){
  #UiUjUij is a vector c(Ui,Uj,Uij)
  #function that generates the X, IVs Z and G s.th G=1 for now
  U1=UiUjUij[1]
  U2=UiUjUij[2]
  U12=UiUjUij[3]
  p<-2*pi
  N1=cos(p*U12) #orthogonal to all other random variables
  N2=cos(3*p*U12)
  eps1 = 1*(U1*U2/3+0.66)*N1 
  eps2 = 1*(U1*U2/3+0.66)*N2
  eps3 = N1+N2 #Correlated to eps1, eps2 but not Z. Also heteroscedastic
  Z1=U1
  Z2=U2
  X1=t(Beta_1)%*%matrix(c(1,Z1,Z2))+eps1
  X2=t(Beta_2)%*%matrix(c(1,Z1,Z2))+eps2
  Y =t(Beta_0)%*%c(1,X1,X2)+eps3 #Y explained by X with endogeneity and heteroscedasticity
  r=c(Y,1,Z1,Z2,X1,X2,eps1,eps2,eps3)
  return (r)
} 


# data generation
#----
M=array(0,c(n,n,3))
M[,,1] = matrix(rep(Ui,n),c(n,n), byrow=FALSE) 
M[,,2] = matrix(rep(Uj,n),c(n,n), byrow=TRUE)
M[,,3] = Uij
A = apply(M,MARGIN = c(1,2),FUN = g)
A=aperm(A,c(2,3,1))

Y = A[,,1]
G = A[,,2]
Z = A[,,3:4]
X = A[,,5:6]
#----
e1=as.vector(A[,,7])
e2=as.vector(A[,,8])
e3=as.vector(A[,,9])
x1=as.vector(X[,,1])
x2=as.vector(X[,,2])
z1=as.vector(Z[,,1])
z2=as.vector(Z[,,2])
cor(e1,z2)
cor(e2,x2)
reg1=FS_OLS(G,Z,X) #regression of X on G and Z
reg1
B = reg1$Betas
D=predIV(G,Z,X,first_Betas=B)
G=array(G,dim=c(dim(G),1))
GD=array(c(G,D),dim=c(dim(G)[1:2],dim(D)[3]+dim(G)[3]))
GX=array(c(G,X),dim=c(dim(G)[1:2],dim(X)[3]+dim(G)[3])) #to compare with X

reg2=OLS(GD,Y,model="BASIC") #regression of Y on G and D. Basic is for calculation speed.
noIV<-OLS(GX,Y,model="BASIC")

reg2
noIV

#Case of 2SLS with control variable----
source('C:/Users/tayoy/Documents/GitHub/Econometrics/Jointly Exchangeable Dissociated Array OLS.R')

Beta_0 = matrix(c(1,2,3,1)) #constant coefficient included :(cst,G,X1,X2)
Beta_1 = matrix(c(10,10,5,7)) #(cst,Z1, Z2) #boost cst coef to max mean(x) and the endogeneity error (on cst coef mainly)
Beta_2 = matrix(c(10,10,4,8)) #(cst,Z1, Z2)
n=100
Ui=runif(n,0,1)
Uj=runif(n,0,1)
Uij=matrix(runif(n^2,0,1),ncol=n)

g2<-function(UiUjUij){
  #UiUjUij is a vector c(Ui,Uj,Uij)
  #function that generates the X, IVs Z and G s.th G=1 for now
  U1=UiUjUij[1]
  U2=UiUjUij[2]
  U12=UiUjUij[3]
  p<-2*pi
  N1=cos(p*U12) #orthogonal to all other random variables
  N2=cos(3*p*U12)
  G = U1+3*cos(7*p*U12) #orthogonal to U2, U12 and their cos(2pi n . ) for n !=7 but not U1 
  eps1 = 30*(U1*U2/3+0.66)*N1 
  eps2 = 30*(U1*U2/3+0.66)*N2
  eps3 = N1+N2 #Correlated to eps1, eps2 but not Z. Also heteroscedastic
  Z1=U1
  Z2=U2
  X1=t(Beta_1)%*%matrix(c(1,G,Z1,Z2))+eps1 #correlate G and X to have impact on G coef in OLS but not so much in 2SLS
  X2=t(Beta_2)%*%matrix(c(1,G,Z1,Z2))+eps2
  Y =t(Beta_0)%*%c(1,G,X1,X2)+eps3 #Y explained by X with endogeneity and heteroscedasticity
  r=c(Y,1,G,Z1,Z2,X1,X2,eps1,eps2,eps3)
  return (r)
}


# data generation
#----
M=array(0,c(n,n,3))
M[,,1] = matrix(rep(Ui,n),c(n,n), byrow=FALSE) 
M[,,2] = matrix(rep(Ui,n),c(n,n), byrow=TRUE) #test with the same i and not i,j
M[,,3] = Uij
A = apply(M,MARGIN = c(1,2),FUN = g2)
A=aperm(A,c(2,3,1))
for (i in 1:7){
  A=dispatch_na(A,i)
}

A=diag_na(A)
Y = A[,,1]
G = A[,,2:3]
Z = A[,,4:5]
X = A[,,6:7]
#----
e1=as.vector(A[,,8])
e2=as.vector(A[,,9])
e3=as.vector(A[,,10])
x1=as.vector(X[,,1])
x2=as.vector(X[,,2])
z1=as.vector(Z[,,1])
z2=as.vector(Z[,,2])
k1=as.vector(G[,,2])
y =as.vector(Y)
cor(e3,x1)
cor(e3,x2)
cor(e3,z1)
cor(e3,z2)
sd(e3)
reg1=FS_OLS(G,Z,X,built_in_reg = FALSE ) #regression of X on G and Z
B2 = reg1$Betas
D=predIV(G,Z,X,first_Betas=B2)
GD=array(c(G,D),dim=c(dim(G)[1:2],dim(D)[3]+dim(G)[3]))
GX=array(c(G,X),dim=c(dim(G)[1:2],dim(X)[3]+dim(G)[3])) #to compare with X
GZ=array(c(G,Z),dim=c(dim(G)[1:2],dim(Z)[3]+dim(G)[3]))

reg2_built_out=OLS(GD,Y,model="JEDA", built_in_reg = FALSE) #regression of Y on G and D
reg2_built_in=OLS(GD,Y,model="JEDA", built_in_reg = TRUE) #regression of Y on G and D
noIV<-OLS(GX,Y,model="BASIC")
reg3_built_in<-IV_LS(G,Z,X,Y,built_in_reg=TRUE)
reg3_built_out2<-IV_LS(G,Z,X,Y,built_in_reg=FALSE)
noIV
#Check wald.test 
#---- 
Beta_hat = matrix(reg2_built_out$coefs[,1])
asvar = reg2_built_out$var
F_= 100 * t(Beta_hat)%*%matrix.inverse(asvar)%*%Beta_hat 
W=wald.test(Sigma=asvar/100, b=Beta_hat, L=diag(length(Beta_hat)),verbose=TRUE)$result$chi2[1]
W-F_ #~1e-10
#----

#comparaison with built-in 2SLS
#----

iv<-ivreg(formula = y ~k1 + x1 + x2  | k1+ z1 + z2 )

reg3_built_in$SLS$coefs[,1]-iv$coefficients #1e-13 even with nas
reg3_built_out$SLS$coefs[,1]-iv$coefficients #1e-10 if no Nas


cz=cbind(g1,z1,z2)
GoodBeta=b$reg_OLS(cz,matrix(x1))$coefs[,1]
BetaOls2=Beta(GZ,X[,,1])
J(GZ,X[,,1])

