# Script for MCO in the case of jointly echangeable and dissociated arrays.
library(matrixcalc) #for matrix inversion and others
library(MASS)#to generate multivariate data
library(lmtest) #to perform test with specified indications with coeftest
library(sandwich) # to compute heteroscedasticityvariance-covariance.
# Use vcovHC(reg, type = "HC0") to have variance covariance matrix calculated with the sandwich formula.
library(aod) #contains functions such that wald test, whiwh we use in 2SLS.
#source("C:/Users/tayoy/Documents/GitHub/Econometrics/Basic Regression.R", local = b <- new.env())

#Computations ----
spot_na<-function(A){
  any(is.na(A))
}
map_na<-function(Arr){
  apply(Arr,c(1,2),spot_na) #applies spot_na to all dimensions having dimensions 1 and 2 fixed.
}

J<-function(X,Y){ #returnns J NORMALISED by nuber of observations used in the sum
  if (dim(X)[1]!=dim(X)[2]){
    stop("Matrix X is not square !")
  }
  not_nax=!map_na(X)
  not_nay=!map_na(Y)
  not_na=not_nax&not_nay
  nb_obs=sum(not_na)
  M=matrix(0,dim(X)[3],dim(X)[3])
  
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[1]){
      if (not_na[i,j]){
        M=M+matrix(X[i,j,])%*%t(matrix(X[i,j,]))
      }
      
    }
  }
  M=M/nb_obs
  return (M)
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
} #Not to be used. Was expected to be better than Beta, but is not.

Beta<-function(X,Y, J_hat=NULL){ #Does the match with lm results (which simply drops all observations) and seems more precise...
  #add H_hat in argument if it was already calculated.
  #Safety tests
  #----
  
  if (any(dim(Y)!=dim(X)[1:2])){
    messageX = paste("2 first dimensions of X", dim(X)[1], dim(X)[2])
    messageY = paste("\n 2 first dimensions of Y", dim(Y)[1], dim(Y)[2])
    stop(paste(messageX,messageY))
  }
  if (dim(X)[1]!=dim(X)[2]){
    stop("Matrix X is not square !")
  }
  #----
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
  if (is.matrix(J_hat)&&is.square.matrix(J_hat)&&dim(J_hat)[1]==dim(X)[3]){
      print('J given in argument is used')
    }
  else{
    J_hat=J(X,Y)
  }
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

estimate <-function(X,Beta_hat){
  apply(X,MARGIN=c(1,2), FUN =(function(x){t(x)%*%Beta_hat}))
}

eps<-function(X,Y,b){
  ep=Y-estimate(X,b)
  return (ep)
}

last_multip<-function(Arr, d=-1){
  Arr[1:d-1]*Arr[d]
}

H_basic_OLS<-function(X,Y,b){
  #returns a list(H, moments) that contains H and moments that should be equal to 0
  #Safety tests
  #----
  if (any(dim(Y)!=dim(X)[1:2])){
    messageX = paste("2 first dimensions of X", dim(X)[1], dim(X)[2])
    messageY = paste("\n 2 first dimensions of Y", dim(Y)[1], dim(Y)[2])
    stop(paste(messageX,messageY))
  }
  if (dim(X)[1]!=dim(X)[2]){
    stop("Matrix X is not square !")
  }
  d=dim(X)[3]
  nx=dim(X)[1]
  if (length(b)!=d){
    stop("b and dim(X)[3] should have the same length")
  }
  #----
  ep=eps(X,Y,b)
  not_naX=!map_na(X)
  not_naY=!map_na(Y)
  viable = not_naX&not_naY #ep is already NA if X or Y is NA
  #tests that diagonal is NA
  if (sum(diag(viable))!=0){
    warn=paste(sum(diag(viable)),"terms in the diagonal. Should be 0")
    warning(warn)
  }
  #----
  Xeps = array(0,c(nx,nx,d+1))
  Xeps[,,1:d]=X
  Xeps[,,d+1]=ep
  Xeps=apply(Xeps,MARGIN=c(1,2), last_multip, d=d+1)
  Xeps=aperm(Xeps,c(2,3,1))
  #all(!map_na(Xeps)==viable) #test if all NAs are at the right place
  #moment = apply(Xeps,MARGIN= 3, FUN= sum, na.rm=TRUE) # Moment condition of least squares ~1e-12 to have an idea of calculation precision loss.
  not_zeros = sum(viable) #counts the number of not calculated terms in the sum over i.
  #----
  eXXe = apply(Xeps,MARGIN = c(1,2), FUN=function(x){
    m=ifelse(any(is.na(x))*matrix(1,d,d),matrix(NA,d,d),matrix(x)%*%t(x))
    # yields d*d matrices taken as d^2 length vectors. Does not change anything for the sum and is re-matrixed after.
    return (m)})
  print(dim(eXXe))
  #print(all(!map_na(eXXe[1,,])==viable)) #tests that all NA values are placed correctly for index 1 of matrices.
  mat=apply(eXXe,MARGIN = 1, FUN = sum, na.rm=TRUE) #is a d^2 vector
  mat=matrix(mat, ncol = d)/not_zeros #normalisation by nb of valid observations
  print ('nb obs :')
  print (not_zeros)
  return(mat)
 #---- 
  #heavy unvectorized calculations
  # for (i in 1:nx){ #Sum over i
  #   count=0 #counts the number of viable j terms
  #   V=rep(0,d) #initialize the sum of X_ij*eps_hat_ij
  #   for (j in 1:nx){
  #     if (viable[i,j]&&viable[j,i]){
  #       count=count+1
  #       V=V+Xeps[i,j,]+Xeps[j,i,]
  #     }
  #   }
  #   if (count==0){
  #     nb_zeros=nb_zeros+1
  #   }
  #   count = max(count,1)
  #   S=S+(matrix(V)%*%t(V)/count^2)
  # }
}


H_JEDA<-function(X,Y,b){
  
  #Safety tests
  #----
  if (any(dim(Y)!=dim(X)[1:2])){
    messageX = paste("2 first dimensions of X", dim(X)[1], dim(X)[2])
    messageY = paste("\n 2 first dimensions of Y", dim(Y)[1], dim(Y)[2])
    stop(paste(messageX,messageY))
  }
  if (dim(X)[1]!=dim(X)[2]){
    stop("Matrix X is not square !")
  }
  d=dim(X)[3]
  nx=dim(X)[1]
  if (length(b)!=d){
    stop("b and dim(X)[3] should have the same length")
  }
  #----
  ep=eps(X,Y,b)
  not_naX=!map_na(X)
  not_naY=!map_na(Y)
  viable = not_naX&not_naY #ep is already NA if X or Y is NA
  #tests that diagonal is NA
  if (sum(diag(viable))!=0){
    warn=paste(sum(diag(viable)),"terms in the diagonal. Should be 0")
    warning(warn)
  }
  Xeps = array(0,c(nx,nx,d+1))
  Xeps[,,1:d]=X
  Xeps[,,d+1]=ep
  Xeps=apply(Xeps,MARGIN=c(1,2), last_multip, d=d+1)
  Xeps=aperm(Xeps,c(2,3,1))
  #all(!map_na(Xeps)==viable) #test if all NAs are at the right place
  #moment = apply(Xeps,MARGIN= 3, FUN= sum, na.rm=TRUE) # Moment condition of least squares ~1e-12 to have an idea of calculation precision loss.
  nb_zeros = 0 #counts the number of not calculated terms in the sum over i.
  S=matrix(0,d,d)
  #heavy unvectorized calculations
  for (i in 1:nx){ #Sum over i
    count=0 #counts the number of viable j terms
    V=rep(0,d) #initialize the sum of X_ij*eps_hat_ij+X_ji*eps_hat_ji
    for (j in 1:nx){
      if (viable[i,j]&&viable[j,i]){
        count=count+1
        V=V+Xeps[i,j,]+Xeps[j,i,]
      }
    }
    if (count==0){
      nb_zeros=nb_zeros+1
    }
    count = max(count,1)
    S=S+(matrix(V)%*%t(V)/count^2)
  }
  
  S=S/(nx-nb_zeros)
  return (S)
}

#Small functions for inference
#----


asymptotic_variance_OLS<-function(J,H){
  invJ=matrix.inverse(J)
  return(invJ %*% H %*% invJ)
}

standard_errors<-function (var,nx){
  s=matrix(sqrt(diag(var)/nx))
}

student_t<-function(Beta_hat,sd,hyp = 0){
  t=(Beta_hat-hyp)/sd
  return(t)
}

p_values<-function(Beta_hat, sd ,nx , hyp = 0){
  deg_freedom = nx-length(Beta_hat)
  t=abs(student_t(Beta_hat, sd, hyp))
  p=2*(1-(pt(t,deg_freedom)))
  return(p)
}

#----

OLS<-function(X,Y,hyp=0,model="JEDA"){
  #safety tests
  #----
  if (!any(model==(c('JEDA','BASIC')))){
    stop("Model not recognised. Try \'JEDA\' or \'BASIC\' instead.")
  }
  if (!is.array(X)){
    stop("X should be an array")
  }
  if (!is.matrix(Y)){
    stop("Y should be a matrix or a 2d array")
  }
  dimx=dim(X)
  dimy=dim(Y)
  if (any(dimx[1:2]!=dimy)){
    message=paste("\ndimension of X :", dimx, "\nimension of Y :", dimy)
    stop(message)
  }
  if (!is.square.matrix(Y)){
    stop('Y and X should be square ! check dim(X) and dim(Y)')
  }
  #----
  not_naX=!map_na(X)
  not_naY=!map_na(Y)
  viable = not_naX&not_naY
  nb_obs = sum(viable)
  message("Start Beta Calculation...")
  J_hat=J(X,Y)
  if (det(J_hat)==0){
    print('J:')
    print (J_hat)
    stop ("J is singular !")
  }
  Beta_hat = Beta(X,Y, J=J_hat)
  if (model == "JEDA"){
    H_hat=H_JEDA(X,Y,Beta_hat)
  }
  if (model == "BASIC"){
    H_hat=H_basic_OLS(X,Y,Beta_hat)
  }
  asvar=asymptotic_variance_OLS(J=J_hat, H=H_hat)
  std=standard_errors(asvar,nb_obs)
  ts=student_t(Beta_hat,std,hyp = hyp)
  p_val=p_values(Beta_hat, std ,nb_obs , hyp = 0)
  if (sum(abs(hyp))==0){hyp = rep(0,length(Beta_hat))}
  values = matrix(c(Beta_hat,std,ts,p_val,hyp), ncol=5)
  colnames(values)=c('Beta_hat', 'Std_Err','Student_t','p-values', 'H0_hyp')
  F_test = wald.test(Sigma=asvar/nb_obs, b=Beta_hat ,df=length(Beta_hat), L=diag(length(Beta_hat)))$result$Ftest[1]
  Y_hat = estimate(X,Beta_hat)
  ep = eps(X,Y,Beta_hat)
  my=mean(Y, na.rm=TRUE)

  SCE = sum((Y_hat-my)^2, na.rm=TRUE)
  SCT = sum((Y-my)^2, na.rm = TRUE)
  R2<-SCE/SCT
  return(list ('coefs'=values, 'var'=asvar, 'F'=F_test,"R2"= R2))
}

#Data generation----
#generating data as jointly exchangeable dissociated arrays.

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


#Case of OLS -----

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



#Case of 2SLS----

Betas_1SLS<-function(G,Z,X){
  #----
  if (!(dim(X)[1:2]==dim(Z)[1:2]&&dim(Z)[1:2]==dim(G)[1:2])){
    stop('\nFirst 2 dimensions of X,Z and G should be the same')
  }

  if (!dim(X)[1]==dim(X)[2]){
    stop('Arrays should be square in their first 2 dimensions ')
  }
  m=matrix(1,dim(G)[1],dim(G)[2])
  if(any(G[,,1]!=m)){
    stop('Array G should contain the vector 1 in first position')
  }
  #----
  d=dim(X)[3]
  GZ=array(0,c(dim(X)[1:2],dim(Z)[3]+dim(G)[3]))
  GZ[,,1:dim(G)[3]]=G
  GZ[,,(1+dim(G)[3]):(dim(G)[3]+dim(Z)[3])]=Z
  
  Betas = matrix(0,dim(Z)[3]+dim(G)[3],d)
  J_hat=J(GZ)
  for (i in 1:d){
    print(paste('Regression on IV and G number :',i ))
    beta_i=Beta(Y=X[,,i],X=GZ, J_hat = J_hat)
    Betas[,i]=beta_i
  }
  return (Betas)
}


D_SLS<-function(G,Z,X,first_Betas=NULL){
  #----
  if (!is.null(first_Betas)){
    if(is.matrix(first_Betas)){
      if (!(all((dim(first_Betas)==c(dim(G)[3]+dim(Z)[3],dim(X)[3]))))){
        mess=paste("first_Betas given has dimension : ", dim(first_Betas)[1],dim(first_Betas)[2],"\n function expected", dim(G)[3]+dim(Z)[3],dim(X)[3])
        stop(mess)
      }
    }
    else {stop('first_Betas is not a matrix !')}
    
  }
  if (!(dim(X)[1:2]==dim(Z)[1:2]&&dim(Z)[1:2]==dim(G)[1:2])){
    stop('\nFirst 2 dimensions of X,Z and G should be the same')
  }

  if (!dim(X)[1]==dim(X)[2]){
    stop('Arrays should be square in their first 2 dimensions ')
  }
  if (!dim(X)[1]==dim(X)[2]){
    stop('Arrays should be square in their first 2 dimensions ')
  }
  m=matrix(1,dim(G)[1],dim(G)[2])
  if(any(G[,,1]!=m)){
    stop('Array G should contain the vector 1 in first position')
  }
  m=matrix(1,dim(G)[1],dim(G)[2])
  #----
  if (is.null(first_Betas)){
    first_Betas = Betas_1SLS(Z,G,X)
  }
  GZ = array(0,dim=c(dim(Z)[1:2],dim(Z)[3]+dim(G)[3]))
  GZ[,,1:dim(G)[3]]=G
  GZ[,,(1+dim(G)[3]):(dim(G)[3]+dim(Z)[3])]=Z
  D=array(0, dim(X))
  for(i in 1:dim(X)[3]){
    b=first_Betas[,i]
    print('b :')
    print(b)
    D[,,i] = apply(GZ,MARGIN=c(1,2), FUN =(function(x){t(x)%*%b}))
    print('sd(X) :')
    print(sd(X[,,i]))
    print('sd(D) :')
    print(sd(D[,,i]))
    
    print(sd(D[,,i],na.rm=TRUE)/sd(X[,,i],na.rm=TRUE))
  }
  return (D)
}

Beta_2SLS<-function(D,G,Y){
  #----
  if (any(dim(G)[1:2]!=dim(D)[1:2])){
    stop("First 2 dimensions of G and D should be equal")
  }
  #----
  GD=array(0,c(dim(G)[1:2],dim(G)[3]+dim(D)[3]))
  GD[,,1:dim(G)[3]]=G
  GD[,,(1+dim(G)[3]):(dim(G)[3]+dim(D)[3])]=D
  B=Beta(GD,Y)
  return (B)
}


Beta_0 = matrix(c(1,2,4,8)) #constant coefficient included :(cst, X1,X2,G)
Beta_1 = matrix(c(3,5,7)) #(cst,Z1, Z2)
Beta_2 = matrix(c(0.2,0.1,0.1)) #(cst,Z1, Z2)
n=10
Ui=runif(n,0,1)
Uj=runif(n,0,1)
Uij=matrix(runif(n^2,0,1),ncol=n)

g<-function(UiUjUij){
  #UiUjUij is a vector c(Ui,Uj,Uij)
  #function that generates the X, IVs Z and G sth cov (X,G)=0 to make 1SLS verification easier.
  Ui=UiUjUij[1]
  Uj=UiUjUij[2]
  Uij=UiUjUij[3]
  p<-2*pi
  a1=1
  a2=1
  b1=1
  b2=1
  e1=1
  e2=1
  N=rnorm(1,1,1)
  Z1=cos(a1*p*Ui)
  Z2=cos(a2*p*Uj)
  r=c(Z1,Z2)
  return (r)
} 
M=array(0,c(n,n,3))
M[,,1] = matrix(rep(Ui,n),c(n,n), byrow=FALSE) 
M[,,2] = matrix(rep(Uj,n),c(n,n), byrow=TRUE)
M[,,3] = Uij

A = apply(M,MARGIN = c(1,2),FUN = g)
print(dim(A))

A=array(0,c(n,n,7),dimnames=list(rep("",n),rep("",n),c('Y','X1','X2','cst','G','Z1','Z2')))
for(i in 1:n){
  for(j in 1:n){
    A[i,j,]=g(Ui[i],Uj[j],Uij[i,j])
  }
}
# A=dispatch_na(A)
# A=diag_na(A)
Yij  = A[,,1]
Xkij = A[,,2:3]
Gkij = A[,,4:5]
Zkij = A[,,6:7]

first_Betas2=Betas_1SLS(Gkij,Zkij,Xkij) #checked : Betas work
first_Betas2

D=D_SLS(Gkij,Zkij,Xkij, first_Betas = first_Betas)
mean((D[,,2]-Xkij[,,2])^2, na.rm=TRUE)

M=array(c(Gkij,Xkij),dim=c(n,n,4))
M2=array(c(Gkij,D),dim=c(n,n,4))

B2SLS = Beta_2SLS(D,Gkij,Yij)
B2SLS
Beta2(M,Yij)
Beta2(M2,Yij)





