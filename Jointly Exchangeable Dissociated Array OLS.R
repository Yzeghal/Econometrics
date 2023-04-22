# Script for MCO in the case of jointly echangeable and dissociated arrays.
library(matrixcalc) #for matrix inversion and others
library(aod) #contains functions such that wald test we use in 2SLS. Homemade version is in Sandbox_tester
#source("C:/Users/tayoy/Documents/GitHub/Econometrics/Basic Regression.R", local = b <- new.env())
library(AER) #for a 2SLS built-in estimator (much better than homemade one in case of NAs)
library (R.utils) #to wrap array if we use built-in methods instead of homemade ones

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
B_hat<-function(X,Y){ #returnns B NORMALISED by nuber of observations used in the sum
  #----
  if (dim(X)[1]!=dim(X)[2]){
    stop("Matrix X is not square !")
  }
  if (dim(Y)[1]!=dim(Y)[2]){
    stop("Matrix Y is not square !")
  }
  if(dim(X)[1]!=dim(Y)[1]){
    stop("And Y should have the same size")
  }
  #----
  not_nax=!map_na(X)
  not_nay=!map_na(Y)
  not_na=not_nax&not_nay
  nb_obs=sum(not_na)
  M=matrix(0,dim(X)[3],dim(Y)[3])
  
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[1]){
      if (not_na[i,j]){
        M=M+matrix(X[i,j,])%*%t(matrix(Y[i,j,]))
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
  message(paste(sum(viable),'viable terms'))
  second_term=second_term/sum(viable)
  #print(second_term)
  matrix.inverse(J_hat)%*%second_term
} #Not to be used. Was expected to be better than Beta, but is not.

Beta<-function(X,Y, J_hat=NULL){ #Does the match with lm results (which simply drops all NA observations) and seems more precise...
  #add J_hat in argument if it was already calculated.
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
      message('J given in argument is used')
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
  message(paste(sum(viable),'viable terms'))
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


H_JEDA<-function(X,Y,b=NULL, ep=NULL){
  #give either b or diretly ep the residuals
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
  if (!is.null(b)){
    if (length(b)!=d){
      stop("b and dim(X)[3] should have the same length")
    }
    
  }
  else{
    if(is.null(ep)){
      stop("At least one of b or ep should not be NULL")
    }
    if (!is.matrix(ep)){
      stop("ep should be a matrix")
    }
    if(any(dim(ep)!=dim(X)[1:2])){
      stop("ep dimensions do not match X dimensions")
    }
  }
  
  #----
  if (!is.null(b)){
    ep=eps(X,Y,b)
  }
  
  not_naX=!map_na(X)
  not_naY=!map_na(Y)
  not_naEP=!map_na(ep)
  viable = (not_naX&not_naY)&not_naEP

  #tests that diagonal is NA
  paste(sum(viable),'viable ')
  
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
  if (any(sd==0)){stop("sd is 0. Sth to correct")}
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

OLS<-function(X,Y,hyp=0,model="JEDA",built_in_reg=TRUE){
  #safety tests
  #----
  if (!any(model==(c('JEDA','BASIC','NO ASVAR')))){
    stop("Model not recognised. Try \'JEDA\' or \'BASIC\' or \'NO ASVAR\'instead.")
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
    message=paste("\ndimension of X :", dimx, "\ndimension of Y :", dimy)
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
  J_hat=J(X,Y)
  if (det(J_hat)==0){
    print('J:')
    print (J_hat)
    stop ("J is singular !")
  }
  if (built_in_reg ==TRUE){
    x=wrap(X, map = list(NA,3))
    y=as.vector(Y)
    message("Beta Calculation with built-in method...")
    Beta_hat=matrix(lm(y~x -1)$coefficients)
  }
  else{
    message("Beta Calculation with homemade method...")
    Beta_hat = Beta(X,Y, J=J_hat)}
  
  if (model == "JEDA"){
    H_hat=H_JEDA(X,Y,b=Beta_hat)
    asvar=asymptotic_variance_OLS(J=J_hat, H=H_hat)
    std=standard_errors(asvar,dimx[1])
    ts=student_t(Beta_hat,std,hyp = hyp)
    p_val=p_values(Beta_hat, std ,dimx[1] , hyp = 0)
    if (sum(abs(hyp))==0){hyp = rep(0,length(Beta_hat))}
    values = matrix(c(Beta_hat,std,ts,p_val,hyp), ncol=5)
    colnames(values)=c('Beta_hat', 'Std_Err','Student_t','p-values', 'H0_hyp')
    F_test = wald.test(Sigma=asvar/dimx[1], b=Beta_hat, L=diag(length(Beta_hat)))$result$chi2[1] #normalized by nb of rows of X,
    Y_hat = estimate(X,Beta_hat)
    ep = eps(X,Y,Beta_hat)
    my=mean(Y, na.rm=TRUE)
  } 
  if (model == "BASIC"){
    H_hat=H_basic_OLS(X,Y,Beta_hat)
    asvar=asymptotic_variance_OLS(J=J_hat, H=H_hat)
    std=standard_errors(asvar,nb_obs)
    ts=student_t(Beta_hat,std,hyp = hyp)
    p_val=p_values(Beta_hat, std ,nb_obs, hyp = 0)
    if (sum(abs(hyp))==0){hyp = rep(0,length(Beta_hat))}
    values = matrix(c(Beta_hat,std,ts,p_val,hyp), ncol=5)
    colnames(values)=c('Beta_hat', 'Std_Err','Student_t','p-values', 'H0_hyp')
    F_test = wald.test(Sigma=asvar/nb_obs, b=Beta_hat, L=diag(length(Beta_hat)))$result$chi2[1] #normalized by nb of rows of X,
    Y_hat = estimate(X,Beta_hat)
    ep = eps(X,Y,Beta_hat)
    my=mean(Y, na.rm=TRUE)
  }
  if (model =="NO ASVAR"){
    values = matrix(Beta_hat)
    Y_hat = estimate(X,Beta_hat)
    my=mean(Y, na.rm=TRUE)
    asvar = NA
    F_test=NA
  }

  if (sd(Y,na.rm=TRUE)==0){
    R2=1
  }
  else{
    SCE = sum((Y_hat-my)^2, na.rm=TRUE)
    SCT = sum((Y-my)^2, na.rm = TRUE)
    R2<-SCE/SCT
  }

  return(list ('coefs'=values, 'var'=asvar, 'F'=F_test,"R2"= R2))
}

#2SLS functions
#----
FS_OLS<-function(G,Z,X,built_in_reg=TRUE){
  #G contains control variables
  #Z contains IV
  #X contains the variables that will be regressed on G and Z
  #----
  if(is.matrix(G)){
    G=array(G,dim=c(dim(G),1))
  }
  if(is.matrix(Z)){
    Z=array(Z,dim=c(dim(Z),1))
  }
  if(is.matrix(X)){
    X=array(X,dim=c(dim(X),1))
  }
  if (!(dim(X)[1:2]==dim(Z)[1:2]&&dim(Z)[1:2]==dim(G)[1:2])){
    stop('\nFirst 2 dimensions of X,Z and G should be the same')
  }
  
  if (!dim(X)[1]==dim(X)[2]){
    stop('Arrays should be square in their first 2 dimensions ')
  }
  if(any(replace(G[,,1],is.na(G[,,1]),1)!=1)){
    stop('Array G should contain the vector 1 in first position')
  }
  #----
  Betas = matrix(0,dim(G)[3]+dim(Z)[3],dim(X)[3])
  R2=rep(0,dim(X)[3])
  F_=rep(0,dim(X)[3])
  GZ=array(c(G,Z),dim=c(dim(G)[1:2], dim(G)[3]+dim(Z)[3]))
  for (i in 1:dim(X)[3]){
    reg=OLS(GZ,X[,,i],model = "JEDA", built_in_reg = built_in_reg) #We only need the R2 and usual F-tests
    Betas[,i]=reg$coefs[,1]
    R2[i]=reg$R2
    F_[i]=reg$F
  }
  return(list('Betas'=Betas,'R2'=R2,'F'=F_))
}

predIV<-function(G,Z,X,first_Betas=NULL){
  #----
  if(is.matrix(G)){
    G=array(G,dim=c(dim(G),1))
  }
  if(is.matrix(Z)){
    Z=array(Z,dim=c(dim(Z),1))
  }
  if(is.matrix(X)){
    X=array(X,dim=c(dim(X),1))
  }
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
  m=matrix(1,dim(G)[1],dim(G)[2])
  if(any(replace(G[,,1],is.na(G[,,1]),1)!=1)){
    stop('Array G should contain the vector 1 in first position')
  }
  #----
  if (is.null(first_Betas)){
    first_Betas = FS_OLS(G,Z,X)$Betas
  }
  GZ=array(c(G,Z),dim=c(dim(G)[1:2],dim(Z)[3]+dim(G)[3]))
  D=array(0,dim=dim(X))
  for (i in 1:dim(X)[3]){
    D[,,i]=estimate(GZ,first_Betas[,i])
  }
  return (D)
}

IV_LS<-function(G,Z,X,Y,hyp=0,built_in_reg=TRUE){
  #G : controls
  #Z : instruments
  #X : endogenous variable
  #Y : explained variable
  #----
  if(is.matrix(G)){
    G=array(G,dim=c(dim(G),1))
  }
  if(is.matrix(Z)){
    Z=array(Z,dim=c(dim(Z),1))
  }
  if(is.matrix(X)){
    X=array(X,dim=c(dim(X),1))
  }
  if (!(dim(X)[1:2]==dim(Z)[1:2]&&dim(Z)[1:2]==dim(G)[1:2]&&dim(Z)[1:2]==dim(Y))){
    stop('\nFirst 2 dimensions of G,Z,X,and Y should be the same')
  }
  
  if (!dim(X)[1]==dim(X)[2]){
    stop('Arrays should be square in their first 2 dimensions ')
  }
  if(any(replace(G[,,1],is.na(G[,,1]),1)!=1)){
    stop('Array G should contain the vector 1 in first position')
  }
  #----
  GX=array(c(G,X),dim=c(dim(G)[1:2],dim(X)[3]+dim(G)[3]))
  GZ=array(c(G,Z),dim=c(dim(G)[1:2],dim(Z)[3]+dim(G)[3])) 
  reg1=FS_OLS(G,Z,X,built_in_reg = FALSE)#regression of X on G and Z
  D=predIV(G,Z,X,first_Betas=reg1$Betas)
  GD=array(c(G,D),dim=c(dim(G)[1:2],dim(D)[3]+dim(G)[3]))
  
  if(built_in_reg==FALSE){
    Beta_2SLS=OLS(GD,Y, model = "NO ASVAR",built_in_reg=built_in_reg)$coefs[,1] #only Beta_hat matters here
  }
  else{
    y=as.vector(Y)
    g1=wrap.array(G,map = list(NA,3))[,-1]
    x1=wrap.array(X,map = list(NA,3))
    z1=wrap.array(Z,map = list(NA,3))
    Beta_2SLS = ivreg(y~g1+x1|g1+z1)$coefficients #puts coefs in the order used in the whole script.
  }

  ep=eps(GX,Y,b=Beta_2SLS)
  
  A=J(GZ,Y)
  B=B_hat(GZ,GX)
  H=H_JEDA(X=GZ,Y=Y,ep=ep)

  invA = matrix.inverse(A)
  bread = matrix.inverse(t(B)%*%invA%*%B)
  salad = t(B)%*%invA
  asvar = bread%*%salad%*%H%*%t(salad)%*%bread #Formula
  
  print("final regression")
  print(Beta_2SLS)
  dimx=dim(X)

  std=standard_errors(asvar,dimx[1])

  ts=student_t(Beta_2SLS,std,hyp = hyp)
  p_val=p_values(Beta_2SLS, std ,dimx[1] , hyp = 0)
  if (sum(abs(hyp))==0){hyp = rep(0,length(Beta_2SLS))}
  values = matrix(c(Beta_2SLS,std,ts,p_val,hyp), ncol=5)
  colnames(values)=c('Beta_hat', 'Std_Err','Student_t','p-values', 'H0_hyp')
  F_test = wald.test(Sigma=asvar/dimx[1], b=Beta_2SLS, L=diag(length(Beta_2SLS)))$result$chi2[1] #normalized by nb of rows of X,
  Y_hat = estimate(GD,Beta_2SLS)
  ep2 = eps(GD,Y,Beta_2SLS)
  
  my=mean(Y, na.rm=TRUE)
  SCE = sum((Y_hat-my)^2, na.rm=TRUE)
  SCT = sum((Y-my)^2, na.rm = TRUE)
  R2<-SCE/SCT
  
  return(list('SLS'=list('coefs'=values, 'var'=asvar, 'F'=F_test,"R2"= R2),
              'FLS'=list('R2'=reg1$R2,'F'=reg1$F)))
}