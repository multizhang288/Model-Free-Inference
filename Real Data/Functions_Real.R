#Load Packages
library(pracma)
library(dplyr)
library(splines)
library(stringr)
library(MASS)
library(ncvreg)
################
FDP = function(value,cutoff)
{
  I = sum(ifelse(value>=cutoff,1,0))
  quantile = 1-pchisq(cutoff,5)
  FDP = p*quantile/max(I,1)
  list(FDP = FDP)
}
#################################
Highdim.Real = function(X,Y,h,penalty,a)
  {
  ################
  #X: Predictors
  #Y: Response Variable
  #h: Numbers of Transformation functions
  #penalty: "lasso" or "SCAD"
  #a: the scale coefficient of standard error
  ################
  p = dim(X)[2]
  n = dim(X)[1]
  index.ini = 1:p
  set.seed(0)
  temp = bs(Y,df=h,degree = 2)
  model.save = rep(list(1),h)
  FYK =rep(list(1),h)
  for(k in 1:h)
  {
    fyk = temp[,k]
    X1 = cbind(1,X)
    factor = c(0,rep(1,p))
    cvfit = cv.ncvreg(X,fyk,penalty = penalty,nlambda = 100)
    lambda0 = cvfit$lambda
    lambda1 = seq(max(lambda0),min(lambda0),length = 100)
    cvfit = cv.ncvreg(X,fyk,lambda = lambda1,penalty = penalty)  
    fit = cvfit$fit
    cv1 = cvfit$cve[cvfit$min]+ a*sd(cvfit$cve)
    cv2 = cv1+1
    kk = 1
    while(cv2>=cv1)
    { 
      cv2 = cvfit$cve[kk]
      lambda = cvfit$lambda[kk]
      kk=kk+1
    }
    fit = cvfit$fit
    beta = fit$beta[,kk-1]
    FYK[[k]] = fyk
    model.save[[k]] = as.numeric(beta)
  }
  record = rep(-1,p)
  for(i in 1:p)
  {
    time0 = Sys.time()
    index = index.ini[-i]
    Z = X[,index]
    Xi = X[,i]
    intp1 = rep(1,n)
    Z1 = cbind(1,Z)
    factor = c(0,rep(1,p-1))
    cvfit = cv.ncvreg(Z,Xi,penalty = penalty,seed = i,nlambda = 100)
    lambda0 = cvfit$lambda
    lambda1 = seq(max(lambda0),min(lambda0),length = 100)
    cvfit = cv.ncvreg(Z,Xi,lambda = lambda1,seed = i,penalty = penalty)  
    fit = cvfit$fit
    cv1 = cvfit$cve[cvfit$min]+ a*sd(cvfit$cve)
    cv2 = cv1+1
    kk = 1
    while(cv2>=cv1)
    { 
      cv2 = cvfit$cve[kk]
      lambda = cvfit$lambda[kk]
      kk=kk+1
    }
    hattheta = fit$beta[,kk-1]
    res = Xi - Z1%*%hattheta
    eta = matrix(rep(0,n*h),n,h)
    for(k in 1:h)
    {
      fyk = temp[,k]
      beta = model.save[[k]]
      hatgammak = as.numeric(beta)[-c(i+1)]
      eta.temp = fyk-Z1%*%hatgammak
      eta[,k] = eta.temp*res
    }
    omega = t(eta)%*%eta/n
    eta.mean = apply(eta,2,mean)
    Wni = n*t(eta.mean)%*%ginv(omega)%*%eta.mean #sometimes omega is not invertible, we use ginv here.
    record[i] =  1-pchisq(Wni,h)#
    print(c(i,record[i]))
  }
  list(record = record)
}
##############################
FDR.Real = function(pvalue,p,alpha)
{
  bp = 2*log(p) + 2.75*log(log(p))
  ap = 2*log(p) + 3.5*log(log(p))
  testvalue = qchisq(1-pvalue,h)
  first = 0
  end = bp
  u0 = FDP(testvalue,first)$FDP
  u1 = FDP(testvalue,end)$FDP
  criterion = abs(u1-u0)
  step = 0
  while(criterion>=10^(-5)&step<=2000)
  {
    temp0 = first
    temp1 = end
    middle = 0.5*first + 0.5*end
    tvalue = FDP(testvalue,middle)$FDP
    if(tvalue<alpha)
    {
      end = middle
    }
    if(tvalue>alpha)
    {
      first = middle
    }
    if(tvalue==alpha)
    {
      break
    }
    criterion = abs(tvalue-alpha)
    cutoff = middle
    step = step+1
  }
  if(abs(cutoff-bp)<=0.001)
  {
    cutoff = ap
  } 
  list(cutoff = cutoff)
}
