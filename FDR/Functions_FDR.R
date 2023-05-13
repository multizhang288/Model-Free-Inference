#Load Packages
#install.packages(c("pracma","dplyr","splines","stringr","MASS","ncvreg"))
library(pracma)
library(dplyr)
library(splines)
library(stringr)
library(MASS)
library(ncvreg)
library(mvtnorm)
################
gen_error<-function(N,p,rho,seed){
  X = matrix(NA,N,p)
  X[,1] = rnorm(N)
  for(ii in 2:p){
    X[,ii] = rho*X[,(ii-1)] + sqrt(1-rho^2)*rnorm(N)
  }
  return(X)
}

data_seq = function(x,h)
{
  value = seq(0,1,by=1/(h))
  then = quantile(x,probs=value)
  return(as.numeric(then))
}
#########################################
FDP = function(value,cutoff,active,inactive)
{
  I0 = sum(ifelse(value[inactive]>=cutoff,1,0))
  I1 = sum(ifelse(value[active]>=cutoff,1,0))
  I = I0 + I1
  quantile = 1-pchisq(cutoff,5)
  FDP = p*quantile/max(I,1)
  power = I1/length(active)
  list(FDP = FDP, power = power)
}

FDPtrue = function(value,cutoff,active,inactive)
{
  I0 = sum(ifelse(value[inactive]>=cutoff,1,0))
  I1 = sum(ifelse(value[active]>=cutoff,1,0))
  I = I0 + I1
  FDP = I0/max(I,1)
  power = I1/length(active)
  list(FDP = FDP, power = power)
}
################
#n: sample size
#p: dimension of covariates
#h: degrees of freedom for B-spline Basis
#s: sparsity level
#rho: Sigma~AR(1) with rho
#penalty: Penalty function, you can choose one between lasso and SCAD
#model: data-generation model
################
Highdim = function(n,p,s,rho,h,penalty,model,q)
{
  index.ini = 1:p
  if(model==4)
  {
    set.seed(q)
    X = gen_error(n,p,rho,q)
    beta.true = c(rep(1,s),rep(0,p-s))
    error <- rnorm(n = n,0,sigma) 
    Y = X%*%beta.true + 1*error
    Y = as.vector(Y)
  }
  if(model == 5)
  {
    set.seed(q)
    X = gen_error(n,p,rho,q)
    beta.true = c(rep(1,s-2),rep(0,p-s+2))
    error <- rnorm(n = n,0,sigma) 
    Y = (X%*%beta.true)/(0.5+(1.5+X[,p-1]+X[,p])^2)+0.1*error  
  }
  ################################
  record_pvalue = record_pvalue1 = record_stat = rep(-1,p)
  if(h>=2)
  {
    temp = bs(Y,df=h,degree = 2)
  }
  if(h==1)
  {
    percens = data_seq(Y,h)
    temp = NULL
    for(k in 1:h)
    {
      fyk = ifelse(Y>=percens[k]&Y<=percens[k+1],Y,0)
      temp = cbind(temp,fyk)
    }
  }
  #####################Decorrelated score part
  model2 = cv.ncvreg(X,Y,penalty = penalty,seed = q,nlambda = 50)
  fit = model2$fit
  beta = fit$beta[,model2$min]
  hatdecor = as.numeric(beta)
  ###############################Reduce computational cost
  model.save = rep(list(1),h)
  FYK =rep(list(1),h)
  for(k in 1:h)
  {
    fyk = temp[,k]
    model1 = cv.ncvreg(X,fyk,penalty = penalty)
    fit = model1$fit
    beta = fit$beta[,model1$min]
    FYK[[k]] = fyk
    model.save[[k]] = as.numeric(beta)#[-1]
  }
  ##########################
  for(i in 1:p)
  {
    index = index.ini[-i]
    Z = X[,index]
    Xi = X[,i]
    cvfit = cv.ncvreg(Z,Xi,penalty=penalty)
    fit = cvfit$fit
    beta <- fit$beta[,cvfit$min]
    beta = as.numeric(beta)
    Z1 = cbind(1,Z)
    res = Xi - Z1%*%beta
    fisher = mean(res^2)
    eta = eta1  = matrix(rep(0,n*h),n,h) 
    for(k in 1:h)
    {
      fyk = FYK[[k]]
      model2 = model.save[[k]]
      beta = as.numeric(model2)
      inpt = beta[1]
      beta = beta[-1]
      hattheta = beta[-i]
      hatgammak = hattheta#as.numeric(model1)
      eta.temp = fyk-Z%*%hatgammak - inpt
      eta.temp1 = fyk -Z%*%hatgammak - inpt
      eta[,k] = eta.temp*res
    }
    omega = t(eta)%*%eta/n
    eta.mean = apply(eta,2,mean)
    Wni = n*t(eta.mean)%*%ginv(omega)%*%eta.mean #sometimes omega is not invertible, we use ginv here.
    record_stat[i] = Wni
    record_pvalue[i] = 1-pchisq(Wni,h)
    time1 = Sys.time()
    ###########################Decorrelated score
    resdecor = Y-Z1%*%hatdecor[-c(i+1)]
    sigmahat = t(resdecor)%*%resdecor/n
    etadecor = resdecor*res
    Un = sum(etadecor)/sqrt(sigmahat*n*fisher)
    record_pvalue1[i] = 2*(1-pnorm(abs(Un)))
    #############################
  }
  list(pvalue = record_pvalue, pvalue1 = record_pvalue1, stat = record_stat)
}
################################FDR estimation
FDR.Est = function(testvalue,bp,ap,active,inactive)
{
  first = 0
  end = bp
  u0 = FDP(testvalue,first,active,inactive)$FDP
  u1 = FDP(testvalue,end,active,inactive)$FDP
  criterion = abs(u1-u0)
  step = 0
  while(criterion>=10^(-5)&step<=30)
  {
    temp0 = first
    temp1 = end
    middle = 0.5*first + 0.5*end
    tvalue = FDP(testvalue,middle,active,inactive)$FDP
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
  if(abs(cutoff-bp)<=0.0001)
  {
    cutoff = ap
  }
  list(cutoff = cutoff)
}