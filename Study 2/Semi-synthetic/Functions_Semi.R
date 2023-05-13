######################################
gen_error<-function(N,p,rho){
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
######################################Synthetic data by adding noise term
Highdim.Semi = function(X0,h,hfix,s,location,Trans,penalty,model,delta,seed)
{
  n = dim(X0)[1]
  p = 1000 - dim(X0)[2]
  set.seed(seed)
  error = rnorm(n,0,1)
  Xa =gen_error(n,p,0.5)
  X = cbind(X0,Xa)
  p = dim(X)[2]
  index.ini = 1:p
  interest = c(32,64,128,256,location)
  if(model==-1)
  {
    Y = mean_fun + error
  }
  if(model==1)
  {
    set.seed(seed)
    error <- rnorm(n = n,0,1) 
    beta.true = rep(0,p)
    beta.true[location] = delta
    Y =(X%*%beta.true)+ 1*error
    Y = as.vector(Y)
    
    d=1
    if(hfix==-1)
    {
      h=7
    }
  }
  if(model==2)
  {
    set.seed(seed)
    error <- rnorm(n = n,0,1) 
    beta1 = beta2 = rep(0,p)
    ll = length(location)
    beta1[location[1:(floor(ll*0.8))]] = delta
    beta2[location[-c(1:(floor(ll*0.8)))]] = delta
    Y = (X%*%beta1)/(0.5+(1.5+X%*%beta2)^2)+0.1*error  
    
    d=2
    if(hfix==-1)
    {
      h=7
    }
  }
  if(model == 3)
  {
    error <- rnorm(n = n,0,1) 
    beta1 = beta2 = beta3 = rep(0,p)
    ll = length(location)
    beta1[location[1:floor(ll/3)]] = delta
    beta2[location[c((floor(ll/3)+1):(2*floor(ll/3)))]] = delta
    beta3[location[-c(1:(2*floor(ll/3)))]] = delta
    Y = sinh(X%*%beta1) + sinh(X%*%beta2) + sinh(X%*%beta3) +0.5*error
    
    d=3
    if(hfix==-1)
    {
      h=7
    }
  }
  if(model == 4)
  {
    error <- rnorm(n = n,0,1) 
    beta1 = beta2 = beta3 = beta4 = rep(0,p)
    ll = length(location)
    beta1[location[1:floor(ll/4)]] = delta
    beta2[location[c((floor(ll/4)+1):(2*floor(ll/4)))]] = delta
    beta3[location[c((2*floor(ll/4)+1):(3*floor(ll/4)))]] = delta
    beta4[location[-c(1:(3*floor(ll/4)))]] = delta
    Y = X%*%beta1 + sinh(X%*%beta2) + 0.3*exp(X%*%beta3) + abs(X%*%beta4+3) +1*error
    
    d=4
    if(hfix==-1)
    {
      h=d+3
    }
  }
  if(model == 5)
  {
    error <- rnorm(n = n,0,1) 
    beta1 = beta2 = beta3 = beta4 = beta5 = rep(0,p)
    ll = length(location)
    beta1[location[1:floor(ll/5)]] = delta
    beta2[location[c((floor(ll/5)+1):(2*floor(ll/5)))]] = delta
    beta3[location[c((2*floor(ll/5)+1):(3*floor(ll/5)))]] = delta
    beta4[location[c((3*floor(ll/5)+1):(4*floor(ll/5)))]] = delta
    beta5[location[-c(1:(4*floor(ll/5)))]] = delta
    Y = X%*%beta1 + sinh(X%*%beta2) + exp(X%*%beta3) + abs(X%*%beta4+3) + 4*sin(X%*%beta5) + 1*error
    
    d=5
    if(hfix==-1)
    {
      h=d+2
    }
  }
  ##############################Create the basis function with h = max(d+3,5)
  ##############################Create the basis function with h = max(d+3,5)
  if(Trans == "Bspline")
  {
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
  }
  if(Trans == "SIR")
  {
    percens = data_seq(Y,h)
    temp = NULL
    for(k in 1:(h-1))
    {
      fyk = ifelse(Y>=percens[k]&Y<percens[k+1],1,0)
      temp = cbind(temp,fyk)
    }
    fyk = ifelse(Y>=percens[h]&Y<=percens[h+1],1,0)
    temp = cbind(temp,fyk)
    
  }
  if(Trans == "SIR2")
  {
    percens = data_seq(Y,h)
    temp = NULL
    for(k in 1:(h-1))
    {
      fyk = ifelse(Y>=percens[k]&Y<percens[k+1],Y,0)
      temp = cbind(temp,fyk)
    }
    fyk = ifelse(Y>=percens[h]&Y<=percens[h+1],Y,0)
    temp = cbind(temp,fyk)
  }
  if(Trans == "Poly")
  {
    temp = NULL
    T = max(abs(Y))
    for(k in 1:h)
    {
      fyk1 = (Y/T)^k
      temp = cbind(temp,fyk1)
    }
  }
  ###############################
  model.save = rep(list(1),h)
  FYK =rep(list(1),h)
  for(k in 1:h)
  {
    fyk = temp[,k]
    model1 = cv.ncvreg(X,fyk,penalty = penalty,nfolds = 5)
    cvfit = model1$fit
    beta = cvfit$beta[,model1$min] 
    FYK[[k]] = fyk
    model.save[[k]] = as.numeric(beta)
  }
  ##############################Decorrelated score
  model2 = cv.ncvreg(X,Y,nfolds=5,penalty=penalty)
  cvfit  = model2$fit
  beta = cvfit$beta[,model2$min]
  hatdecor = as.numeric(beta)
  ##############################
  record.stat = record.pvalue = record.pvalue1= rep(-1,p)
  ###############################################
  for(i in c(interest))
  {
    time0 = Sys.time()
    index = index.ini[-i]
    Z = X[,index]
    Xi = X[,i]
    cvfit = cv.ncvreg(Z,Xi,penalty=penalty,nfolds = 5)
    fit = cvfit$fit
    beta <- fit$beta[,cvfit$min]
    beta = as.numeric(beta)
    hattheta = as.numeric(beta)[-1]
    Z1 = cbind(1,Z)
    res = Xi - Z1%*%beta
    fisher = mean(res^2)
    eta = eta1  = matrix(rep(0,n*h),n,h) 
    for(k in 1:h)
    {
      fyk = temp[,k]
      fit = model.save[[k]]
      beta = fit
      hatgammak = as.numeric(beta)[-c(i+1)]
      eta.temp = fyk-Z1%*%hatgammak
      eta.temp1 = fyk -Z1%*%hatgammak
      eta[,k] = eta.temp
      eta1[,k] = eta.temp1*res
    }
    omega = t(eta1)%*%eta1/n
    omega0 = t(eta)%*%eta/n
    
    eta.mean = apply(eta1,2,mean)
    Wni = n*t(eta.mean)%*%ginv(omega)%*%eta.mean#/as.numeric(fisher) #sometimes omega is not invertible, we use ginv here.
    record.stat[i] = Wni
    record.pvalue[i] = 1-pchisq(Wni,h)#ifelse(Wni>=qchisq(0.95,h),1,0)
    time1 = Sys.time()
    ################
    resdecor = Y-Z1%*%hatdecor[-c(i+1)]
    sigmahat = t(resdecor)%*%resdecor/n
    etadecor = resdecor*res
    Un = sum(etadecor)/sqrt(sigmahat*n*fisher)
    record.pvalue1[i] = 2*(1-pnorm(abs(Un)))
  }
  list(pvalue = record.pvalue[interest], pvalue.decor = record.pvalue1[interest],stat = record.stat[interest],h = h,d =d)
}
