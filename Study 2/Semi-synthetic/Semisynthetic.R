#set the path of this script
root_path = getwd()#"D:/Project_Revision/MF_Inference_Simulations/Study 2/Semi-synthetic"
source(paste(root_path,"/Functions_Semi.R",sep=""))
########
X <- t(read.table(paste(root_path,"/RatEyeExpression.txt",sep="")))
data = matrix(rep(0,(dim(X)[1]-1)*dim(X)[2]), dim(X)[1]-1,dim(X)[2])
colnames(data) = X[1,]
for(i in 1:dim(X)[2])
{
  data[,i] = as.numeric(X[2:121,i])
}
X = data
Y = as.numeric(X[,which(colnames(X)=="1389163_at")])
X = X[,-which(colnames(X)=="1389163_at")]
colmeans = apply(X,2,mean)
colsd  = apply(X,2,sd)
select = names(sort(colsd,decreasing = TRUE))[1:20]
X = scale(X)
X0 = X[,colnames(X) %in% select]
p = dim(X)[2]
n = dim(X)[1]
mu_y = mean(Y)
std_Y = sqrt(var(Y))
for(i in 1:n)
{
  Y[i] = (Y[i]-mu_y)/std_Y
}
###########################Run 500 times
s = 10
set.seed(1)
location = 1:s #index of active variables
h = 7 #
hfix = -1#hfix!=-1, set your own h. Otherwise h is fixed as 7
delta = 1#signal strength
for(model in c(1:5))
{
Record = Record0 = Record1 = Record_decor = NULL#
for(k in 1:1000)
{
  time0 = Sys.time()
  Result0 = Highdim.Semi(X0,h,hfix,s,location,"Bspline","lasso",model,delta,seed = k)#X,h,hfix,s,d,location,choice2,active,penalty,model,delta,seed
  Score_save = Result0$pvalue
  Sigma_save = Result0$stat
  Decor_save = Result0$pvalue.decor
  Record_decor_add = ifelse(Result0$pvalue.decor<=0.05,1,0)
  Record_add = ifelse(Result0$pvalue<=0.05,1,0)
  Record = rbind(Record,Record_add)
  Record0 = rbind(Record0, Score_save)
  Record1 = rbind(Record1, Sigma_save)
  Record_decor = rbind(Record_decor, Record_decor_add)
  time1 = Sys.time()
  Type1 = apply(Record,2,mean)#Empirical rejection rate of Wn
  Type2 = apply(Record_decor,2,mean)#Empirical rejection rate of TNL
  print(c(time1-time0,model,k))
  print(Type1)
  h1 = Result0$h
  d1 = Result0$d
  Result = c(k,hfix,model,delta,"Bspline",s,d1,h1,Type1,Type2)
  write.csv(Result,paste("New-Semi-Impact_Model_",model,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE)
}
}
  