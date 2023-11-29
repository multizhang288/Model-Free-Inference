library(splines)
library(MASS)
library(ncvreg)
##################
#setwd("D:/Project_Revision/MF_Inference_Simulations/Real Data") #Set the path of this script.
root_path = getwd()
source(paste(root_path,"/Functions_Real.R",sep=""))
#########################################
X <- t(read.table(paste(root_path,"/Data/RatEyeExpression.txt",sep="")))
data = matrix(rep(0,(dim(X)[1]-1)*dim(X)[2]), dim(X)[1]-1,dim(X)[2])
colnames(data) = X[1,]
for(i in 1:dim(X)[2])
{
  data[,i] = as.numeric(X[2:121,i])
}
X = data
# Standardise
# Remove predictor gene, which is tightly related to Ro131 as done in Segal et al
Y = as.numeric(X[,which(colnames(X)=="1389163_at")])
X = X[,-which(colnames(X)=="1389163_at")]
colmeans =apply(X,2,mean)
colsd = apply(X,2,sd)
select = names(sort(colsd,decreasing = TRUE))[1:3000]
X = scale(X)
X = X[,colnames(X) %in% select]
#Save the ID of genes
IDs = colnames(X)
#write.csv(IDs,"Gene_ID.csv",row.names = FALSE)
p = dim(X)[2]
n = dim(X)[1]
mu_y = mean(Y)
std_Y = sqrt(var(Y))
for(i in 1:n)
{
  Y[i] = (Y[i]-mu_y)/std_Y
}
#####
#X: Predictors
#Y: Response Variable
#h: Numbers of Transformation functions
#penalty: "lasso" or "SCAD"
#a: the scale coefficient of standard error
#####
penalty = "lasso"
h = 5 
a = 0.75
Result = Highdim.Real(X,Y,h,penalty,a)
record = Result$record #pvalues of all variables
#write.csv(record,paste("Real_Data_Pvalue_",100*a,"se",".csv",sep=""),row.names = FALSE)
data = record
################################
################################
################################
#FDR control
#data = read.csv(paste(root_path,"/Intermediate Data/Real_Data_Pvalue_75se.csv",sep=""))
#IDs = read.csv(paste(root_path,"/Intermediate Data/Gene_ID.csv",sep=""))
#IDs = t(IDs[,2])
#pvalue = data[,2]
pvalue = data #If do not load the saving data
testvalue = qchisq(1-pvalue,h) #statistic values
alpha = 0.05 #FDR level
Result_FDR = FDR.Real(pvalue,p,alpha)
cutoff = Result_FDR$cutoff
names = IDs[(testvalue>=cutoff)]
names[order(pvalue[testvalue>=cutoff])]
sort(pvalue[testvalue>=cutoff])
Result = rbind(c(a,names[order(pvalue[testvalue>=cutoff])]),c(a,sort(pvalue[testvalue>=cutoff])))
Result
#write.csv(Result,paste("Probes_Selected_",100*a,"se.csv",sep=""),row.names = FALSE)
