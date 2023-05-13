####This script generate the Figure 2
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
penalty = "lasso"
h = 5 
a = 0.75 #scale coefficient of standard error
index.ini = 1:p
#####
i = 2
index = index.ini[-i]
Z = X[,index]
Xi = X[,i]
cvfit = cv.ncvreg(Z,Xi,penalty=penalty,nlambda = 100)
lambda0 = cvfit$lambda
lambda1 = seq(max(lambda0),min(lambda0),length = 100)
cvfit = cv.ncvreg(Z,Xi,lambda = lambda1,penalty = penalty) 
cv1 = cvfit$cve[cvfit$min]+0.75*sd(cvfit$cve)
cv2 = cv1+1
kk = 1
while(cv2>=cv1)
{ 
  cv2 = cvfit$cve[kk]
  lambda = cvfit$lambda[kk]
  kk=kk+1
}
##################
plotdata = data.frame(lambda = lambda1, CV_error = cvfit$cve)
ggplot(plotdata,aes(x= lambda,y = CV_error)) +  theme_bw() + geom_point(size=1.5,color = "red") + 
  geom_hline(aes(yintercept = cvfit$cve[cvfit$min]),linetype = "dashed") +
  geom_hline(aes(yintercept = cv1),linetype = "dashed") +
  geom_vline(aes(xintercept = cvfit$lambda[cvfit$min]),linetype = "dotted",linewidth = 1.25) +
  geom_vline(aes(xintercept = cvfit$lambda[kk-1]),linetype = "dotted",linewidth = 1.25) +
  annotate("text", x=0.18, y=cv1*1.025, label="+0.75 se", angle=0,size = 7.5)+
  theme(text = element_text(size = 20),
        axis.text.x=element_text(size=rel(1.5), angle=0),
        axis.text.y=element_text(size=rel(1.5), angle=90))
