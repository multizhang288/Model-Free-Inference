#This script will load the intermediate data results in Data Folder to estimate the FDR and Power
#Around 6-8 hours in local desktop.
#setwd("D:/Project_Revision/MF_Inference_Simulations/FDR")#set the path of this script
source("Functions_FDR.R")
root_path = getwd()
#####
#################
file_names <- dir(paste(root_path,"/Data",sep=""))
d = length(file_names)
for(case in 1:d)
{
  selects = file_names[case]
  split = as.vector(str_split(selects,"_"))[[1]]
  model = as.numeric(split[3])
  n = as.numeric(split[5])
  s = as.numeric(str_split(split[7],"\\.")[[1]])[1]
  data = read.csv(paste(root_path,"/Data/",selects,sep=""))
  ######################Size compute
  testvalue = data
  N = dim(data)[1]
  p = dim(data)[2]
  ###########################################################################
  bp = 2*log(p) +  2.75*log(log(p)) 
  ap = 2*log(p) +  3.5*log(log(p)) 
  ################################
  index = 1:p
  if(model == 4)
  {
    active = c(1:s)
    inactive = index[-active]
  }
  if(model == 5)
  {
    active = c(1:(s-2),(p-1):p)
    inactive = index[-active]
  }
  ###################################
  for(alpha in c(0.1,0.2))
  {
    record = matrix(rep(0,2*N),N,2)
    mrecord = c(-1,-1)
    cutoffs = rep(-1,N)
    ############################################
    for(i in 1:N)
    {
      FDR_cutoff = FDR.Est(testvalue[i,],bp,ap,active,inactive)
      cutoff = FDR_cutoff$cutoff
      record[i,1] = FDPtrue(testvalue[i,],cutoff,active,inactive)$FDP
      record[i,2] = FDPtrue(testvalue[i,],cutoff,active,inactive)$power
      cutoffs[i] = cutoff
   #   write.csv(cutoffs,paste("ZCutoff_Model_",model,"_n_",n,"_s_",s,"_",alpha*10,".csv",sep=""),row.names = FALSE,col.names = FALSE)
    }
    mrecord = apply(record,2,mean)
    print(c(model,n,s,alpha,mrecord))
   # write.csv(record,paste("ZEFDR_Model_",model,"_n_",n,"_s_",s,"_",alpha*10,".csv",sep=""),row.names = FALSE,col.names = FALSE)
    write.csv(mrecord,paste("ZTFDR_Model_",model,"_n_",n,"_s_",s,"_",alpha*10,".csv",sep=""),row.names = FALSE,col.names = FALSE)
  }
}
