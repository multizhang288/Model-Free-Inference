library(tidyr)
library(stringr)
library(ggplot2)
############################################
#Set your own file path firstly
setwd("D:/Project_Revision/MF_Inference_Simulations/Study 1/Results_h_5/Intermediate Data")
file_names <- dir()
d = length(file_names)
for(penalty in c("lasso","SCAD"))
{
Record = Record1 = NULL
for(i in 1:12)
{
      selects = file_names[i]#BS 200 400
      split =  as.vector(str_split(selects,"_"))[[1]]
      if(split[5]=="lasso")
      {
      data = do.call(rbind,lapply(selects,read.csv)) 
      data = t(data)
      Wn = c("Wn",data[c(2,4,5,6:15)])
      TNL = c("TNL",data[c(2,4,5,16:25)])
      model = as.numeric(Wn[3])
      n = as.numeric(Wn[2])
      penalty = Wn[4]
      if(model==1&n==200&penalty=="lasso")
      {
        ZZ = do.call(rbind,lapply(file_names[13],read.csv))[,-1]
        TZZ = c("TZZ",data[c(2,4,5)],apply(ZZ,2,mean)[c(1:5,1996:2000)])
        Size = rbind(Wn,TNL,TZZ)
        Record = rbind(Record,Size)
      }
      if(model==1&n==400&penalty=="lasso")
      {
        ZZ = do.call(rbind,lapply(file_names[14],read.csv))[,-1]
        TZZ = c("TZZ",data[c(2,4,5)],apply(ZZ,2,mean)[c(1:5,1996:2000)])
        Size = rbind(Wn,TNL,TZZ)
        Record = rbind(Record,Size)
      }
      if(model==2&n==200&penalty=="lasso")
      {
        ZZ = do.call(rbind,lapply(file_names[15],read.csv))[,-1]
        TZZ = c("TZZ",data[c(2,4,5)],apply(ZZ,2,mean)[c(1:5,1996:2000)])
        Size = rbind(Wn,TNL,TZZ)
        Record = rbind(Record,Size)
      }
      if(model==2&n==400&penalty=="lasso")
      {
        ZZ = do.call(rbind,lapply(file_names[16],read.csv))[,-1]
        TZZ = c("TZZ",data[c(2,4,5)],apply(ZZ,2,mean)[c(1:5,1996:2000)])
        Size = rbind(Wn,TNL,TZZ)
        Record = rbind(Record,Size)
      }
      if(model==3&n==200&penalty=="lasso")
      {
        ZZ = do.call(rbind,lapply(file_names[17],read.csv))[,-1]
        TZZ = c("TZZ",data[c(2,4,5)],apply(ZZ,2,mean)[c(1:5,1996:2000)])
        Size = rbind(Wn,TNL,TZZ)
        Record = rbind(Record,Size)
      }
      if(model==3&n==400&penalty=="lasso")
      {
        ZZ = do.call(rbind,lapply(file_names[18],read.csv))[,-1]
        TZZ = c("TZZ",data[c(2,4,5)],apply(ZZ,2,mean)[c(1:5,1996:2000)])
        Size = rbind(Wn,TNL,TZZ)
        Record = rbind(Record,Size)
      }
      }
      if(split[5]=="SCAD")
      {
        data = do.call(rbind,lapply(selects,read.csv)) 
        data = t(data)
        Wn = c("Wn",data[c(2,4,5,6:15)])
        TNL = c("TNL",data[c(2,4,5,16:25)])
        Size = rbind(Wn,TNL)
        Record1 = rbind(Record1,Size)
      }
}
Record = data.frame(Record)
Record[,5:14] = sapply(Record[,5:14],as.numeric)
Record[,5:14] = sapply(Record[,5:14],function(x) round(x,3))
colnames(Record)[1:4] = c("Method","n","model","penalty")
Record1 = data.frame(Record1)
Record1[,5:14] = sapply(Record1[,5:14],as.numeric)
Record1[,5:14] = sapply(Record1[,5:14],function(x) round(x,3))
colnames(Record1)[1:4] = c("Method","n","model","penalty")
#write.csv(Record,paste("D:/Project_Revision/MF_Inference_Simulations/Study 1/Results_h_5/Result/Lasso","_h_",5,".csv",sep=""))
#write.csv(Record1,paste("D:/Project_Revision/MF_Inference_Simulations/Study 1/Results_h_5/Result/SCAD","_h_",5,".csv",sep=""))
}
Record
Record1

