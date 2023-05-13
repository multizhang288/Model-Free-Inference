library(tidyr)
library(stringr)
library(ggplot2)
############################################
setwd("D:/Project_Revision/MF_Inference_Simulations/Study 2/Transformation/Intermediate Data")#Set the path
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
selects = file_names[i]#BS 200 400
split =  as.vector(str_split(selects,"_"))[[1]]
model = as.numeric(split[3])
s = as.numeric(split[9])
data = do.call(rbind,lapply(selects,read.csv)) 
data = t(data)
Trans = data[7]
Inact_Wn = as.numeric(data[16])
d = ifelse(model==1,1,0) + ifelse(model==2,2,0) + ifelse(model==3,3,0) + ifelse(model==4,4,0) + ifelse(model==5,5,0)
if(model==1)
{
  Index_Wn = c(17)
  Index_TNL = c(17+s+4)
  Size_Wn = rep("$***$",5)
   Size_Wn[1:d] = as.numeric(data[Index_Wn])
}

if(model == 2)
{
  intv = floor(0.8*s)
  Index_Wn = c(17, 17+intv)
  Size_Wn = rep("$***$",5)
  Size_Wn[1:d] = as.numeric(data[Index_Wn])
}

if(model>2)
{
  intv = floor(s/d)
  location = c(0:(d-1))*(intv) + 1
  Index_Wn = 16 + location
  Size_Wn = rep("$***$",5)
  Size_Wn[1:d] = as.numeric(data[Index_Wn])
}
Size = c(model,Trans,Inact_Wn,Size_Wn)
Record = rbind(Record,Size)
}
Record = data.frame(Record)
write.csv(Record,"D:/Project_Revision/MF_Inference_Simulations/Study 2/Transformation/Results/Impact_Transformation.csv")#Check the path


