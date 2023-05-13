library(tidyr)
library(stringr)
library(ggplot2)
############################################
setwd("D:/Project_Revision/MF_Inference_Simulations/Study 2/Semi-synthetic/Intermediate Data")#set the path of the Intermediate Data
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
selects = file_names[i]#BS 200 400
split =  as.vector(str_split(selects,"_"))[[1]]
model = as.numeric(split[3])
data = do.call(rbind,lapply(selects,read.csv)) 
data = data[9:36,]
data = as.numeric(t(data))
Wn = c(model,data[4:14])
TNL = c(model,data[18:28])
Size = rbind(Wn,TNL)
Record = rbind(Record,Size)
}
write.csv(Record,"D:/Project_Revision/MF_Inference_Simulations/Study 2/Semi-synthetic/Result/Semi-synthetic_Result.csv")

