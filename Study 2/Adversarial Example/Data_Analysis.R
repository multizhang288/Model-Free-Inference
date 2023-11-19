library(tidyr)
library(stringr)
library(ggplot2)
############################################
#set the path of the intermediate results
setwd("D:/Project_Revision/Revision_1116/JASARevision/Model_Free_Inference_1118/Study 2/Adversarial Example/Intermediate Data")
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
  selects = file_names[i]#BS 200 400
  split =  as.vector(str_split(selects,"_"))[[1]]
  model = as.numeric(split[3])
  data = do.call(rbind,lapply(selects,read.csv)) 
  data = t(data)
  df = as.numeric(data[6])
  ERR = as.numeric(data[c(11,13)])
  Size = c(df,ERR)
  Record = rbind(Record,Size)
}
Record = data.frame(Record)
colnames(Record) = c("df","X2","X1")
#######################
Power = Record[order(as.numeric(Record$df)),]
df = 1:11
data1 = data.frame(df = df, Power = as.numeric(Power$X2), index = rep("X2",11))
data2 = data.frame(df = df, Power = as.numeric(Power$X1), index = rep("X1",11))
plotdata = rbind(data1,data2)
Figure =  ggplot(plotdata,aes(x = df ,Power,group=index, color = index,shape = index,linetype = index)) + ggtitle("Mixture distribution with respect to df")+ theme_bw()+  theme(plot.title = element_text(hjust = 0.5)) + geom_point(size = 1.5)+ geom_line(linewidth = 0.5)
Figure
#################################



