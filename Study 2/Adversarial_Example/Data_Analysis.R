library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)
############################################
#set the path of the intermediate results without MGT
setwd("D:/Project_Revision/Revision_1116/JASARevision/Model_Free_Inference_1118/Study 2/Adversarial_Example/Intermediate Data Nomgt")
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
  selects = file_names[i]#BS 200 400
  data = do.call(rbind,lapply(selects,read.csv)) 
  data = t(data)
  df = as.numeric(data[6])
  ERR = as.numeric(data[c(11,13:16)])
  Size = c(df,ERR)
  Record = rbind(Record,Size)
}
Record = data.frame(Record)
colnames(Record) = c("df","X5","X1","X2","X3","X4")
#######################
Power = Record[order(as.numeric(Record$df)),]
df = c(0.01,0.05,0.1,seq(0.2,2,0.2))
data1 = data.frame(df = df, Power = as.numeric(Power$X5), index = rep("X5",13))
data2 = data.frame(df = df, Power = as.numeric(Power$X1), index = rep("X1",13))
data3 = data.frame(df = df, Power = as.numeric(Power$X2), index = rep("X2",13))
data4 = data.frame(df = df, Power = as.numeric(Power$X3), index = rep("X3",13))
data5 = data.frame(df = df, Power = as.numeric(Power$X4), index = rep("X4",13))
plotdata = rbind(data1,data2,data3,data4,data5)
Figure1 =  ggplot(plotdata,aes(x = df ,Power,group=index, color = index,shape = index,linetype = index)) + scale_x_continuous(breaks=c(0.01,0.2,0.6,1,1.4,2)) + ggtitle("No MGT")+ theme_bw()+  theme(text = element_text(size = 16),plot.title = element_text(hjust = 0.5)) + geom_point(size = 1.5)+ geom_line(linewidth = 0.5)
Figure1 = Figure1 + labs(y="Empirical Rejection Rate") +  geom_hline(yintercept=0.05, linetype="dashed", color = "red") + geom_hline(yintercept=0.064, linetype="dashed", color = "blue") + geom_hline(yintercept=0.036, linetype="dashed", color = "blue") 
Figure1 = Figure1 + scale_y_continuous(breaks = c(0,0.05,0.25,0.5,0.75,1))
#################################
#set the path of the intermediate results with MGT
setwd("D:/Project_Revision/Revision_1116/JASARevision/Model_Free_Inference_1118/Study 2/Adversarial_Example/Intermediate Data MGT")
file_names <- dir()
d = length(file_names)
Record  = NULL
for(i in 1:d)
{
  selects = file_names[i]#BS 200 400
  data = do.call(rbind,lapply(selects,read.csv)) 
  data = t(data)
  df = as.numeric(data[6])
  ERR = as.numeric(data[c(11,13:16)])
  Size = c(df,ERR)
  Record = rbind(Record,Size)
}
Record = data.frame(Record)
colnames(Record) = c("df","X5","X1","X2","X3","X4")
#######################
Power = Record[order(as.numeric(Record$df)),]
df = c(0.01,0.05,0.1,seq(0.2,2,0.2))
data1 = data.frame(df = df, Power = as.numeric(Power$X5), index = rep("X5",13))
data2 = data.frame(df = df, Power = as.numeric(Power$X1), index = rep("X1",13))
data3 = data.frame(df = df, Power = as.numeric(Power$X2), index = rep("X2",13))
data4 = data.frame(df = df, Power = as.numeric(Power$X3), index = rep("X3",13))
data5 = data.frame(df = df, Power = as.numeric(Power$X4), index = rep("X4",13))
plotdata = rbind(data1,data2,data3,data4,data5)
Figure2 =  ggplot(plotdata,aes(x = df ,Power,group=index, color = index,shape = index,linetype = index)) +scale_x_continuous(breaks=c(0.01,0.2,0.6,1,1.4,2)) + ggtitle("MGT")+ theme_bw()+  theme(text = element_text(size = 16),plot.title = element_text(hjust = 0.5)) + geom_point(size = 1.5)+ geom_line(linewidth = 0.5)
Figure2 = Figure2 + labs(y="Empirical Rejection Rate") + labs(y="Empirical Rejection Rate") +  geom_hline(yintercept=0.05, linetype="dashed", color = "red") + geom_hline(yintercept=0.064, linetype="dashed", color = "blue") + geom_hline(yintercept=0.036, linetype="dashed", color = "blue")
Figure2 = Figure2 + scale_y_continuous(breaks = c(0,0.05,0.25,0.5,0.75,1))

grid.arrange(Figure1, Figure2, nrow = 1)


