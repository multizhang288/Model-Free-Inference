library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(gridExtra)
########################
###############
#setwd("D:/Project_Revision/MF_Inference_Simulations/Study 1/Impact of h/Codes")#Set the path of this script
root_path = dirname(getwd())
setwd(paste(root_path,"/Intermediate Data",sep=""))
file_names <- dir(paste(root_path,"/Intermediate Data",sep=""))
d = length(file_names)
Record1 = Record2 = Record3 = NULL
for(i in 1:d)
{
  selects = file_names[i]#BS 200 400
  split = as.vector(str_split(selects,"_"))[[1]]
  model = as.numeric(split[3])
  penalty = split[5]
  n = as.numeric(split[7])
  h =  str_split(split[9],"\\.")[[1]]
  h = as.numeric(h[1])
  data = t(do.call(rbind,lapply(selects,read.csv)))
  s = ifelse(model==1,2,0) + ifelse(model==2,4,0) + ifelse(model==3,3,0)
  ERR = c(penalty,data[c(4,2,3,6:(7+s))])
  if(model == 1)
  {
    Record1 = rbind(Record1,ERR)
  }
  if(model == 2)
  {
    Record2 = rbind(Record2,ERR)
  }
  if(model == 3)
  {
    Record3 = rbind(Record3,ERR)
  }
}
Record1 = data.frame(Record1)
colnames(Record1)[1:4] = c("penalty","model","n","h")
Record2 = data.frame(Record2)
colnames(Record2)[1:4] = c("penalty","model","n","h")
Record3 = data.frame(Record3)
colnames(Record3)[1:4] = c("penalty","model","n","h")
#########################################Let's plot
#Record_400: save 400 for lasso and scad
#Record_lasso_200: save 200 of lasso
#Record_scad_200: save 200 of scad
for(penalty in c("lasso","SCAD"))
{
Plot_save = NULL
my_plots <- list()
i = 1
for(model in 1:3)
{
  for(n in c(200,400))
  {
    if(model==1)
    {
      Power = Record1[Record1$n==n&Record1$penalty==penalty,]
      Power = Power[order(as.numeric(Power$h)),]
      h = 1:20
      data1 = data.frame(h = h, Power = as.numeric(Power$X5), index = rep("1",20))
      data2 = data.frame(h = h, Power = as.numeric(Power$X6), index = rep("2",20))
      data3 = data.frame(h = h, Power = as.numeric(Power$X7), index = rep("t1",20))
      data4 = data.frame(h = h, Power = as.numeric(Power$X8), index = rep("t2",20))
      plotdata = rbind(data1,data2,data3,data4)
      t="I"
      plot_add =  ggplot(plotdata,aes(h,Power,group=index, color = index,shape = index,linetype = index)) + ggtitle(paste("Model ",t," (n=",n,")",sep=""))+ theme_bw()+  theme(plot.title = element_text(hjust = 0.5)) + geom_point(size = 1.5)+ geom_line(size = 0.5)
      }
    
    if(model==2)
    {
      Power = Record2[Record2$n==n&Record2$penalty==penalty,]
      Power = Power[order(as.numeric(Power$h)),]
      h = 1:20
      data1 = data.frame(h = h, Power = as.numeric(Power$X5), index = rep("1",20))
      data2 = data.frame(h = h, Power = as.numeric(Power$X6), index = rep("2",20))
      data3 = data.frame(h = h, Power = as.numeric(Power$X7), index = rep("3",20))
      data4 = data.frame(h = h, Power = as.numeric(Power$X8), index = rep("4",20))
      data5 = data.frame(h = h, Power = as.numeric(Power$X9), index = rep("t1",20))
      data6 = data.frame(h = h, Power = as.numeric(Power$X10), index = rep("t2",20))
      plotdata = rbind(data1,data2,data3,data4,data5,data6)
      t="II"
      plot_add =  ggplot(plotdata,aes(h,Power,group=index, color = index,shape = index,linetype = index)) + ggtitle(paste("Model ",t," (n=",n,")",sep=""))+ theme_bw()+  theme(plot.title = element_text(hjust = 0.5)) + geom_point(size = 1.5)+ geom_line(size = 0.5)
    }
    
    if(model==3)
    {
      Power = Record3[Record3$n==n&Record3$penalty==penalty,]
      Power = Power[order(as.numeric(Power$h)),]
      h = 1:20
      data1 = data.frame(h = h, Power = as.numeric(Power$X5), index = rep("1",20))
      data2 = data.frame(h = h, Power = as.numeric(Power$X6), index = rep("3",20))
      data3 = data.frame(h = h, Power = as.numeric(Power$X7), index = rep("p",20))
      data4 = data.frame(h = h, Power = as.numeric(Power$X8), index = rep("t1",20))
      data5 = data.frame(h = h, Power = as.numeric(Power$X9), index = rep("t2",20))
      plotdata = rbind(data1,data2,data3,data4,data5)
      t="III"
      plot_add =  ggplot(plotdata,aes(h,Power,group=index, color = index,shape = index,linetype = index)) + ggtitle(paste("Model ",t," (n=",n,")",sep=""))+ theme_bw()+  theme(plot.title = element_text(hjust = 0.5)) + geom_point(size = 1.5)+ geom_line(size = 0.5)
    }
    Plot_save = Plot_save + plot_add
    my_plots[[i]] = plot_add
    i = i + 1
  }
}
Plots = my_plots[[1]] + my_plots[[2]] + my_plots[[3]] + my_plots[[4]] + my_plots[[5]] + my_plots[[6]] + plot_layout((ncol=2))
#grid.arrange(grobs=my_plots,ncol=2)
setwd(paste(root_path,"/Result",sep=""))
ggsave(Plots, file=paste("h_impact_",penalty,".eps",sep=""), device="eps",width = 9,height = 9)
}
