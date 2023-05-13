#setwd("D:/Project_Revision/MF_Inference_Simulations/FDR")#set the path of this script
source("Functions_FDR.R")
#####
repeatnum = 1000 #Do not recommend running in local desktop
n=n1= 600
p=6
model= 4 #Can select model 4 or 5
h=5
penalty = "lasso"
s= 2 #Sparsity level
Record = Record0 = Stat = NULL
#############
for(q in 1:repeatnum)
{
  tryCatch({
    seed = q  
    time0 = Sys.time()
    Result0 =  Highdim(n,p,s,rho,h,penalty,model,q)
    Pval1_save = Result0$pvalue
    Pval2_save = Result0$pvalue1
    Stat1_save = Result0$stat
    Record0 = rbind(Record0,c(Pval1_save,Pval2_save))
    Stat = rbind(Stat, Stat1_save)
    time1 = Sys.time()
    xx = Stat[,6]
    err = apply(Stat,2,function(x) ifelse(x>=qchisq(0.95,5),1,0))
    apply(err,2,mean)
    plot(density(xx))
    xfit<-seq(min(xx),max(xx),length=40)
    yfit<-dchisq(xfit,5)
    lines(xfit, yfit, col="blue", lwd=2)
    print(q) #Time counter
    #write.csv(c(q),paste("Time_Model_",model,"_penalty_",penalty,"_n_",n,"_s_",s,"_part_",part,".csv",sep=""),row.names = FALSE,col.names = FALSE)
    write.csv(Stat,paste("Stat_Model_",model,"_n_",n,"_s_",s,".csv",sep=""),row.names = FALSE,col.names = FALSE)
    #write.csv(Record0,paste("Pvalue_Model_",model,"_penalty_",penalty,"_n_",n,"_s_",s,"_part_",part,".csv",sep=""),row.names = FALSE,col.names = FALSE)
    #write.csv(Record1,paste("Stat_Model_",model,"_n_",n,"_h_",h1,"_s_",s,"_d_",d1,"_delta_",delta,"_choice1_",choice1,"_choice2_",choice2,"_active_",active,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE)
  },error=function(e){cat("ERROR","\n")})
}

