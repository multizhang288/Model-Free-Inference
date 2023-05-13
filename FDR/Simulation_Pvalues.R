#setwd("D:/Project_Revision/MF_Inference_Simulations/FDR")#set the path of this script
source("Functions_FDR.R")
#####
repeatnum = 1000 #Do not recommend running in local desktop
n=n1= 400
p=20
model= 4 #Can select model 4 or 5
h=5 
penalty = "lasso" #SCAD or lasso
s= 4 #Sparsity level
Record = Record0 = Stat = NULL
rho = 0.5#correlation coefficient in AR(1) Sigma
sigma = 1# std of noise term 
#############
for(q in 1:repeatnum)
{
  tryCatch({
    seed = q  
    time0 = Sys.time()
    Result0 =  Highdim(n,p,s,rho,h,penalty,model,q)
    Pval1_save = Result0$pvalue#pvalues of our proposed method
  #  Pval2_save = Result0$pvalue1#pvalues of decorrelated score method
    Stat1_save = Result0$stat #statistics of our proposed method
  #  Record0 = rbind(Record0,c(Pval1_save,Pval2_save))
    Stat = rbind(Stat, Stat1_save)
    time1 = Sys.time()
    print(q) #Time counter
    #write.csv(c(q),paste("Time_Model_",model,"_penalty_",penalty,"_n_",n,"_s_",s,"_part_",part,".csv",sep=""),row.names = FALSE,col.names = FALSE) #Output of the number of finished repetitions. 
    write.csv(Stat,paste("Stat_Model_",model,"_n_",n,"_s_",s,".csv",sep=""),row.names = FALSE,col.names = FALSE)#Output of the statistics of finished repetitions.
    #write.csv(Record0,paste("Pvalue_Model_",model,"_penalty_",penalty,"_n_",n,"_s_",s,"_part_",part,".csv",sep=""),row.names = FALSE,col.names = FALSE)
  },error=function(e){cat("ERROR","\n")})
}

