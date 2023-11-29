######################################
#This is a demo script of the adversarial examples
#You can set the parameters as the paper to reproduce the results, the computational cost is high.
######################################
#Set the path of this script
#setwd("D:/Project_Revision/MF_Inference_Simulations/Study 2/Distribution")
#root_path = dirname(getwd())
source("Function-Adversarial_MGT.R")#
#####
repeatnum = 1000 #number of simulations
n=n1= 400
p=2000
penalty = "lasso"
hfix= -1 #if hfix !=-1, you can choose h, otherwise h is fixed as 7 when model is 1-5.
h=7
delta=1 #signal strength
df= 0 # 
Trans ="Bspline" #Poly, SIR, SIR2 or Bspline
Inactive = c(5, 365)#The indices of inactive variables.
sigma = 1#std of noise
Record = Record0 = Record1 = Record_decor = NULL
#####
for(q in 1:repeatnum)
{
  tryCatch({
    seed = q  
    time0 = Sys.time()
    Result0 = Highdim.Ad.MGT(n,p,h,hfix,df,Trans,penalty,delta,Inactive,seed)
    Pvalue_save = Result0$pvalue #pvalus of Wn
    Stat_save = Result0$stat #statistics of Wn
    Decor_save = Result0$pvalue.decor #pvalues of decorrelated score
    Record_decor_add = ifelse(Result0$pvalue.decor<=0.05,1,0)
    Record_add = ifelse(Result0$pvalue<=0.05,1,0)
    Record = rbind(Record,Record_add)
    Record0 = rbind(Record0, Pvalue_save)
    Record1 = rbind(Record1, Stat_save)
    Record_decor = rbind(Record_decor, Record_decor_add)
    time1 = Sys.time()
    Type1 = apply(Record,2,mean)#
    Type2 = apply(Record_decor,2,mean)
    h1 = Result0$h
    d1 = Result0$d
    Result = c(q,hfix,n,p,delta,df,Trans,d1,h1,penalty,Type1,Type2)
    print(c(q,df,time1-time0,Type1))
    print(c(q,df,time1-time0,Type2))
    #Output of empirical rejection rate until q-th replications
    write.csv(Result,paste("Impact_Adversarial_","_n_",n,"_h_",h1,"_d_",d1,"_delta_",delta,"_df_",df,"_Trans_",Trans,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE)
    #Output of pvalues of parameter of interest until q-th replications
    write.csv(Record0,paste("PvalueAdversarial_","_n_",n,"_h_",h1,"_d_",d1,"_delta_",delta,"_df_",df,"_Trans_",Trans,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE) 
    #Output of statistics of parameter of interest until q-th replications
    write.csv(Record1,paste("Stat_Adversarial_","_n_",n,"_h_",h1,"_d_",d1,"_delta_",delta,"_df_",df,"_Trans_",Trans,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE)
    
  },error=function(e){cat("ERROR","\n")})
}

