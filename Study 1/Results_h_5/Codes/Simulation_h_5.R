######################################
#This is a demo script
#You can set the parameters as the paper to reproduce the results, the computational cost is high.
######################################
#set the path of this script
#setwd("D:/Project_Revision/MF_Inference_Simulations/Study 1/Results_h_5/Interminate Data")
root_path = dirname(getwd())
source(paste(root_path,"/Codes/Functions_Study1.R",sep=""))
#####
repeatnum = 1000 #Num of replications
n= 200 # sample size
p=2000  #dimension of X
model= 1 # 1: Model I, 2: Model II and 3: Model III
h= 5 #Number of h 
penalty="lasso" # "lasso" or "SCAD"
#Define the index of parameters of interest
interest = c(1:5,1996:2000)
rho = 0.5#Covariance matrix with AR(1) structure and rho = 0.5
sigma = 1 #standard error of noise
#############
#The results are default saved in the root_path, you can set your own saving address
#setwd(dir)
#############
Record = Record0 = Record1 = Record_decor = NULL
for(q in 1:repeatnum)
{
  tryCatch({
    Result0 =  Highdim.Fast(n,p,rho,h,penalty,model,q,interest)
    Pval1_save = Result0$pvalue #pvalues of Wn
    Pval2_save = Result0$pvalue1#pvalues of TNL i.e., decorrelated score 
    Record_decor_add = ifelse(Result0$pvalue1<=0.05,1,0)
    Record_add = ifelse(Result0$pvalue<=0.05,1,0)
    Record = rbind(Record,Record_add)
    Record_decor = rbind(Record_decor, Record_decor_add)
    Record0 = rbind(Record0,c(Pval1_save,Pval2_save))
    Type1 = apply(Record,2,mean)#
    Type2 = apply(Record_decor,2,mean)#Empirical rejection rate of decorrelated score method
    Result = c(q,n,h,model,penalty,Type1,Type2)
    print(c(q,Type1))#Output of empirical rejection rates of Wn of q replications
    #Return the empirical rejection rate until q replications
    write.csv(Result,paste("Impact_Model_",model,"_penalty_",penalty,"_n_",n,"_h_",h,".csv",sep=""),row.names = FALSE,col.names = FALSE)
    #Return the Pvalues of interested variables
    write.csv(Record0,paste("Pvalue_Model_",model,"_penalty_",penalty,"_n_",n,"_h_",h,".csv",sep=""),row.names = FALSE,col.names = FALSE)
    #    write.csv(Record1,paste("Stat_Model_",model,"_n_",n,"_h_",h1,"_s_",s,"_d_",d1,"_delta_",delta,"_choice1_",choice1,"_choice2_",choice2,"_active_",active,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE)
  },error=function(e){cat("ERROR","\n")})
}
