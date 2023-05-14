######################################
#This is a demo script
#You can set the parameters as the paper to reproduce the results, the computational cost is high.
######################################
#Set the path of this script
setwd("D:/Project_Revision/MF_Inference_Simulations/Study 2/Distribution")
root_path = dirname(getwd())
source(paste(root_path,"/Function-Study2.R",sep=""))
#####
repeatnum = 1000 #number of simulations
n=n1= 200
p=20
CovMatrix = matrix(rep(0,p^2),p,p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    CovMatrix[i,j] = rho^(abs(i-j))  
  }
}
penalty = "lasso"
model= 1 #1: model I', 2: model II', 3: model III', 4: model IV' and 5: model V'
s= 10#sparsity level
hfix= -1 #if hfix !=-1, you can choose h, otherwise h is fixed as 7 when model is 1-5.
h=7
delta=1 #signal strength
#Dist= "MVN" # "MVN","MSN","Mixture","Tdist","Chi"
active="random" #random or consecutive. random means the indices of active variables are evenly distributed.
Trans ="Bspline" #Poly, SIR, SIR2 or Bspline
Inactive = c(19)#The indices of inactive variables.
sigma = 1#std of noise
Record = Record0 = Record1 = Record_decor = NULL
#####
for(Dist in c("MVN","MSN","Mixture","Tdist","Chi"))
{
for(q in 1:repeatnum)
{
  tryCatch({
    seed = q  
    time0 = Sys.time()
    Result0 = Highdim.Study2(n,p,h,hfix,s,d,Dist,Trans,active,penalty,CovMatrix,model,delta,Inactive,seed)
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
    Result = c(q,hfix,n,model,delta,Dist,Trans,active,s,d1,h1,penalty,Type1,Type2)
    print(c(q,Dist,time1-time0,Type1))
    #Output of empirical rejection rate until q-th replications
    write.csv(Result,paste("Impact_Model_",model,"_n_",n,"_h_",h1,"_s_",s,"_d_",d1,"_delta_",delta,"_Dist_",Dist,"_Trans_",Trans,"_active_",active,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE)
    #Output of pvalues of parameter of interest until q-th replications
    write.csv(Record0,paste("Pvalue_Model_",model,"_n_",n,"_h_",h1,"_s_",s,"_d_",d1,"_delta_",delta,"_Dist_",Dist,"_Trans_",Trans,"_active_",active,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE) 
    #Output of statistics of parameter of interest until q-th replications
    write.csv(Record1,paste("Stat_Model_",model,"_n_",n,"_h_",h1,"_s_",s,"_d_",d1,"_delta_",delta,"_Dist_",Dist,"_Trans_",Trans,"_active_",active,"_end.csv",sep=""),row.names = FALSE,col.names = FALSE)
    
  },error=function(e){cat("ERROR","\n")})
}
}
