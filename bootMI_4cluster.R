# Bootstrap Multiple Imputation + Transporting effect using Inverse Odds of Sampling Weights
#
# Author: Steven Xu, Reviewer: KR Mollan
# Edited to run on a cluster, Ann Marie Weideman
#
# Main analysis code that performs bootstrap multiple imputation to construct nonparamteric
# confidence intervals for hazard ratio (HR), incidence rate (IR) and incidence rate differ-
# ence (IRD) in the target population. Also included codes for calculating inverse odds of 
# sampling weights to transport trial effect to a target population. 
#
# The histograms of HR and IRD are plotted in a separate script.

library(haven) #for importing SAS file
library(tidyr) #for piping (e.g., %in%)
library(mice) #for mice() function
library(boot) #for boot() function
library(survey) #for svydesign() and svyglm()
library(survival) #for coxph() and Surv()

#Read in the combined dataset (cohort & trials stacked)
gen = data.frame(read_sas("/PATH_HERE/gen.sas7bdat"))

#Some preprocessing
group = gen$study

group[which(group=="")] = "CNICS"

gen_s = data.frame(age = as.numeric(gen$AGE),cd4 = as.numeric(gen$bcd4),lrna = as.numeric(gen$lrna),
                   sex = factor(gen$SEX,1:2,labels = c("Male","Female")),
                   race = factor(gen$raceth2,c(1:3,9),c("White","Black","Hispanic","Other")),
                   idu = factor(gen$ivdrugb,0:1,c("No","Yes")), aids = factor(gen$hxaids,0:1,c("No","Yes")),
                   hep = factor(gen$hep,0:1,c("Negative/Missing/Indeterminate","Positive")),
                   dep = factor(gen$hxdeprx,0:1,c("No","Yes")), cartyr = as.numeric(gen$cartyr), site = factor(gen$site,2:8), 
                   risk_hetero = factor(gen$risk_hetero,0:1,c("No","Yes")), risk_msm = factor(gen$risk_msm,0:1,c("No","Yes")), 
                   risk_other = factor(gen$risk_other,0:1,c("No","Yes")),
                   risk_hier = factor(gen$risk_hier,1:4),rand_efv = factor(gen$rand_efv,0:1,c("EFV-free","EFV-containing")),
                   s = factor(gen$S,0:1,c("No","Selected")),study = factor(gen$study,c("A5095","A5142","A5175","A5202")),
                   isuicyrd = as.numeric(gen$isuicyr_d),ln_suicyrd = gen$ln_suicyrd,isuicwk_d=as.numeric(gen$isuicwk_d),
                   isuic_d= as.numeric(gen$isuic_d),
                   datsuic = as.numeric(gen$datsuic),datsuicw = as.numeric(gen$datsuicw), datsuicy=as.numeric(gen$datsuicy),
                   ln_datsuicy = gen$ln_datsuicy,efv_indicator = factor(gen$efv_indicator,0:1,c("EFV-free","EFV-containing")),
                   group)

#Modify the log-time of subjects with zero follow-up
gen_s$ln_datsuicy[gen_s$datsuicy == 0 ] = log(1/365.25)
gen_s$ln_suicyrd[gen_s$isuicyrd == 0 ] = log(1/365.25)

#Setting up column names for multiple imputation
column.name = c(".imp",colnames(gen_s))

#Strata for bootstrap
num_group = factor(group,c("A5095","A5142","A5175","A5202","CNICS"),labels=1:5)

#Specify bootstrap sample size (jobarrays will run this 20x each)
B = 5

impute.lvl = levels(gen_s$group)

number.lvl = length(impute.lvl)

#Multiple imputation function (that will be input to boot function)
mult_impute = function(data,indices){
  
  source("/PATH_HERE/rqspline.R")
  
  current.df = data[indices,]
  
  column.name = c(".imp",colnames(data))
  
  impute.lvl = c("A5095","A5142","A5175","A5202","CNICS")
  
  number.lvl = length(impute.lvl)
  
  #Specify imputation iterations
  mi.num = 30
  
  impute.df = vector("list",mi.num)
  
  for(i in 1:mi.num){
    
    impute.df[[i]] = vector("list",number.lvl)
    
  }
  
  full.predictors = colnames(current.df)
  
  for(i in 1:number.lvl){
    
    if(impute.lvl[i]=="CNICS"){
      
      #Predictors to impute CNICS covariates
      predictors = c("age","cd4","lrna","cartyr","sex","race","idu","aids","hep","dep",
                     "site","risk_hetero","risk_msm","risk_other")
      
    }else{
      
      #Predictors to impute trial covariates
      predictors = c("age","cd4","lrna","sex","race","idu","aids","hep","dep","isuic_d")
      
    }
    
    current.group <- subset(current.df,group==impute.lvl[i])
    
    if(sum(is.na(current.group[,predictors])>0)){
      
      current.impute <- current.group[,predictors]
      
      current.remain <- rep(list(current.group[,!(full.predictors%in%predictors)]),mi.num)
      
      if(impute.lvl[i]=="CNICS"){
        
        predMat = make.predictorMatrix(current.impute)
        
        predMat["risk_hetero",c("risk_msm","risk_other")] = 0
        
        predMat["risk_msm",c("risk_hetero","risk_other")] = 0
        
        predMat["risk_other",c("risk_msm","risk_hetero")] = 0
        
        #MICE using PMM and random forest
        current.MICE = mice(current.impute,m=mi.num,
      	                      method=c("","pmm",rep("",3),"rf","rf","","","rf","",rep("rf",3)), maxit = 20, printFlag = F)
      }else{
        
        current.MICE = mice(current.impute,m=mi.num,
                            method=c("",rep("pmm",2),"","rf",rep("",4),""), maxit = 20, printFlag = F)        
      }
      
      
      current.complete = complete(current.MICE,"long")[,-2]
      
      for(j in 1:mi.num){
        
        impute.df[[j]][[i]] = subset(data.frame(current.complete,do.call("rbind",current.remain)),.imp==j)[,column.name]
        
      }
      
    }else{
      
      for(j in 1:mi.num){
        
        impute.df[[j]][[i]] = data.frame(.imp=j,current.group)[,column.name]
        
      }
    }
  }
  
  ans_mat = matrix(0,nrow=mi.num,ncol=4)
  
  for(i in 1:mi.num){
    
    #Construct spline variables
    comp_sub_df = rqspline(do.call("rbind",impute.df[[i]]),c("age","cd4","lrna"),k=4,equal=T,slice = ".imp")[,-1]
    
    #Marginal probability
    marg.p = glm(s~1,data=comp_sub_df,family = binomial())$fitted.values
    
    #Conditional probability
    cond.p = glm(s~sex+race+age+cd4+lrna+aids+idu+dep+hep+age1+age2+age3+cd41+cd42+cd43+lrna1+lrna2+lrna3+
                   sex:race+sex:age+sex:cd4+sex:lrna+sex:aids+sex:idu+sex:dep+sex:hep+
                   race:age+race:cd4+race:lrna+race:aids+race:idu+race:dep+race:hep+
                   age:cd4+age:lrna+age:aids+age:idu+age:dep+age:hep+
                   cd4:lrna+cd4:aids+cd4:idu+cd4:dep+cd4:hep+
                   lrna:aids+lrna:idu+lrna:dep+lrna:hep+
                   aids:idu+aids:dep+aids:hep+
                   idu:dep+idu:hep+
                   dep:hep,data=comp_sub_df,family = binomial())$fitted.values
    
    w = as.numeric((comp_sub_df$s=="Selected")*(marg.p/(1-marg.p))/(cond.p/(1-cond.p)))
    
    logw = as.numeric(ifelse(comp_sub_df$s=="Selected",log(w),0))
    
    comp_sub_df = data.frame(comp_sub_df,marg.p,cond.p,w,logw)
    
    #Weighted stratified Cox regression for estimating HR
    cph_fit = coxph(Surv(isuicwk_d,isuic_d) ~ rand_efv+strata(study),ties = "efron",weights = w,data = comp_sub_df, subset = (s=="Selected"))
    
    #Weighted stratified Poisson regreesion for estimating IR
    design.ps = svydesign(ids=~1, weights=~w, strata=~study, data=subset(comp_sub_df, s == "Selected"))
    
    poi_fit = svyglm(isuic_d ~ rand_efv + offset(ln_suicyrd), design=design.ps, family=quasipoisson())
    
    v.ir.t = matrix(c(1,1),nrow=1)
    
    v.ir.c = matrix(c(1,0),nrow=1)
    
    ans_mat[i,1] = v.ir.t%*%poi_fit$coefficients
    
    ans_mat[i,2] = v.ir.c%*%poi_fit$coefficients
    
    ans_mat[i,3] = exp(ans_mat[i,1])-exp(ans_mat[i,2])
    
    ans_mat[i,4] = cph_fit$coefficients #est of log HR
    
  }
  
  #Pooling by Rubin's rule
  result.vec= colMeans(ans_mat)
  
  return(result.vec)
  
}

#Generate a random seed with **very** low prob of replication in future runs
seedling <- as.integer(runif(1)*2e9)
set.seed(seedling)

#Stratified bootstrap with parallel computing
boot.out <- boot(data=gen_s,mult_impute,R=B,strata=num_group,ncpus=1)

batch_comb = cbind(boot.out$t,rep(seedling,B))

colnames(batch_comb) <- c("LIR_t","LIR_c","IRD","LHR","Seed")

batch_comb = as.data.frame(batch_comb)

#!file.exists prevents column names from repeating every 5 bootstraps
write.table(as.matrix(batch_comb),file="/PATH_HERE/bootMI_output.csv",sep=",",row.names=F,col.names=!file.exists("/PATH_HERE/bootMI_output.csv"),append=T)
