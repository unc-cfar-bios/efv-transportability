# Author: Steven Xu
# Modifications: KR Mollan (KM plot axes)
# Purpose: Kaplan-Meier plot of IOSW risks using MI


#library(tidyr)
#library(tibble)
#library(boot)
#library(ggplot2)
#library(gridExtra)
#library(Cairo)
#library(coxed)

library(haven)
library(survival)
library(mice)
library(parallel)
library(randomForest)
library(tictoc)


#Restricted quadratic spline function
source("G:/projects/cfar_intern/steven_xu/DR014/code/rqspline.R")

setwd("G:/projects/cfar_intern/steven_xu/DR014/data")

gen <- data.frame(read_sas("gen.sas7bdat"))

group = gen$study

group[which(group=="")] = "CNICS"

num_group = factor(group,c("A5095","A5142","A5175","A5202","CNICS"),labels=1:5)

#Reconstruct the dataset to include only variables of interest
gen_s <- data.frame(age = as.numeric(gen$AGE),cd4 = as.numeric(gen$bcd4),lrna = as.numeric(gen$lrna),
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
                    isuic_d= gen$isuic_d,
                    datsuic = gen$datsuic,datsuicw = as.numeric(gen$datsuicw), datsuicy=as.numeric(gen$datsuicy),
                    ln_datsuicy = gen$ln_datsuicy,efv_indicator = factor(gen$efv_indicator,0:1,c("EFV-free","EFV-containing")),
                    group)

gen_s$ln_datsuicy[gen_s$datsuicy == 0 ] = log(1/365.25)
gen_s$ln_suicyrd[gen_s$isuicyrd == 0 ] = log(1/365.25)

gen_rct = subset(gen_s,s=="Selected")

risktix = seq(0,192,by=48)

efv_cont = subset(gen_rct,rand_efv=="EFV-containing")

efv_fr = subset(gen_rct,rand_efv=="EFV-free")

t1_itt = efv_cont[,"isuicwk_d"]

e1_itt = efv_cont[,"isuic_d"]

t2_itt = efv_fr[,"isuicwk_d"]

e2_itt = efv_fr[,"isuic_d"]

fit1_itt = survfit(Surv(t1_itt,e1_itt)~1)

fit2_itt = survfit(Surv(t2_itt,e2_itt)~1)


n1_itt <- n2_itt <- rep(NA,length(risktix))

for(j in 1:length(risktix)){
  n1_itt[j] = max(fit1_itt$n.risk[fit1_itt$time >= risktix[j]])
  n2_itt[j] = max(fit2_itt$n.risk[fit2_itt$time >= risktix[j]])

}

cn1_itt = as.character(n1_itt)
cn2_itt = as.character(n2_itt)


mi.num = 30

impute.lvl = levels(gen_s$group)

number.lvl = length(impute.lvl)
  
column.name = c(".imp",colnames(gen_s))

impute.df = vector("list",mi.num)

for(i in 1:mi.num){
  
  impute.df[[i]] = vector("list",number.lvl)
  
}

full.predictors = colnames(gen_s)

## The following MI code takes considerable time to run (for plotting style edits skip below to load("KM_data.rda"))

tic()

for(i in 1:number.lvl){
  
  if(impute.lvl[i]=="CNICS"){
    
    predictors = c("age","cd4","lrna","cartyr","sex","race","idu","aids","hep","dep",
                   "site","risk_hetero","risk_msm","risk_other")
    
  }else{
    
    predictors = c("age","cd4","lrna","sex","race","idu","aids","hep","dep","isuic_d")
    
  }
  
  current.group <- subset(gen_s,group==impute.lvl[i])
  
  if(sum(is.na(current.group[,predictors])>0)){
    
    current.impute <- current.group[,predictors]
    
    current.remain <- rep(list(current.group[,!(full.predictors%in%predictors)]),mi.num)
    
    if(impute.lvl[i]=="CNICS"){
      
      predMat = make.predictorMatrix(current.impute)
      
      predMat["risk_hetero",c("risk_msm","risk_other")] = 0
      
      predMat["risk_msm",c("risk_hetero","risk_other")] = 0
      
      predMat["risk_other",c("risk_msm","risk_hetero")] = 0
      
      current.MICE = parlmice(current.impute,
                          method=c("","pmm",rep("",3),"rf","rf","","","rf","",rep("rf",3)),maxit = 30,cluster.seed=1234,n.imp.core=15,n.core=2,printFlag = F)
      
    }else{
      current.MICE = parlmice(current.impute,
                          method=c("",rep("pmm",2),"","rf",rep("",4),""),maxit = 20,cluster.seed=1234,n.imp.core=15,n.core=2, printFlag = F)
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

toc()

save(file="KM_data.rda",impute.df)

# ---> Skip to here if the MI step does not need updated 

load("KM_data.rda")

itt_mat = {}

for(i in 1:mi.num){

  
  comp_sub_df = rqspline(do.call("rbind",impute.df[[i]]),c("age","cd4","lrna"),k=4,equal=T,slice = ".imp")[,-1]

  #marginal probability
  marg.p = glm(s~1,data=comp_sub_df,family = binomial())$fitted.values
  
  #conditional probability
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
  
  cph_fit_itt = survfit(Surv(isuicwk_d,isuic_d) ~ rand_efv,type="kaplan-meier",weights = w,data = comp_sub_df, subset = (s=="Selected"))
  
  itt_mat = rbind(itt_mat,cph_fit_itt$surv)
  
  
}

#Pooling by Rubin's rule
itt_vec = colMeans(itt_mat)


cph_fit_itt$surv = itt_vec


tiff(file="KM_ITT_v2.tiff",
      width=6,height=5,units = "in",res=300)

par(mfrow=c(1,1),oma=c(1,1,1,1),mar=c(8,6,1,1))

# Probability of suicidal thoughts or behaviors

plot(cph_fit_itt,mark.time=F,fun="event",col = c("gray49","black"),lty=c(5,1),xmax=192,axes=F,ylim=c(0,0.05),lwd=3,
     xlab="Weeks since study entry", ylab="Probability of suicidal thoughts/behaviors") ## , main="Intention-to-treat")
axis(1,at=seq(0,192,by=48),cex=1)
axis(2,at=seq(0,0.05,by=0.01),cex=1)
myleg = c(paste("EFV-containing:",sum(na.omit(e1_itt)),"events"),paste("EFV-free:",sum(na.omit(e2_itt)),"events"))
legend(0,0.05,myleg,lty=c(1,5),lwd=3,col=c("black","gray49"),cex=1,bty="n")
mtext("No. at risk:",side=1,at=-14,line=4,adj=1,cex=1)
mtext("EFV-containing",side=1,at=-14,line=5,adj=1,cex=1)
mtext(cn1_itt,side=1,at=risktix,line=5,cex=1)
mtext("EFV-free",side=1,at=-14,line=6,adj=1,cex=1)
mtext(cn2_itt,side=1,at=risktix,line=6,cex=1)

dev.off()



