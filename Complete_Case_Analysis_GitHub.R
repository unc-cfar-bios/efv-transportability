# Authors: Steven Xu & KR Mollan
# Purpose: Conduct complete-case IOSW analyses

library(haven)
library(tidyr)
library(tibble)
library(survival)
library(tictoc)
library(boot)
library(ggplot2)
library(gridExtra)
library(Cairo)
library(coxed)

# Restricted quadratic spline function
source("G:/projects/cfar_intern/steven_xu/DR014/code/rqspline.R")

setwd("G:/projects/cfar_intern/steven_xu/DR014/data")

gen <- data.frame(read_sas("gen.sas7bdat"))

gen = subset(gen, Zobs==1)

group = gen$study

group[which(group=="")] = "CNICS"

num_group = factor(group,c("A5095","A5142","A5175","A5202","CNICS"),labels=1:5)

#Reconstruct the dataset to include only variables of interest
gen_s <- data.frame(age = as.numeric(gen$AGE),cd4 = as.numeric(gen$bcd4),lrna = as.numeric(gen$lrna),
                    sex = factor(gen$SEX,1:2,labels = c("Male","Female")),
                    race = factor(gen$raceth2,c(1:3,9),c("White","Black","Hispanic","Other")),
                    idu = factor(gen$ivdrugb,0:1,c("No","Yes")), aids = factor(gen$hxaids,0:1,c("No","Yes")),
                    hep = factor(gen$hep,0:1,c("Negative/Missing/Indeterminate","Positive")),
                    dep = factor(gen$hxdeprx,0:1,c("No","Yes")),rand_efv = factor(gen$rand_efv,0:1,c("EFV-free","EFV-containing")),
                    s = factor(gen$S,0:1,c("No","Selected")),study = factor(gen$study,c("A5095","A5142","A5175","A5202")), z = gen$Zobs,
                    isuicyrd = as.numeric(gen$isuicyr_d),ln_suicyrd = gen$ln_suicyrd,isuicwk_d=as.numeric(gen$isuicwk_d),isuic_d=as.numeric(gen$isuic_d),
                    datsuic = as.numeric(gen$datsuic),datsuicw = as.numeric(gen$datsuicw), datsuicy=as.numeric(gen$datsuicy),
                    ln_datsuicy = gen$ln_datsuicy,efv_indicator = factor(gen$efv_indicator,0:1,c("EFV-free","EFV-containing")))

gen_s$ln_datsuicy[gen_s$datsuicy == 0 ] = log(1/365.25)
gen_s$ln_suicyrd[gen_s$isuicyrd == 0 ] = log(1/365.25)

tic("Keep track of time")
complete.analysis = function(data,indices){
  
    library(survival)
    library(survey)
    source("G:/projects/cfar_intern/steven_xu/DR014/code/rqspline.R")
  
  
    current.df = data[indices,]
    
    current.df = rqspline(current.df,c("age","cd4","lrna"),k=4,equal=T)
    
    #marginal probability
    marg.p = glm(s~1,data=current.df,family = binomial())$fitted.values
    
    #conditional probability
    cond.p = glm(s~sex+race+age+cd4+lrna+aids+idu+dep+hep+age1+age2+age3+cd41+cd42+cd43+lrna1+lrna2+lrna3+
                   sex:race+sex:age+sex:cd4+sex:lrna+sex:aids+sex:idu+sex:dep+sex:hep+
                   race:age+race:cd4+race:lrna+race:aids+race:idu+race:dep+race:hep+
                   age:cd4+age:lrna+age:aids+age:idu+age:dep+age:hep+
                   cd4:lrna+cd4:aids+cd4:idu+cd4:dep+cd4:hep+
                   lrna:aids+lrna:idu+lrna:dep+lrna:hep+
                   aids:idu+aids:dep+aids:hep+
                   idu:dep+idu:hep+
                   dep:hep,data=current.df,family = binomial())$fitted.values
    

    w = as.numeric((current.df$s=="Selected")*(marg.p/(1-marg.p))/(cond.p/(1-cond.p)))
    
    logw = as.numeric(ifelse(current.df$s=="Selected",log(w),0))
    
    current.df = data.frame(current.df,marg.p,cond.p,w,logw)
    
    design.ps = svydesign(ids=~1, weights=~w, strata=~study, data=subset(current.df, s == "Selected"))
    
    poi_fit = svyglm(isuic_d ~ rand_efv + offset(ln_suicyrd), design=design.ps, family=quasipoisson())
    
    cph_fit = coxph(Surv(isuicwk_d,isuic_d) ~ rand_efv+strata(study),ties = "efron",weights = w, data = subset(current.df, s == "Selected"))
    
    v.ir.t = matrix(c(1,1),nrow=1)
    
    v.ir.c = matrix(c(1,0),nrow=1)
    
    result.vec = numeric(4)
    
    result.vec[1] = v.ir.t%*%poi_fit$coefficients
    
    result.vec[2] = v.ir.c%*%poi_fit$coefficients
    
    result.vec[3] = exp(result.vec[1])-exp(result.vec[2])
    
    result.vec[4] = cph_fit$coefficients #est of log HR
  
    return(result.vec)
}

#Stratified bootstrap wrapper

set.seed(1234,kind = "L'Ecuyer-CMRG")

boot.out <- boot(data=gen_s,complete.analysis,R=15000,strata=num_group,parallel = "snow",ncpus = 3) #strata = 5 levels

boot.mat = boot.out$t

colnames(boot.mat) <-  c("LIR_t","LIR_c","IRD","LHR")

#save(boot.mat,file="complete_boot.rda")
#load(file = "complete_boot.rda")

result.df = data.frame(boot.mat)


#95% percentile interval
round(exp(quantile(result.df$LIR_t,c(0.025,0.975)))*1000, 1)  #IR 95%CI, EFV-cont

round(exp(quantile(result.df$LIR_c,c(0.025,0.975)))*1000, 1)  #IR 95%CI, EFV-free

round(quantile(result.df$IRD,c(0.025,0.975))*1000, 1) #IRD 95%CI

exp(quantile(result.df$LHR,c(0.025,0.975))) #HR 95%CI

toc()