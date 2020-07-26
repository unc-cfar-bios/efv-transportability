#################################################
# does: produce bootstrap estimates for 95% CI
#       for IR, IR difference, and HR
# by: Steven Xu, Ann Marie Weideman
# date: 05/20/20
#################################################

#################################################
#import data
#################################################
batch_comb = read.csv("PATH_HERE/bootMI_output.csv")

#################################################
# 95% CIs
#################################################

#95th percentile interval for:
#Incidence rate in treatment group

CI_LIR_t<-exp(quantile(batch_comb$LIR_t,c(0.025,0.975)))*1000

#Incidence rate in control group

CI_LIR_c<-exp(quantile(batch_comb$LIR_c,c(0.025,0.975)))*1000

#Incidence rate difference

CI_IRD<-quantile(batch_comb$IRD,c(0.025,0.975))*1000

#Hazard ratio

CI_HR<-exp(quantile(batch_comb$LHR,c(0.025,0.975)))

df_CI<-data.frame(CI_LIR_t,CI_LIR_c,CI_IRD, CI_HR)

write.csv(df_CI,"PATH_HERE/CI_output.csv")

#################################################
#Plotting
#################################################

require(ggplot2)

tiff(file="PATH_HERE/IRD.tiff",
     width=4.5,height=4,units = "in",res=300)

ggplot(data=batch_comb)+geom_histogram(aes(x=IRD,y=..density..),bins=30,fill="Steel Blue",alpha=0.5)+
   labs(x=NULL,y="Density")+geom_vline(xintercept = 0,colour="red")+theme_bw()+
   scale_x_continuous(breaks=c(-0.005,0,0.005,0.01,0.015))+theme(axis.text=element_text(size=12))

dev.off()

tiff(file="PATH_HERE/HR.tiff",
     width=4.5,height=4,units = "in",res=300)

ggplot(data=batch_comb)+geom_histogram(aes(x=LHR,y=..density..),bins=30,fill="Steel Blue",alpha=0.5)+
   labs(x=NULL,y="Density")+geom_vline(xintercept = 0,colour="red")+theme_bw()+theme(axis.text=element_text(size=12))
 
dev.off()
