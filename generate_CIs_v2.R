#################################################
# does: produce bootstrap estimates for 95% CI
#       for IR, IR difference, and HR
# by: Steven Xu, Ann Marie Weideman
# date created: 05/20/20
# date modified: 03/26/21
# modifications: modify figures for publication
#################################################

#################################################
# Import data
#################################################
batch_comb = read.csv("G:/projects/cfar/Mollan/AnnW/efv_transport_paper/data/batch_comb_withseed_06_21_20.csv")

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

write.csv(df_CI,"G:/projects/cfar/Mollan/AnnW/efv_transport_paper/output/bootstrap_ci_withseed.csv")

#################################################
# Generate ggplots
#################################################

require(ggplot2)

#---------------------------
#Incidence Rate Difference
#---------------------------

gg_ird<-ggplot(data=batch_comb)+
   geom_histogram(aes(x=IRD,y=..ncount..),bins=30,fill="#0072B2",
                  color="#0072B2",alpha=0.4)+
   labs(x="\nIncidence Rate Difference",y="Density\n")+
   geom_vline(xintercept=0, linetype="dashed", colour="#D55E00", size=1)+
   xlim(-0.01, 0.02)+
   theme(plot.margin = unit(c(1, 1, 2, 1), "lines"), #margin(t, r, b, l)
         axis.text=element_text(colour="black",size=12),
         axis.title=element_text(colour="black", size=12),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black",size=1) #weight axes
   )+ 
   coord_cartesian(clip = "off") #turn off margin clipping so A) is visible

#---------------------------
# Natural-log Hazard Ratio
#---------------------------

gg_lhr<-ggplot(data=batch_comb)+
   geom_histogram(aes(x=LHR,y=..ncount..),bins=30,fill="#0072B2",
                  color="#0072B2",alpha=0.4)+
   labs(x="\nNatural-log Hazard Ratio",y="Density\n")+
   geom_vline(xintercept=0, linetype="dashed", colour="#D55E00", size=1)+
   xlim(-2, 3)+
   theme(plot.margin = unit(c(1, 1, 2, 1), "lines"), #margin(t, r, b, l)
         axis.text=element_text(colour="black",size=12),
         axis.title=element_text(colour="black", size=12),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black",size=1) #weight axes
         )+ 
   coord_cartesian(clip = "off") #turn off margin clipping so B) is visible

#################################################
# Export as panel pdf and individual eps
#################################################

# Panel pdf
library(cowplot)
pdf("G:/projects/cfar/Mollan/AnnW/efv_transport_paper/output/Figure3.pdf",
    width=11, height=4)
plot_grid(gg_ird, gg_lhr, labels = c("A)","B)"),nrow=1, ncol=2,
          label_fontface="plain",label_fonfamily="serif",
          hjust=0,vjust=1.5, size=12, greedy=T, scale=0.95)
dev.off()

#Individual eps
cairo_ps("G:/projects/cfar/Mollan/AnnW/efv_transport_paper/output/Figure3A.eps",
         width=5.5, height=4,fallback_resolution = 300)
plot_grid(gg_ird, labels = "A)",nrow=1, ncol=1,
          label_fontface="plain",label_fonfamily="serif",
          hjust=0, vjust=1.5, size=12, greedy=T, scale=0.95)
dev.off()

cairo_ps("G:/projects/cfar/Mollan/AnnW/efv_transport_paper/output/Figure3B.eps",
         width=5.5, height=4,fallback_resolution = 300)
plot_grid(gg_lhr, labels = "B)",nrow=1, ncol=1,
          label_fontface="plain",label_fonfamily="serif",
          hjust=0, vjust=1.5, size=12, greedy=T, scale=0.95)
dev.off()
