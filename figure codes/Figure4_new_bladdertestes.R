library(ggplot2)
library(lme4)


## read in results
load(file='GTex_bladder_testes_exact/results/nofilter_results.rda')

############################# significant for miRglmm only
miRNA_plot="hsa-miR-26b-5p"

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_groupTestis+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

#pull out fixed effect estimates
groupA=fixef(f1)[1]-fixef(f1)[1]
groupB=fixef(f1)[1]+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_overall=rbind(groupA, groupB)

miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
groupA=coef(miRNA_sub)[1]-coef(miRNA_sub)[1]
groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]-coef(miRNA_sub)[1]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_miRNA=rbind(groupA, groupB)
ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="NB GLM estimates"), size=2)+xlab("tissue")+theme_minimal()+
  scale_y_continuous(breaks=c(-4, -2, 0, 2), labels=round(exp(c(-4,-2,0,2)),3))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "NB GLM estimates"="blue"))+
  theme(legend.position="bottom", legend.direction="horizontal", 
        plot.title=element_text(size=20, hjust=0.5), axis.title=element_text(size=20), axis.text=element_text(size=15), legend.text=element_text(size=15))



#ggsave("figures/figure4_legend.tif", plot=last_plot(), device="tiff", width=12.8, height=4.3, units="in", dpi=320, bg="white")
  

p1=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+theme_minimal()+
  scale_y_continuous(limits=c(log(0.01), log(100)), breaks=log(c(100,10,1,0.1, 0.01, 0.001)), labels=c(expression(paste("10"^"2")), 10,1,0.1,expression(paste("10"^"-2")), expression(paste("10"^"-3"))))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-26b-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))


############################# significant for miRglmm only
miRNA_plot="hsa-let-7i-5p"#"hsa-miR-100-5p"#

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_groupTestis+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

#pull out individual sequence estimates
groupA=fixef(f1)[1]-fixef(f1)[1]
groupB=fixef(f1)[1]+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_overall=rbind(groupA, groupB)

miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
groupA=coef(miRNA_sub)[1]-coef(miRNA_sub)[1]
groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]-coef(miRNA_sub)[1]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_miRNA=rbind(groupA, groupB)

p2=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+theme_minimal()+
  scale_y_continuous(limits=c(log(0.01), log(100)),breaks=log(c(1000, 100,10,1,0.1, 0.01)), labels=c(expression(paste("10"^"3")),expression(paste("10"^"2")), 10,1,0.1, expression(paste("10"^"-2"))))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-let-7i-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))


############################# significant for miRglmm only
miRNA_plot="hsa-miR-24-3p"#"hsa-miR-100-5p"#

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_groupTestis+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

#pull out individual sequence estimates
groupA=fixef(f1)[1]-fixef(f1)[1]
groupB=fixef(f1)[1]+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_overall=rbind(groupA, groupB)

miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
groupA=coef(miRNA_sub)[1]-coef(miRNA_sub)[1]
groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]-coef(miRNA_sub)[1]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_miRNA=rbind(groupA, groupB)

p3=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+theme_minimal()+
 scale_y_continuous(breaks=log(c(1000, 100,10,1,0.1, 0.01)), labels=c(expression(paste("10"^"3")),expression(paste("10"^"2")), 10,1,0.1, expression(paste("10"^"-2"))))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-24-3p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))



############################# miRglmm only method not sig
miRNA_plot="hsa-miR-30a-5p"#"hsa-miR-100-5p"#

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_groupTestis+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

#pull out individual sequence estimates
groupA=fixef(f1)[1]-fixef(f1)[1]
groupB=fixef(f1)[1]+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_overall=rbind(groupA, groupB)

miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
groupA=coef(miRNA_sub)[1]-coef(miRNA_sub)[1]
groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]-coef(miRNA_sub)[1]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_miRNA=rbind(groupA, groupB)

p4=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+theme_minimal()+
  scale_y_continuous(breaks=log(c(1000, 100,10,1,0.1, 0.01, 0.001)), labels=c(expression(paste("10"^"3")),expression(paste("10"^"2")), 10,1,0.1,expression(paste("10"^"-2")), expression(paste("10"^"-3"))))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-30a-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))



# ##### study 46 results
# ## read in results
# load(file='study46_exact/results/filter_neg1_results.rda')
# 
# ############################# significant for miRglmm only
# miRNA_plot="hsa-miR-26b-5p"
# 
# f1=fits[["miRglmm"]][[miRNA_plot]]
# 
# #pull out individual sequence estimates
# groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
# groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
# groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
# groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df=rbind(groupA, groupB)
# 
# #pull out fixed effect estimates
# groupA=fixef(f1)[1]-fixef(f1)[1]
# groupB=fixef(f1)[1]+fixef(f1)[2]-fixef(f1)[1]
# groupA=data.frame("estimate"=groupA)
# groupB=data.frame("estimate"=groupB)
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df_overall=rbind(groupA, groupB)
# 
# miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
# groupA=coef(miRNA_sub)[1]-coef(miRNA_sub)[1]
# groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]-coef(miRNA_sub)[1]
# groupA=data.frame("estimate"=groupA)
# groupB=data.frame("estimate"=groupB)
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df_miRNA=rbind(groupA, groupB)
# 
# 
# p5=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
#   geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
#   geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+theme_minimal()+
#   scale_y_continuous(limits=c(log(0.1), log(10)), breaks=log(c(100,10,1,0.1, 0.01, 0.001)), labels=c(expression(paste("10"^"2")), 10,1,0.1,expression(paste("10"^"-2")), expression(paste("10"^"-3"))))+
#   ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
#   scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
#   ggtitle("hsa-miR-26b-5p")+
#   theme(legend.position="none",
#         plot.title=element_text(size=20, hjust=0.5), 
#         axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))
# 
# #ggsave("figures/figure5_A.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")
# 
# ############################# significant for miRglmm only
# miRNA_plot="hsa-let-7i-5p"#"hsa-miR-100-5p"#
# 
# f1=fits[["miRglmm"]][[miRNA_plot]]
# 
# #pull out individual sequence estimates
# groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
# groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
# groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
# groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df=rbind(groupA, groupB)
# 
# #pull out individual sequence estimates
# groupA=fixef(f1)[1]-fixef(f1)[1]
# groupB=fixef(f1)[1]+fixef(f1)[2]-fixef(f1)[1]
# groupA=data.frame("estimate"=groupA)
# groupB=data.frame("estimate"=groupB)
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df_overall=rbind(groupA, groupB)
# 
# miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
# groupA=coef(miRNA_sub)[1]-coef(miRNA_sub)[1]
# groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]-coef(miRNA_sub)[1]
# groupA=data.frame("estimate"=groupA)
# groupB=data.frame("estimate"=groupB)
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df_miRNA=rbind(groupA, groupB)
# 
# p6=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
#   geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
#   geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+theme_minimal()+
#   scale_y_continuous(limits=c(log(0.1), log(10)),breaks=log(c(1000, 100,10,1,0.1, 0.01)), labels=c(expression(paste("10"^"3")),expression(paste("10"^"2")), 10,1,0.1, expression(paste("10"^"-2"))))+
#   ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
#   scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
#   ggtitle("hsa-let-7i-5p")+
#   theme(legend.position="none",
#         plot.title=element_text(size=20, hjust=0.5), 
#         axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))
# 
# #ggsave("figures/figure5_B.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")
# 
# ############################# significant for miRglmm only
# miRNA_plot="hsa-miR-24-3p"#"hsa-miR-100-5p"#
# 
# f1=fits[["miRglmm"]][[miRNA_plot]]
# 
# #pull out individual sequence estimates
# groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
# groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
# groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
# groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df=rbind(groupA, groupB)
# 
# #pull out individual sequence estimates
# groupA=fixef(f1)[1]-fixef(f1)[1]
# groupB=fixef(f1)[1]+fixef(f1)[2]-fixef(f1)[1]
# groupA=data.frame("estimate"=groupA)
# groupB=data.frame("estimate"=groupB)
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df_overall=rbind(groupA, groupB)
# 
# miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
# groupA=coef(miRNA_sub)[1]-coef(miRNA_sub)[1]
# groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]-coef(miRNA_sub)[1]
# groupA=data.frame("estimate"=groupA)
# groupB=data.frame("estimate"=groupB)
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df_miRNA=rbind(groupA, groupB)
# 
# p7=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
#   geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
#   geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+theme_minimal()+
#   scale_y_continuous(breaks=log(c(1000, 100,10,1,0.1, 0.01)), labels=c(expression(paste("10"^"3")),expression(paste("10"^"2")), 10,1,0.1, expression(paste("10"^"-2"))))+
#   ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
#   scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
#   ggtitle("hsa-miR-24-3p")+
#   theme(legend.position="none",
#         plot.title=element_text(size=20, hjust=0.5), 
#         axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))
# 
# #ggsave("figures/figure5_C.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")
# 
# 
# ############################# miRglmm only method not sig
# miRNA_plot="hsa-miR-30a-5p"#"hsa-miR-100-5p"#
# 
# f1=fits[["miRglmm"]][[miRNA_plot]]
# 
# #pull out individual sequence estimates
# groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
# groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
# groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
# groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df=rbind(groupA, groupB)
# 
# #pull out individual sequence estimates
# groupA=fixef(f1)[1]-fixef(f1)[1]
# groupB=fixef(f1)[1]+fixef(f1)[2]-fixef(f1)[1]
# groupA=data.frame("estimate"=groupA)
# groupB=data.frame("estimate"=groupB)
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df_overall=rbind(groupA, groupB)
# 
# miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
# groupA=coef(miRNA_sub)[1]-coef(miRNA_sub)[1]
# groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]-coef(miRNA_sub)[1]
# groupA=data.frame("estimate"=groupA)
# groupB=data.frame("estimate"=groupB)
# groupA$Pool="Bladder"
# groupB$Pool="Testes"
# comb_df_miRNA=rbind(groupA, groupB)
# 
# p8=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.7)+geom_point(alpha=0.7)+
#   geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
#   geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+theme_minimal()+
#   scale_y_continuous(breaks=log(c(1000, 100,10,1,0.1, 0.01, 0.001)), labels=c(expression(paste("10"^"3")),expression(paste("10"^"2")), 10,1,0.1,expression(paste("10"^"-2")), expression(paste("10"^"-3"))))+
#   ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
#   scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
#   ggtitle("hsa-miR-30a-5p")+
#   theme(legend.position="none",
#         plot.title=element_text(size=20, hjust=0.5), 
#         axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

p5=p6=p7=p8=NULL
library(ggpubr)
ggarrange(p1, p5, p2, p6, p3, p7, p4, p8, nrow=4, ncol=2)
ggsave("figures/resub/figure4.tif", plot=last_plot(), device="tiff", width=6.8, height=18, units="in", dpi=320, bg="white")

ggarrange(p1, p2, p3, p4, nrow=2, ncol=2)
ggsave("figures/resub/figure4_2by2.tif", plot=last_plot(), device="tiff", width=12, height=12, units="in", dpi=320, bg="white")
