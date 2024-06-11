library(ggplot2)
library(lme4)


## read in results
load(file='bladder_testes_results/filter_neg1_results.rda')

############################# significant for miRglmm only
miRNA_plot="hsa-miR-100-5p"

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
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
ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.5)+geom_point(alpha=0.5)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="NB GLM estimates"), size=2)+xlab("tissue")+
  scale_y_continuous(breaks=c(-4, -2, 0, 2), labels=round(exp(c(-4,-2,0,2)),3))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "NB GLM estimates"="blue"))+
  theme(legend.position="bottom", legend.direction="horizontal", 
        plot.title=element_text(size=20, hjust=0.5), axis.title=element_text(size=20), axis.text=element_text(size=15), legend.text=element_text(size=15))

#ggsave(paste0("figures/visabstract_Fig2_legend.tif"), plot=last_plot(), device="tiff", scale=2.5, width=83, height=50, units="mm", dpi=320, bg="white")
#ggsave(paste0("figures/visabstract_Fig2.tif"), plot=last_plot(), device="tiff", scale=3, width=62, height=50, units="mm", dpi=320, bg="white")
ggsave("figures/figure6_legend.tif", plot=last_plot(), device="tiff", width=12.8, height=4.3, units="in", dpi=320, bg="white")
  

p1=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.5)+geom_point(alpha=0.5)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+
  scale_y_continuous(breaks=log(c(10,1,0.1, 0.01)), labels=c(10,1,0.1, 0.01))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-100-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

#ggsave("figures/figure5_A.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")

############################# significant for miRglmm only
miRNA_plot="hsa-miR-143-5p"#"hsa-miR-100-5p"#

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
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

p2=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.5)+geom_point(alpha=0.5)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+
  scale_y_continuous(breaks=log(c(2,1,0.5)), labels=c(2,1,0.5))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-143-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

#ggsave("figures/figure5_B.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")

############################# significant for miRglmm only
miRNA_plot="hsa-miR-222-3p"#"hsa-miR-100-5p"#

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
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

p3=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.5)+geom_point(alpha=0.5)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+
  scale_y_continuous(breaks=log(c(1000, 100,10,1,0.1, 0.01)), labels=c(expression(paste("10"^"3")),expression(paste("10"^"2")), 10,1,0.1, 0.01))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-222-3p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

#ggsave("figures/figure5_C.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")


############################# miRglmm only method not sig
miRNA_plot="hsa-miR-25-3p"#"hsa-miR-100-5p"#

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
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

p4=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.5)+geom_point(alpha=0.5)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+
  scale_y_continuous(breaks=log(c(10,1,0.1, 0.01)), labels=c(10,1,0.1, 0.01), limits=c(log(0.09),NA))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-25-3p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

#ggsave("figures/figure5_D.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")



############################# miRglmm only method not sig
miRNA_plot="hsa-miR-423-5p"#

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
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

p5=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.5)+geom_point(alpha=0.5)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+xlab("tissue")+
  scale_y_continuous(breaks=log(c(10,1,0.1, 0.01)), labels=c(10,1,0.1, 0.01))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-423-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5),
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

#ggsave("figures/figure5_E.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")


############################# miRglmm  and agg very diff
miRNA_plot="hsa-miR-664a-5p.SNPC"#

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
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

p6=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm isomiR estimates (random effects)"))+geom_line(alpha=0.5)+geom_point(alpha=0.5)+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm estimates (fixed effects)"), size=2)+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="aggregated method estimates"), size=2)+
  scale_y_continuous(breaks=log(c(10,1,0.1, 0.01)), labels=c(10,1,0.1, 0.01), limits=c(-3, NA))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.1,0))+
  scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle("hsa-miR-664a-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(),
        axis.title=element_text(size=20), axis.text=element_text(size=20))

#ggsave("figures/figure5_F.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")

library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3)

ggsave("figures/figure6.tif", plot=last_plot(), device="tiff", width=16.8, height=9, units="in", dpi=320, bg="white")
