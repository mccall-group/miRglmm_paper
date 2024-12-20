library(ggplot2)
library(ggpubr)
library(lme4)
library(scales)
miRNA_plot="hsa-let-7b-5p"

#####study46 isomiRs

load(file='study46_exact/results/filter_neg1_results.rda')


f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

p1=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence))+geom_line()+geom_point()+
  #scale_y_continuous(limits=c(log(0.05), log(200)), breaks=log(c(1000, 100,10,1,0.1, 0.01, 0.001)), labels=c(expression(paste("10"^"3")), expression(paste("10"^"2")), 10,1,0.1,expression(paste("10"^"-2")), expression(paste("10"^"-3"))))+
  ylab('isomiR-level effect estimates')+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  #scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle(paste(miRNA_plot, "isomiRs"))+theme_minimal()+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))



############# study 46 estimates
## read in results
load(file='study46_exact/results/filter_neg1_processed_results.rda')
beta_hat=results[["beta_hat"]]
p_df=data.frame(results[["pvals"]])
beta_hat=beta_hat[, -c(2,4)]
p_df=p_df[, -c(2,4)]
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2])
for (ind in seq(1,dim(p_df)[2])){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
colnames(sig_df)=colnames(p_df)
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))



#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
miRNA_in$sig=t(sig_df[which(rownames(sig_df)==miRNA_plot),])
colnames(miRNA_in)=c("estimate", "method", "sig")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))

ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  #geom_point(data=comb_df[comb_df$Pool=="Testes" & comb_df[,3]==TRUE,], aes(2.015, estimate), shape="*", size=10)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  #scale_y_continuous(limits=c(0.5, 1.5), breaks=c(0.5, 1, 1.5), labels=c(0.5, 1, 1.5))+
  ylab('Expression relative to bladder')+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle(paste(miRNA_plot, "miRNA-level effect"))+theme_minimal()+
  theme(legend.text=element_text(size=16), plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

p2=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  geom_point(data=comb_df[comb_df$Pool=="Testes" & comb_df[,3]==TRUE,], aes(2.015, estimate), shape="*", size=10)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  #scale_y_continuous(limits=c(0.5, 1.5), breaks=c(0.5, 1, 1.5), labels=c(0.5, 1, 1.5))+
  ylab('Expression relative to bladder')+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle(paste(miRNA_plot, "miRNA-level effect"))+theme_minimal()+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

##################################next miRNA
miRNA_plot="hsa-let-7a-5p"

#####study46 isomiRs

load(file='study46_exact/results/filter_neg1_results.rda')



f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

p3=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence))+geom_line()+geom_point()+
  #scale_y_continuous(limits=c(log(0.05), log(200)), breaks=log(c(1000, 100,10,1,0.1, 0.01, 0.001)), labels=c(expression(paste("10"^"3")), expression(paste("10"^"2")), 10,1,0.1,expression(paste("10"^"-2")), expression(paste("10"^"-3"))))+
  ylab('isomiR-level effect estimates')+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  #scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle(paste(miRNA_plot, "isomiRs"))+theme_minimal()+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))



############# study 46 estimates
## read in results
load(file='study46_exact/results/filter_neg1_processed_results.rda')
beta_hat=results[["beta_hat"]]
p_df=data.frame(results[["pvals"]])
beta_hat=beta_hat[, -c(2,4)]
p_df=p_df[, -c(2,4)]
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2])
for (ind in seq(1,dim(p_df)[2])){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
colnames(sig_df)=colnames(p_df)
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))



#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
miRNA_in$sig=t(sig_df[which(rownames(sig_df)==miRNA_plot),])
colnames(miRNA_in)=c("estimate", "method", "sig")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))


p4=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  geom_point(data=comb_df[comb_df$Pool=="Testes" & comb_df[,3]==TRUE,], aes(2.015, estimate), shape="*", size=10)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  #scale_y_continuous(limits=c(0.5, 1.5), breaks=c(0.5, 1, 1.5), labels=c(0.5, 1, 1.5))+
  ylab('Expression relative to bladder')+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle(paste(miRNA_plot, "miRNA-level effect"))+theme_minimal()+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

##################################next miRNA
miRNA_plot="hsa-miR-146b-5p"

#####study46 isomiRs

load(file='study46_exact/results/filter_neg1_results.rda')



f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

p5=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence))+geom_line()+geom_point()+
  #scale_y_continuous(limits=c(log(0.05), log(200)), breaks=log(c(1000, 100,10,1,0.1, 0.01, 0.001)), labels=c(expression(paste("10"^"3")), expression(paste("10"^"2")), 10,1,0.1,expression(paste("10"^"-2")), expression(paste("10"^"-3"))))+
  ylab('isomiR-level effect estimates')+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  #scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle(paste(miRNA_plot, "isomiRs"))+theme_minimal()+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))



############# study 46 estimates
## read in results
load(file='study46_exact/results/filter_neg1_processed_results.rda')
beta_hat=results[["beta_hat"]]
p_df=data.frame(results[["pvals"]])
beta_hat=beta_hat[, -c(2,4)]
p_df=p_df[, -c(2,4)]
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2])
for (ind in seq(1,dim(p_df)[2])){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
colnames(sig_df)=colnames(p_df)
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))



#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
miRNA_in$sig=t(sig_df[which(rownames(sig_df)==miRNA_plot),])
colnames(miRNA_in)=c("estimate", "method", "sig")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))


p6=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  geom_point(data=comb_df[comb_df$Pool=="Testes" & comb_df[,3]==TRUE,], aes(2.015, estimate), shape="*", size=10)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  #scale_y_continuous(limits=c(0.5, 1.5), breaks=c(0.5, 1, 1.5), labels=c(0.5, 1, 1.5))+
  ylab('Expression relative to bladder')+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle(paste(miRNA_plot, "miRNA-level effect"))+theme_minimal()+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

##################################next miRNA
miRNA_plot="hsa-miR-29a-3p"

#####study46 isomiRs

load(file='study46_exact/results/filter_neg1_results.rda')



f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`-fixef(f1)[1]
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]-fixef(f1)[1]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

p7=ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence))+geom_line()+geom_point()+
  #scale_y_continuous(limits=c(log(0.05), log(200)), breaks=log(c(1000, 100,10,1,0.1, 0.01, 0.001)), labels=c(expression(paste("10"^"3")), expression(paste("10"^"2")), 10,1,0.1,expression(paste("10"^"-2")), expression(paste("10"^"-3"))))+
  ylab('isomiR-level effect estimates')+scale_x_discrete(expand=c(0.1,0))+xlab('x')+
  #scale_color_manual(name="", values=c("miRglmm estimates (fixed effects)"="red", "miRglmm isomiR estimates (random effects)"="black", "aggregated method estimates"="blue"))+
  ggtitle(paste(miRNA_plot, "isomiRs"))+theme_minimal()+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))



############# study 46 estimates
## read in results
load(file='study46_exact/results/filter_neg1_processed_results.rda')
beta_hat=results[["beta_hat"]]
p_df=data.frame(results[["pvals"]])
beta_hat=beta_hat[, -c(2,4)]
p_df=p_df[, -c(2,4)]
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2])
for (ind in seq(1,dim(p_df)[2])){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
colnames(sig_df)=colnames(p_df)
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))



#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
miRNA_in$sig=t(sig_df[which(rownames(sig_df)==miRNA_plot),])
colnames(miRNA_in)=c("estimate", "method", "sig")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))

p8=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  geom_point(data=comb_df[comb_df$Pool=="Testes" & comb_df[,3]==TRUE,], aes(2.015, estimate), shape="*", size=10)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  #scale_y_continuous(limits=c(0.5, 1.5), breaks=c(0.5, 1, 1.5), labels=c(0.5, 1, 1.5))+
  ylab('Expression relative to bladder')+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle(paste(miRNA_plot, "miRNA-level effect"))+theme_minimal()+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

#ggarrange(p1, p3, p5, p7, p2, p4, p6, p8, nrow=2, ncol=4)
ggarrange(p1, p3, p5, p2, p4, p6, nrow=2, ncol=3)
ggsave(paste0("figures/resub/miRNA_high_rand_slope.tif"), plot=last_plot(), 
       device="tiff", width=25, height=12, units="in", dpi=320, bg="white")
