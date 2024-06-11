#load(file='monocyte_exact_subset_filtered2.rda')
load(file='study 89/study89_data_subset_filtered2.rda')
exact_subset_filtered2=study89_data_subset_filtered2[, study89_data_subset_filtered2$cell_tissue=="Monocyte"]
Y_all=t(assay(exact_subset_filtered2))
total_counts=rowSums(Y_all)

#top 2 expressing sequences

cpm_sub=Y_all/(total_counts/1000000)
median_cpm_sub=apply(cpm_sub, 2, median)
median_cpm_sub_seq=median_cpm_sub
sum_seq=colSums(cpm_sub)

max_cpm=max(sum_seq)
idx_max=which(sum_seq==max_cpm)
rowData(exact_subset_filtered2)[idx_max,]

second_max_cpm=max(sum_seq[sum_seq<max_cpm])
idx_second_max=which(sum_seq==second_max_cpm) 
rowData(exact_subset_filtered2)[idx_second_max,]

library(ggplot2)
df2=data.frame(Y_all[,idx_max], Y_all[, idx_second_max])
df2=df2/total_counts
spear_corr=cor(df2[,1], df2[,2], method="spearman")
print(spear_corr)
names(df2)=c("max sequence", "second max sequence")
p2=ggplot(df2, aes(x=`max sequence`, y=`second max sequence`))+geom_point()+
  ggtitle(label='Proportion of reads mapped to top expressing isomiRs', subtitle='Spearman Correlation = 0.44')+
  xlab(paste('Proportion of reads mapped to', '\n', 'top expressing isomiR'))+
  ylab(paste('Proportion of reads mapped to second ', '\n', 'highest expressing isomiR'))+
  theme(plot.title=element_text(size=15, hjust=0.5), plot.subtitle=element_text(size=12, hjust=0.5),
        axis.title.x = element_text(size=15), axis.text.x=element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.y=element_text(size=12))
  
#ggsave("figures/figure3A.tif", plot=last_plot(), device="tiff", scale=4, width=40, height=40, units="mm", dpi=320, bg="white")



#top two miRNA corr
data_miRNA = t(apply(assay(exact_subset_filtered2), 2, function(x) by(x, rowData(exact_subset_filtered2)$miRNA, sum)))


cpm_sub=data_miRNA/(total_counts/1000000)
median_cpm_sub=apply(cpm_sub, 2,median)
sum_seq=colSums(cpm_sub)

max_cpm=max(sum_seq)
idx_max=which(sum_seq==max_cpm)
colnames(data_miRNA)[idx_max] #"hsa-miR-21-5p"


second_max_cpm=max(sum_seq[sum_seq<max_cpm])
idx_second_max=which(sum_seq==second_max_cpm) #"hsa-miR-191-5p"
colnames(data_miRNA)[idx_second_max]

df2=data.frame(data_miRNA[,idx_max], data_miRNA[, idx_second_max])
df2=df2/total_counts
spear_corr=cor(df2[,1], df2[,2], method="spearman")
print(spear_corr)
names(df2)=c("max sequence", "second max sequence")
p1=ggplot(df2, aes(x=`max sequence`, y=`second max sequence`))+geom_point()+
  ggtitle(label='Proportion of reads mapped to top expressing miRNA', subtitle='Spearman Correlation = -0.65')+
  xlab(paste('Proportion of reads mapped to', '\n', 'top expressing miRNA'))+
  ylab(paste('Proportion of reads mapped to second', '\n', 'highest expressing miRNA'))+
  theme(plot.title=element_text(size=15, hjust=0.5), plot.subtitle=element_text(size=12, hjust=0.5),
        axis.title.x = element_text(size=15), axis.text.x=element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.y=element_text(size=12))

library(gridExtra)
ggarrange(p1,p2, nrow=1, ncol=2)
ggsave("figures/supfigure4_mono.tif", plot=last_plot(), device="tiff", width=14, height=6, units="in", dpi=320, bg="white")

