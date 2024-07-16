library(ggplot2)

load("bladder_testes_data_subset_filtered2.rda")
miRNA_counts = apply(assay(bladder_testes_data_subset_filtered2), 2, function(x) by(x, rowData(bladder_testes_data_subset_filtered2)$miRNA, sum))
total_counts=colSums(assay(bladder_testes_data_subset_filtered2))

#find canonical sequences from ERCC file
ratio_pool=read.csv(file="ERCC files/FINAL_Ratiometric_SynthA_and_SynthB-1_fixlast3.csv", sep="\t")



########################################Figure 1A
miRNA_plot="hsa-let-7g-5p"

#find aggregated counts
agg_choose=data.frame(count=miRNA_counts[rownames(miRNA_counts)==miRNA_plot,])
agg_choose$total_counts=total_counts
agg_choose$cpm=agg_choose$count/(agg_choose$total_counts/1000000)
agg_choose$sequence="aggregated"
agg_choose$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue

# calculate sequence median CPM expression
idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA==miRNA_plot)
Y_all_sub=t(assay(bladder_testes_data_subset_filtered2)[idx,])
Y_seq_labels_sub=rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]
cpm_sub=Y_all_sub/(total_counts/1000000)
median_cpm_sub=apply(cpm_sub, 2, median)
logmedCPM=data.frame(sequence=Y_seq_labels_sub, logmedCPM=log(median_cpm_sub))

#find canonical sequence
ratio_pool_in=subset(ratio_pool, ratio.seqID==miRNA_plot)
data_canon=data.frame(cpm=cpm_sub[, which(Y_seq_labels_sub==ratio_pool_in$sequence)])
data_canon$sequence="canonical" #ratio_pool_in$sequence
data_canon$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue



# from remaining sequences find the top two expressing sequences
data_left=data.frame(as.matrix(cpm_sub[, which(Y_seq_labels_sub!=ratio_pool_in$sequence)]))
colnames(data_left)=Y_seq_labels_sub[which(Y_seq_labels_sub!=ratio_pool_in$sequence)]
logmedCPM_left=subset(logmedCPM, sequence %in% colnames(data_left))
next_two_seq=logmedCPM_left$sequence[order(logmedCPM_left$logmedCPM, decreasing=TRUE)[1:2]]
data_next_one=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[1])])
data_next_one$sequence="isomiR 1"#next_two_seq[1]
data_next_one$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue

data_next_two=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[2])])
data_next_two$sequence="isomiR 2"#next_two_seq[2]
data_next_two$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
data_next_two=rbind(data_next_one, data_next_two)

#boxplots
data_comb=rbind(agg_choose[, c("cpm", "sequence", "col_group")], data_canon[,c("cpm", "sequence", "col_group")],
                data_next_two[, c("cpm", "sequence", "col_group")])
data_comb$sequence=factor(data_comb$sequence, levels=unique(data_comb$sequence))
p1=ggplot(data_comb, aes(x=sequence, y=log(cpm+1)))+geom_boxplot(aes(fill=col_group))+ ylab('log(CPM+1)')+
  theme(strip.text = element_text(size=12), axis.text.x=element_text(size=15), axis.title.x=element_blank(),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.text=element_text(size=15),
        legend.title=element_text(size=15), plot.title=element_text(size=15, hjust=0.5))+
  labs(fill="tissue")+ggtitle(miRNA_plot)


#ggsave(paste0("figures/Fig2A_", miRNA_plot, ".tif"), plot=last_plot(), device="tiff", scale=3, width=60, height=50, units="mm", dpi=320, bg="white")



########################################Figure 1B
miRNA_plot="hsa-miR-26a-5p"

#find aggregated counts
agg_choose=data.frame(count=miRNA_counts[rownames(miRNA_counts)==miRNA_plot,])
agg_choose$total_counts=total_counts
agg_choose$cpm=agg_choose$count/(agg_choose$total_counts/1000000)
agg_choose$sequence="aggregated"
agg_choose$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue

# calculate sequence median CPM expression
idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA==miRNA_plot)
Y_all_sub=t(assay(bladder_testes_data_subset_filtered2)[idx,])
Y_seq_labels_sub=rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]
cpm_sub=Y_all_sub/(total_counts/1000000)
median_cpm_sub=apply(cpm_sub, 2, median)
logmedCPM=data.frame(sequence=Y_seq_labels_sub, logmedCPM=log(median_cpm_sub))

#find canonical sequence
ratio_pool_in=subset(ratio_pool, ratio.seqID==miRNA_plot)
data_canon=data.frame(cpm=cpm_sub[, which(Y_seq_labels_sub==ratio_pool_in$sequence)])
data_canon$sequence="canonical" #ratio_pool_in$sequence
data_canon$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue



# from remaining sequences find the top two expressing sequences
data_left=data.frame(as.matrix(cpm_sub[, which(Y_seq_labels_sub!=ratio_pool_in$sequence)]))
colnames(data_left)=Y_seq_labels_sub[which(Y_seq_labels_sub!=ratio_pool_in$sequence)]
logmedCPM_left=subset(logmedCPM, sequence %in% colnames(data_left))
next_two_seq=logmedCPM_left$sequence[order(logmedCPM_left$logmedCPM, decreasing=TRUE)[1:2]]
data_next_one=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[1])])
data_next_one$sequence="isomiR 1"#next_two_seq[1]
data_next_one$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue

data_next_two=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[2])])
data_next_two$sequence="isomiR 2"#next_two_seq[2]
data_next_two$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
data_next_two=rbind(data_next_one, data_next_two)

#boxplots
data_comb=rbind(agg_choose[, c("cpm", "sequence", "col_group")], data_canon[,c("cpm", "sequence", "col_group")],
                data_next_two[, c("cpm", "sequence", "col_group")])
data_comb$sequence=factor(data_comb$sequence, levels=unique(data_comb$sequence))
p2=ggplot(data_comb, aes(x=sequence, y=log(cpm+1)))+geom_boxplot(aes(fill=col_group))+ ylab('log(CPM+1)')+
  theme(strip.text = element_text(size=12), axis.text.x=element_text(size=15), axis.title.x=element_blank(),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.text=element_text(size=15),
        legend.title=element_text(size=15), plot.title=element_text(size=15, hjust=0.5))+
  labs(fill="tissue")+ggtitle(miRNA_plot)




########################################Figure 1C
miRNA_plot="hsa-let-7a-5p"

#find aggregated counts
agg_choose=data.frame(count=miRNA_counts[rownames(miRNA_counts)==miRNA_plot,])
agg_choose$total_counts=total_counts
agg_choose$cpm=agg_choose$count/(agg_choose$total_counts/1000000)
agg_choose$sequence="aggregated"
agg_choose$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue

# calculate sequence median CPM expression
idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA==miRNA_plot)
Y_all_sub=t(assay(bladder_testes_data_subset_filtered2)[idx,])
Y_seq_labels_sub=rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]
cpm_sub=Y_all_sub/(total_counts/1000000)
median_cpm_sub=apply(cpm_sub, 2, median)
logmedCPM=data.frame(sequence=Y_seq_labels_sub, logmedCPM=log(median_cpm_sub))

#find canonical sequence
ratio_pool_in=subset(ratio_pool, ratio.seqID==miRNA_plot)
data_canon=data.frame(cpm=cpm_sub[, which(Y_seq_labels_sub==ratio_pool_in$sequence)])
data_canon$sequence="canonical" #ratio_pool_in$sequence
data_canon$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue



# from remaining sequences find the top two expressing sequences
data_left=data.frame(as.matrix(cpm_sub[, which(Y_seq_labels_sub!=ratio_pool_in$sequence)]))
colnames(data_left)=Y_seq_labels_sub[which(Y_seq_labels_sub!=ratio_pool_in$sequence)]
logmedCPM_left=subset(logmedCPM, sequence %in% colnames(data_left))
next_two_seq=logmedCPM_left$sequence[order(logmedCPM_left$logmedCPM, decreasing=TRUE)[1:2]]
data_next_one=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[1])])
data_next_one$sequence="isomiR 1"#next_two_seq[1]
data_next_one$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue

data_next_two=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[2])])
data_next_two$sequence="isomiR 2"#next_two_seq[2]
data_next_two$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
data_next_two=rbind(data_next_one, data_next_two)

#boxplots
data_comb=rbind(agg_choose[, c("cpm", "sequence", "col_group")], data_canon[,c("cpm", "sequence", "col_group")],
                data_next_two[, c("cpm", "sequence", "col_group")])
data_comb$sequence=factor(data_comb$sequence, levels=unique(data_comb$sequence))
p3=ggplot(data_comb, aes(x=sequence, y=log(cpm+1)))+geom_boxplot(aes(fill=col_group))+ ylab('log(CPM+1)')+
  theme(strip.text = element_text(size=12), axis.text.x=element_text(size=15), axis.title.x=element_blank(),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.text=element_text(size=15),
        legend.title=element_text(size=15), plot.title=element_text(size=15, hjust=0.5))+
  labs(fill="tissue")+ggtitle(miRNA_plot)

library(gridExtra)
ggarrange(p1,p2,p3, ncol=1, nrow=3)

#ggsave(paste0("figures/Fig2C_", miRNA_plot, ".tif"), plot=last_plot(), device="tiff", scale=3, width=60, height=50, units="mm", dpi=320, bg="white")
ggsave("figures/figure1.tif", plot=last_plot(), device="tiff", width=12.4, height=7.3, units="in", dpi=320, bg="white")
