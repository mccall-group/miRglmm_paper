load("~/scratch/McCall/miRNA/miRglmm_paper_andrea/bladder_testes_data_subset_filtered2.rda")
library(SummarizedExperiment)
library(pheatmap)
######"hsa-let-7g-5p"
idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA=="hsa-let-7g-5p")

total_counts=colSums(assay(bladder_testes_data_subset_filtered2))
Y_all_sub = t(assay(bladder_testes_data_subset_filtered2)[idx, ])
Y_seq_labels_sub = rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]

## compute median CPM for each sequence
cpm_sub = Y_all_sub / (total_counts/1e6)
median_cpm_sub = apply(cpm_sub, 2, median)

idx2 = which(log(median_cpm_sub) > -1)
Y_all_sub = Y_all_sub[ ,idx2]
Y_seq_labels_sub = Y_seq_labels_sub[idx2]

corr_mat=cor(as.matrix(log(cpm_sub[,idx2]+1)))
rownames(corr_mat)=Y_seq_labels_sub

# row_order=hclust(dist(corr_mat))$order
# column_order=hclust(dist(t(corr_mat)))$order
# heatmap(corr_mat[row_order,column_order], labRow=Y_seq_labels_sub, labCol=FALSE, cexRow=0.5, Rowv=NA, Colv=NA)
color=colorRampPalette(c("navy", "white", "red"))(2001)
p1=pheatmap(corr_mat, treeheight_col = 0, treeheight_row = 0, 
            show_rownames = TRUE, color=color, breaks=seq(-1,1,0.001),
            filename = 'figures/Sup1-hsa-let-7g-5p-isomiR-corr.tiff', 
            width=9, height=4.4)

######"hsa-miR-26a-5p"
idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA=="hsa-miR-26a-5p")

total_counts=colSums(assay(bladder_testes_data_subset_filtered2))
Y_all_sub = t(assay(bladder_testes_data_subset_filtered2)[idx, ])
Y_seq_labels_sub = rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]

## compute median CPM for each sequence
cpm_sub = Y_all_sub / (total_counts/1e6)
median_cpm_sub = apply(cpm_sub, 2, median)

idx2 = which(log(median_cpm_sub) > -1)[seq(1,69,2)]
Y_all_sub = Y_all_sub[ ,idx2]
Y_seq_labels_sub = Y_seq_labels_sub[idx2]

corr_mat=cor(as.matrix(log(cpm_sub[,idx2]+1)))
rownames(corr_mat)=Y_seq_labels_sub

# row_order=hclust(dist(corr_mat))$order
# column_order=hclust(dist(t(corr_mat)))$order
# heatmap(corr_mat[row_order,column_order], labRow=Y_seq_labels_sub, labCol=FALSE, cexRow=0.5, Rowv=NA, Colv=NA)
color=colorRampPalette(c("navy", "white", "red"))(2001)
p1=pheatmap(corr_mat, treeheight_col = 0, treeheight_row = 0, 
            show_rownames = TRUE, color=color, breaks=seq(-1,1,0.001),
            filename = 'figures/Sup1-hsa-miR-26a-5p-isomiR-corr.tiff', 
            width=9, height=4.4)

######"hsa-let-7a-5p"
idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA=="hsa-let-7a-5p")

total_counts=colSums(assay(bladder_testes_data_subset_filtered2))
Y_all_sub = t(assay(bladder_testes_data_subset_filtered2)[idx, ])
Y_seq_labels_sub = rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]

## compute median CPM for each sequence
cpm_sub = Y_all_sub / (total_counts/1e6)
median_cpm_sub = apply(cpm_sub, 2, median)

idx2 = which(log(median_cpm_sub) > -1)[seq(1,143,4)]
Y_all_sub = Y_all_sub[ ,idx2]
Y_seq_labels_sub = Y_seq_labels_sub[idx2]

corr_mat=cor(as.matrix(log(cpm_sub[,idx2]+1)))
rownames(corr_mat)=Y_seq_labels_sub

# row_order=hclust(dist(corr_mat))$order
# column_order=hclust(dist(t(corr_mat)))$order
# heatmap(corr_mat[row_order,column_order], labRow=Y_seq_labels_sub, labCol=FALSE, cexRow=0.5, Rowv=NA, Colv=NA)
color=colorRampPalette(c("navy", "white", "red"))(2001)
p1=pheatmap(corr_mat, treeheight_col = 0, treeheight_row = 0, 
            show_rownames = TRUE, color=color, breaks=seq(-1,1,0.001),
            filename = 'figures/Sup1-hsa-let-7a-5p-isomiR-corr.tiff', 
            width=9, height=4.4)

