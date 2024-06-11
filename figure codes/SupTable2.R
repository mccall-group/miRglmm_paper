
load("bladder_testes_data_subset_filtered2.rda")
miRNA_counts = apply(assay(bladder_testes_data_subset_filtered2), 2, function(x) by(x, rowData(bladder_testes_data_subset_filtered2)$miRNA, sum))
total_counts=colSums(assay(bladder_testes_data_subset_filtered2))



samples_plot=c("SRR333658", "SRR333672", "SRR333674", "SRR333678", "SRR333680", "SRR333682")
miRNA_plot="hsa-miR-26a-5p"

idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA==miRNA_plot)
idx_col=which(bladder_testes_data_subset_filtered2$sample_id %in% samples_plot)
Y_out=assay(bladder_testes_data_subset_filtered2)[idx,idx_col]
idx_keep=which(rowSums(Y_out)>0)
Y_out=data.frame(as.matrix(Y_out[idx_keep,]))
colnames(Y_out)=samples_plot
rownames(Y_out)=rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx_keep]
miRNA_counts_in=data.frame('aggregated counts'=miRNA_counts[idx_col])
miRNA_counts_in=t(miRNA_counts_in)
colnames(miRNA_counts_in)=samples_plot
Y_out=rbind(Y_out, miRNA_counts_in)