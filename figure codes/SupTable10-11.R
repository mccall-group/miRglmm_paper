load(results_all, file="sim results/sims_N100_m2_s1_rtruncnorm13_seqDE_results.rda")

#MSE overall
MSE_mat=as.data.frame(matrix(unlist(results_all[["MSE_sim"]]), ncol=3, byrow=TRUE))
colnames(MSE_mat)=colnames(results_all[["MSE_sim"]][[1]])

#MSE by truth
out_miRglmm=data.frame(rbind(t(sapply(seq(1, length(sims)), function(row) results_all[["MSE_sim_by_truth"]][[row]][["est_logFC"]]))))
colnames(out_miRglmm)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
miRglmm_means=data.frame(t(colMeans(out_miRglmm)))
colnames(miRglmm_means)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
miRglmm_means$method="miRglmm"

out_miRglmm2=data.frame(rbind(t(sapply(seq(1, length(sims)), function(row) results_all[["MSE_sim_by_truth"]][[row]][["est_logFC_poisson"]]))))
colnames(out_miRglmm2)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
miRglmm2_means=data.frame(t(colMeans(out_miRglmm2)))
colnames(miRglmm2_means)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
miRglmm2_means$method="miRglmm poisson"

out_DESeq2=data.frame(rbind(t(sapply(seq(1, length(sims)), function(row) results_all[["MSE_sim_by_truth"]][[row]][["DESeq2_estlogFC"]]))))
colnames(out_DESeq2)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
DESeq2_means=data.frame(t(colMeans(out_DESeq2)))
colnames(DESeq2_means)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
DESeq2_means$method="DESeq2"

MSE_by_truth=rbind(miRglmm_means, miRglmm2_means, DESeq2_means)

#var by truth
out_miRglmm=data.frame(rbind(t(sapply(seq(1, length(sims)), function(row) results_all[["var_by_truth"]][[row]][["est_logFC"]]))))
colnames(out_miRglmm)=results_all[["var_by_truth"]][[1]][["true_logFC"]]
miRglmm_means=data.frame(t(colMeans(out_miRglmm)))
colnames(miRglmm_means)=results_all[["var_by_truth"]][[1]][["true_logFC"]]
miRglmm_means$method="miRglmm"

out_miRglmm2=data.frame(rbind(t(sapply(seq(1, length(sims)), function(row) results_all[["var_by_truth"]][[row]][["est_logFC_poisson"]]))))
colnames(out_miRglmm2)=results_all[["var_by_truth"]][[1]][["true_logFC"]]
miRglmm2_means=data.frame(t(colMeans(out_miRglmm2)))
colnames(miRglmm2_means)=results_all[["var_by_truth"]][[1]][["true_logFC"]]
miRglmm2_means$method="miRglmm poisson"

out_DESeq2=data.frame(rbind(t(sapply(seq(1, length(sims)), function(row) results_all[["var_by_truth"]][[row]][["DESeq2_estlogFC"]]))))
colnames(out_DESeq2)=results_all[["var_by_truth"]][[1]][["true_logFC"]]
DESeq2_means=data.frame(t(colMeans(out_DESeq2)))
colnames(DESeq2_means)=results_all[["var_by_truth"]][[1]][["true_logFC"]]
DESeq2_means$method="DESeq2"

var_by_truth=rbind(miRglmm_means, miRglmm2_means, DESeq2_means)
