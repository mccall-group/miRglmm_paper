load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))

TPR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]!=0),]<0.05)/length(which(x[,9]!=0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))
auc_sim=lapply(results[["FDR"]], function(x) t(data.frame(unlist(lapply(x, function(y) auc(roc(as.numeric(x[,9]!=0),y)))))))

TPR_mat=as.data.frame(matrix(unlist(TPR_sim), ncol=9, byrow=TRUE))
colnames(TPR_mat)=colnames(results[["pvals"]][[1]])
TPR_mat=TPR_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]

auc_mat=as.data.frame(matrix(unlist(auc_sim), ncol=9, byrow=TRUE))
colnames(auc_mat)=colnames(results[["pvals"]][[1]])
auc_mat=auc_mat[, 1:8]

mean_TPR=data.frame("100%"=colMeans(TPR_mat))
mean_TNR=data.frame("100%"=colMeans(TNR_mat))
mean_auc=data.frame("100%"=colMeans(auc_mat))

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults_75perlibsize.rda")

results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))

TPR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]!=0),]<0.05)/length(which(x[,9]!=0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))
auc_sim=lapply(results[["FDR"]], function(x) t(data.frame(unlist(lapply(x, function(y) auc(roc(as.numeric(x[,9]!=0),y)))))))

TPR_mat=as.data.frame(matrix(unlist(TPR_sim), ncol=9, byrow=TRUE))
colnames(TPR_mat)=colnames(results[["pvals"]][[1]])
TPR_mat=TPR_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]

auc_mat=as.data.frame(matrix(unlist(auc_sim), ncol=9, byrow=TRUE))
colnames(auc_mat)=colnames(results[["pvals"]][[1]])
auc_mat=auc_mat[, 1:8]

mean_TPR=cbind(data.frame("75%"=colMeans(TPR_mat)), mean_TPR)
mean_TNR=cbind(data.frame("75%"=colMeans(TNR_mat)), mean_TNR)
mean_auc=cbind(data.frame("75%"=colMeans(auc_mat)), mean_auc)

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults_50perlibsize.rda")

results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))

TPR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]!=0),]<0.05)/length(which(x[,9]!=0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))
auc_sim=lapply(results[["FDR"]], function(x) t(data.frame(unlist(lapply(x, function(y) auc(roc(as.numeric(x[,9]!=0),y)))))))

TPR_mat=as.data.frame(matrix(unlist(TPR_sim), ncol=9, byrow=TRUE))
colnames(TPR_mat)=colnames(results[["pvals"]][[1]])
TPR_mat=TPR_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]

auc_mat=as.data.frame(matrix(unlist(auc_sim), ncol=9, byrow=TRUE))
colnames(auc_mat)=colnames(results[["pvals"]][[1]])
auc_mat=auc_mat[, 1:8]

mean_TPR=cbind(data.frame("50%"=colMeans(TPR_mat)), mean_TPR)
mean_TNR=cbind(data.frame("50%"=colMeans(TNR_mat)), mean_TNR)
mean_auc=cbind(data.frame("50%"=colMeans(auc_mat)), mean_auc)


write.csv(mean_TPR, file="figures/SupTable21.csv")
write.csv(mean_TNR, file="figures/SupTable22.csv")
write.csv(mean_auc, file="figures/SupTable23.csv")