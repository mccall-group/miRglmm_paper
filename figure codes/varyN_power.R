load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))

TPR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]!=0),]<0.05)/length(which(x[,9]!=0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))

TPR_mat=as.data.frame(matrix(unlist(TPR_sim), ncol=9, byrow=TRUE))
colnames(TPR_mat)=colnames(results[["pvals"]][[1]])
TPR_mat=TPR_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]


mean_TPR=data.frame("N39"=colMeans(TPR_mat))
mean_TNR=data.frame("N39"=colMeans(TNR_mat))

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults_N15.rda")

results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))

TPR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]!=0),]<0.05)/length(which(x[,9]!=0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))

TPR_mat=as.data.frame(matrix(unlist(TPR_sim), ncol=9, byrow=TRUE))
colnames(TPR_mat)=colnames(results[["pvals"]][[1]])
TPR_mat=TPR_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]

mean_TPR=cbind(data.frame("N30"=colMeans(TPR_mat)), mean_TPR)
mean_TNR=cbind(data.frame("N30"=colMeans(TNR_mat)), mean_TNR)

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults_N10.rda")

results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))

TPR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]!=0),]<0.05)/length(which(x[,9]!=0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))

TPR_mat=as.data.frame(matrix(unlist(TPR_sim), ncol=9, byrow=TRUE))
colnames(TPR_mat)=colnames(results[["pvals"]][[1]])
TPR_mat=TPR_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]

mean_TPR=cbind(data.frame("N20"=colMeans(TPR_mat)), mean_TPR)
mean_TNR=cbind(data.frame("N20"=colMeans(TNR_mat)), mean_TNR)

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults_N5.rda")

results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))

TPR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]!=0),]<0.05)/length(which(x[,9]!=0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))

TPR_mat=as.data.frame(matrix(unlist(TPR_sim), ncol=9, byrow=TRUE))
colnames(TPR_mat)=colnames(results[["pvals"]][[1]])
TPR_mat=TPR_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]

mean_TPR=cbind(data.frame("N10"=colMeans(TPR_mat)), mean_TPR)
mean_TNR=cbind(data.frame("N10"=colMeans(TNR_mat)), mean_TNR)