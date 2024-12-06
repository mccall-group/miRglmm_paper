### summarize simulation results

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

mean_TPR=colMeans(TPR_mat)
std_TPR=apply(TPR_mat,2 ,sd)
median_TPR=apply(TPR_mat, 2, median)
min_TPR=apply(TPR_mat,2,min)
max_TPR=apply(TPR_mat,2,max)

TPR_table=data.frame("mean TPR"=mean_TPR, "SD TPR"=std_TPR, "median TPR"=median_TPR, "minimum TPR"=min_TPR, "maximum TPR"=max_TPR)

mean_TNR=colMeans(TNR_mat)
std_TNR=apply(TNR_mat,2 ,sd)
median_TNR=apply(TNR_mat, 2, median)
min_TNR=apply(TNR_mat,2,min)
max_TNR=apply(TNR_mat,2,max)

TNR_table=data.frame("mean TNR"=mean_TNR, "SD TNR"=std_TNR, "median TNR"=median_TNR, "minimum TNR"=min_TNR, "maximum TNR"=max_TNR)

mean_auc=colMeans(auc_mat)
std_auc=apply(auc_mat,2 ,sd)
median_auc=apply(auc_mat, 2, median)
min_auc=apply(auc_mat,2,min)
max_auc=apply(auc_mat,2,max)

auc_table=data.frame("mean AUC"=mean_auc, "SD AUC"=std_auc, "median AUC"=median_auc, "minimum AUC"=min_auc, "maximum AUC"=max_auc)


write.csv(TPR_table, file="figures/SupTable7.csv")
write.csv(TNR_table, file="figures/SupTable8.csv")
write.csv(auc_table, file="figures/SupTable9.csv")



