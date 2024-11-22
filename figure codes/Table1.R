### summarize simulation results

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

MSE_mat=as.data.frame(matrix(unlist(results[["MSE_sim"]]), ncol=7, byrow=TRUE))
colnames(MSE_mat)=colnames(results[["MSE_sim"]][[1]])

mean_MSE=colMeans(MSE_mat)
std_MSE=apply(MSE_mat,2 ,sd)
median_MSE=apply(MSE_mat, 2, median)
min_MSE=apply(MSE_mat,2,min)
max_MSE=apply(MSE_mat,2,max)
n_min_MSE=apply(MSE_mat,1, which.min)
n_min_MSE=factor(n_min_MSE, levels=c("1", "2", "3", "4", "5", "6", "7"))
out=data.frame(table(n_min_MSE))

MSE_table=data.frame("mean MSE"=mean_MSE*1000, "SD MSE"=std_MSE*1000, "median MSE"=median_MSE*1000, "minimum MSE"=min_MSE*1000, "maximum MSE"=max_MSE*1000, "number of times minimizing MSE"=out[,2])


cov_mat=as.data.frame(matrix(unlist(results[["coverage_prob_sim"]]), ncol=6, byrow=TRUE))
colnames(cov_mat)=colnames(results[["coverage_prob_sim"]][[1]])

mean_cov=colMeans(cov_mat)
std_cov=apply(cov_mat,2 ,sd)
median_cov=apply(cov_mat, 2, median)
min_cov=apply(cov_mat,2,min)
max_cov=apply(cov_mat,2,max)
cov_table=data.frame("mean cov"=mean_cov, "SD cov"=std_cov, "median cov"=median_cov, "minimum cov"=min_cov, "maximum cov"=max_cov)


nullvar_sim=lapply(results[["beta_hat"]], function(x) t(data.frame(variance=apply(x[which(x$true_beta==0),],2, var))))
nullvar_mat=as.data.frame(matrix(unlist(nullvar_sim), ncol=8, byrow=TRUE))
colnames(nullvar_mat)=colnames(nullvar_sim[[1]])
nullvar_mat=nullvar_mat[, 1:7]

mean_nullvar=colMeans(nullvar_mat)
std_nullvar=apply(nullvar_mat,2 ,sd)
median_nullvar=apply(nullvar_mat, 2, median)
min_nullvar=apply(nullvar_mat,2,min)
max_nullvar=apply(nullvar_mat,2,max)

nullvariance_table=data.frame("mean null var"=mean_nullvar, "SD null var"=std_nullvar, 
                          "median null var"=median_nullvar, "minimum null var"=min_nullvar, 
                          "maximum null var"=max_nullvar)

DEvar_sim=lapply(results[["beta_hat"]], function(x) t(data.frame(variance=apply(rbind(-1*x[which(x$true_beta<0),],x[which(x$true_beta>0),]),2, var))))
DEvar_mat=as.data.frame(matrix(unlist(DEvar_sim), ncol=8, byrow=TRUE))
colnames(DEvar_mat)=colnames(DEvar_sim[[1]])
DEvar_mat=DEvar_mat[, 1:7]

mean_DEvar=colMeans(DEvar_mat)
std_DEvar=apply(DEvar_mat,2 ,sd)
median_DEvar=apply(DEvar_mat, 2, median)
min_DEvar=apply(DEvar_mat,2,min)
max_DEvar=apply(DEvar_mat,2,max)

DEvariance_table=data.frame("mean DE var"=mean_DEvar, "SD DE var"=std_DEvar, 
                              "median DE var"=median_DEvar, "minimum DE var"=min_DEvar, 
                              "maximum DE var"=max_DEvar)


# avgCIwidth_sim=lapply(results[["CI_width"]], function(x) t(data.frame(precision=colMeans(x))))
# precision_mat=as.data.frame(matrix(unlist(avgCIwidth_sim), ncol=6, byrow=TRUE))
# colnames(precision_mat)=colnames(avgCIwidth_sim[[1]])
# 
# 
# mean_precision=colMeans(precision_mat)
# std_precision=apply(precision_mat,2 ,sd)
# median_precision=apply(precision_mat, 2, median)
# min_precision=apply(precision_mat,2,min)
# max_precision=apply(precision_mat,2,max)
# 
# precision_table=data.frame("mean precision"=mean_precision, "SD precision"=std_precision, "median precision"=median_precision, "minimum precision"=min_precision, "maximum precision"=max_precision)
results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))

TPR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]!=0),]<0.05)/length(which(x[,9]!=0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))

TPR_mat=as.data.frame(matrix(unlist(TPR_sim), ncol=9, byrow=TRUE))
colnames(TPR_mat)=colnames(results[["pvals"]][[1]])
TPR_mat=TPR_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]


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

table_out=data.frame(cbind("TPR"=TPR_table[,1],  "coverage proportion"=cov_table[match(rownames(TPR_table), rownames(cov_table)),1]))
rownames(table_out)=rownames(TPR_table)
table_out=cbind(table_out,  "null variance"=nullvariance_table[match(rownames(table_out), rownames(nullvariance_table)),1]*1000)
table_out=cbind(table_out,  "DE variance"=DEvariance_table[match(rownames(table_out), rownames(DEvariance_table)),1]*1000)
table_out=cbind(table_out,  "MSE"=MSE_table[match(rownames(table_out), rownames(MSE_table)),1])
table_out=cbind(table_out,  "TNR"=TNR_table[match(rownames(table_out), rownames(TNR_table)),1])
table_out=table_out[, c("MSE", "coverage.proportion", "null variance", "DE variance", "TPR", "TNR")]
write.csv(table_out, file="figures/resub/Table1.csv")

write.csv(MSE_table, file="figures/resub/Table1_MSE.csv")
write.csv(cov_table, file="figures/resub/Table1_cov.csv")
# write.csv(nullvariance_table, file="figures/Table1_nullvariance.csv")
# write.csv(DEvariance_table, file="figures/Table1_DEvariance.csv")
# # write.csv(precision_table, file="figures/Table1_precision.csv")
# write.csv(TPR_table, file="figures/Table1_TPR.csv")
# write.csv(TNR_table, file="figures/Table1_TNR.csv")





mean(MSE_mat$miRglmm[which(n_min_MSE==1)]-MSE_mat$DESeq2[which(n_min_MSE==1)])
mean(MSE_mat$DESeq2[which(n_min_MSE==3)]-MSE_mat$miRglmm[which(n_min_MSE==3)])
median(MSE_mat$miRglmm[which(n_min_MSE==1)]-MSE_mat$DESeq2[which(n_min_MSE==1)])
median(MSE_mat$DESeq2[which(n_min_MSE==3)]-MSE_mat$miRglmm[which(n_min_MSE==3)])