### summarize simulation results

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

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

write.csv(nullvariance_table, file="figures/SupTable5.csv")


