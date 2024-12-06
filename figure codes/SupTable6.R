### summarize simulation results

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

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

write.csv(DEvariance_table, file="figures/SupTable6.csv")


