load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

MSE_mat=as.data.frame(matrix(unlist(results[["MSE_sim"]]), ncol=7, byrow=TRUE))
colnames(MSE_mat)=colnames(results[["MSE_sim"]][[1]])
mean_MSE=data.frame("FC2"=colMeans(MSE_mat))

load(file="sim results/sims_N100_m15_s1_rtruncnorm12_combinedresults.rda")

MSE_mat=as.data.frame(matrix(unlist(results[["MSE_sim"]]), ncol=7, byrow=TRUE))
colnames(MSE_mat)=colnames(results[["MSE_sim"]][[1]])

mean_MSE=cbind(data.frame("FC1.5"=colMeans(MSE_mat)), mean_MSE)

load(file="sim results/sims_N100_m4_s1_rtruncnorm35_combinedresults.rda")

MSE_mat=as.data.frame(matrix(unlist(results[["MSE_sim"]]), ncol=7, byrow=TRUE))
colnames(MSE_mat)=colnames(results[["MSE_sim"]][[1]])

mean_MSE=cbind(mean_MSE, data.frame("FC4"=colMeans(MSE_mat)))



write.csv(mean_MSE, file="figures/SupTable16.csv")
