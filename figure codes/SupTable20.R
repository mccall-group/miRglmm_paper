load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

MSE_mat=as.data.frame(matrix(unlist(results[["MSE_sim"]]), ncol=7, byrow=TRUE))
colnames(MSE_mat)=colnames(results[["MSE_sim"]][[1]])
mean_MSE=data.frame("100%"=colMeans(MSE_mat))


load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults_75perlibsize.rda")

MSE_mat=as.data.frame(matrix(unlist(results[["MSE_sim"]]), ncol=7, byrow=TRUE))
colnames(MSE_mat)=colnames(results[["MSE_sim"]][[1]])

mean_MSE=cbind(data.frame("75%"=colMeans(MSE_mat)), mean_MSE)


load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults_50perlibsize.rda")

MSE_mat=as.data.frame(matrix(unlist(results[["MSE_sim"]]), ncol=7, byrow=TRUE))
colnames(MSE_mat)=colnames(results[["MSE_sim"]][[1]])

mean_MSE=cbind(data.frame("50%"=colMeans(MSE_mat)), mean_MSE)



write.csv(mean_MSE, file="figures/SupTable20.csv")
