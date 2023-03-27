#source('/Users/hkyoung/Desktop/SIC_Paper_Code/SIC_Paper_Function_Library.R')

source('~/SIC_Paper/SIC_Paper_Function_Library.R')
### The above path should point to whatever the path is in your local files


library(abind)
library(readr)
# Setting design parameters
n=14
p=20
k_star = 5
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)





# Creating a folder for the output
output_folder = "./n14_k20_kstar5_allpos_designs_HILS"

if(file.exists(output_folder)){
  print("Output folder already exists, overwritting output")
}else{dir.create(output_folder)}



candidate_designs <- as.matrix(read_csv("./n14_k20_kstar5_allpos_designs_HILS_lam_opt/best_huer_Var_s_design_1.csv"))


for(i in 2:50){
  H_1 = as.matrix(read_csv(paste("./n14_k20_kstar5_allpos_designs_HILS_lam_opt/best_huer_Var_s_design_", paste(i, ".csv", sep=''), sep="")))
  candidate_designs= abind(candidate_designs, H_1, along =3)
}


# d1<-as.matrix(read_in_design("./n14_k20_kstar5_allpos_designs_HILS/d1.txt"))
# candidate_designs<- abind(candidate_designs, d1, along=3)
# 
# 
# 
# d2<-as.matrix(read_in_design("./n14_k20_kstar5_allpos_designs_HILS/d2.txt"))
# candidate_designs<- abind(candidate_designs, d2, along=3)

sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)

submodel_samples = stack_NBIBD(p=p, k=k_star, s_1=64, num_stacks = 15)

A_val = submodel_samples$A_final
A_64 = submodel_samples$A_NBIBD


# only using All positive signs here
start_time <- Sys.time()
cl <- makeCluster(8)
registerDoParallel(cl)



prob_results = foreach(start=1:dim(candidate_designs)[3], .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"),.export = ls(globalenv())) %dopar%{
  F_0 = candidate_designs[,,start]
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  if(any(is.na(F_0_cs))){
    prob=0}
  else{
    
    prob <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half, 
                                            B_mag=rep(3,k_star),k=k_star,
                                            sign_vects=sign_all_pos,
                                            log_lambda = 1,
                                            submodels = A_val,
                                            log_lambda_strategy="integrate",
                                            sigma = 1)$results)$mean
  }
  prob
  
}

stopCluster(cl)
end_time = Sys.time()
print(end_time - start_time)

write.csv(prob_results, file ="./n14_k20_kstar5_allpos_designs_HILS_lam_opt/HILS_eval_results_int.csv", row.names =FALSE, col.names = FALSE)


print(prob_results)

#### Take the top 4 and paralellize the construction from here


print("Index of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1:4])

print("Measure of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$x[1:4])

#starts_for_pairwise = best_huer[,,sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1:4]]

HILS_opt <-  candidate_designs[,,sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1]]


# Save the design
write.csv(HILS_opt, file="./n14_k20_kstar5_allpos_designs_HILS_lam_opt/HILS_opt_n14_k20_kstar5_allpos_integral.csv", row.names = F, col.names=F)


