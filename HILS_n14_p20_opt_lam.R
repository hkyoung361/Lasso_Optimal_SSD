# This script runs the HILS algorithm in the n=14, p=20 (scenario 2 in section 5) case by optimizing the measure over lambda

root_dir = "PATH TO DIRECTORY CONTAINING CODE AND DESIGN CATALOG"

source(paste0(root_dir,'/Lasso_optimal_SSD_function_library.R'))
### The above path should point to whatever the path is in your local files


library(abind)
library(readr)
# Setting design parameters
n=14
p=20
k_star = 5
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)




# Creating a folder for the output
output_folder = "./n14_p20_k5_allpos_designs_HILS_lam_opt"

if(file.exists(output_folder)){
  print("Output folder already exists, overwritting output")
}else{dir.create(output_folder)}



# Generating heuristic designs



start_design_list_Var_s= replicate(50, matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p))
#Generating Var(s+) designs may take some time 
Var_s_values = rep(NA, 50)
output_designs_Var_s = start_design_list_Var_s
for ( i in 1:50){
  print(i)
  Var_S_plus_opt <- coord_exchange_Var_s(n=n, p=p, eff = sample(c(0.5,0.6,0.7)), UES_const= sample(c(0,0.1)))
  print(Var_S_plus_opt$P_0)
  output_designs_Var_s[,,i]<- Var_S_plus_opt$design
  Var_s_values [i]<-Var_S_plus_opt$P_0
  
}


best_var_s <-output_designs_Var_s[,,sort(Var_s_values, index.return = TRUE)$ix[1:50]]



best_huer = best_var_s

candidate_designs = best_huer



for(i in 1:50){write.csv(best_huer[,,i], file = paste("./n14_p20_k5_allpos_designs_HILS_lam_opt/best_huer_Var_s_design_", paste(i, ".csv", sep=''), sep=""), row.names = FALSE, col.names=FALSE)}



#Read in PEDS and add them to the candidate set

d1<-as.matrix(read_in_design(paste0(root_dir, "/Design_Catalog/n14_p20_k5_allpos/d1.txt")))
candidate_designs<- abind(candidate_designs, d1, along=3)



d2<-as.matrix(read_in_design(paste0(root_dir, "/Design_Catalog/n14_p20_k5_allpos/d2.txt")))
candidate_designs<- abind(candidate_designs, d2, along=3)

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
    
    prob <- Joint_prob_general_optimal_lambda(F_0_cs, V_0_half, B_mag=rep(3,k_star), k_star, sign_vects= sign_all_pos, submodels = A_val)
  }
  prob
  
}

stopCluster(cl)
end_time = Sys.time()
print(end_time - start_time)
print(prob_results)

write.csv(prob_results, file ="./n14_p20_k5_allpos_designs_HILS_lam_opt/HILS_eval_results.csv", row.names =FALSE, col.names = FALSE)



#### Take the top 4 and paralellize the construction from here


print("Index of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1:4])

print("Measure of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$x[1:4])


HILS_opt <-  candidate_designs[,,sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1]]


# Save the design
write.csv(HILS_opt, file="./n14_p20_k5_allpos_designs_HILS_lam_opt/HILS_opt_n14_p20_k5_allpos.csv", row.names = F, col.names=F)


