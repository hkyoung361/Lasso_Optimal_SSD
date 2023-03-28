# This script reranks the designs used in the HILS algorithm for the n=9, p=10 (scenario 1 in section 5) case but optimizes over lambda

root_dir = "PATH TO DIRECTORY CONTAINING CODE AND DESIGN CATALOG"

source(paste0(root_dir,'/Lasso_optimal_SSD_function_library.R'))
### The above path should point to whatever the path is in your local files


library(abind)
library(readr)
# Setting design parameters
n=9
p=10
k_star = 3



# Creating a folder for the output
output_folder = "./n9_p10_k3_allsigns_designs_HILS"

if(file.exists(output_folder)){
  print("Output folder already exists, overwritting output")
}else{dir.create(output_folder)}



candidate_designs <- as.matrix(read_csv("./n9_p10_k3_allsigns_designs_HILS/best_huer_UES_sq_design_1.csv"))

for(i in 2:100){
  if (i<= 50){
    H_1 = as.matrix(read_csv(paste("./n9_p10_k3_allsigns_designs_HILS/best_huer_UES_sq_design_", paste(i, ".csv", sep=''), sep="")))
  }
  else{
    H_1 = as.matrix(read_csv(paste("./n9_p10_k3_allsigns_designs_HILS/best_huer_Var_s_design_", paste(i-50, ".csv", sep=''), sep="")))
  }
  
  candidate_designs= abind(candidate_designs, H_1, along =3)
}





all_signs <- gen_all_sign_vects(seq(1,k_star,1))


# Using all signs here
start_time <- Sys.time()
#parallelize
cl <- makeCluster(8)
registerDoParallel(cl)
prob_results = foreach(start=1:dim(candidate_designs)[3], .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"),.export = ls(globalenv())) %dopar%{
  F_0 = candidate_designs[,,start]
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  if(any(is.na(F_0_cs))){
    prob=0}
  else{
   
    prob <- Joint_prob_general_optimal_lambda(F_0_cs, V_0_half, B_mag=rep(3,k_star), k_star, sign_vects= all_signs, submodels = NULL)
  }
  prob
  
}

stopCluster(cl)
end_time = Sys.time()
print(end_time - start_time)

# write the measure results for the heuristic designs to a file. 

write.csv(prob_results, file ="./n9_p10_k3_allsigns_designs_HILS/HILS_eval_results_opt_lam.csv", row.names =FALSE, col.names = FALSE)


print(prob_results)

#### Take the top 4 and paralellize the construction from here


print("Index of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1:4])

print("Measure of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$x[1:4])



HILS_opt <-  best_huer[,,sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1]]


# Save the design
write.csv(HILS_opt, file="./n9_p10_k3_allsigns_designs_HILS/HILS_opt_n9_p10_k3_allsigns_opt_lam.csv", row.names = F, col.names=F)






