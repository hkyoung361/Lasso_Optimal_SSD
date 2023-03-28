# This script runs the HILS algorithm in the n=9, p=10 (scenario 1 in section 5) case by integrating the measure over lambda

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



# Generating heuristic designs


start_design_list_UES_sq = replicate(50, matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p))
output_designs_UES_sq = start_design_list_UES_sq
UES_sq_values = rep(NA, 50)
for ( i in 1:50){
  UES_squared_optimal <- coord_exchange_UES_sq(start_design_list_UES_sq[,,i], max_iter=10000)
  output_designs_UES_sq[,,i] <-UES_squared_optimal$design
  UES_sq_values[i] = UES_squared_optimal$P_0
}

start_design_list_Var_s= replicate(50, matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p))

#Generating Var(s+) designs may take some time given Byran's proposed changes have been implemented
Var_s_values = rep(NA, 50)
output_designs_Var_s = start_design_list_Var_s
for ( i in 1:50){
  print(i)
  Var_S_plus_opt <- coord_exchange_Var_s(n=n, p=p, eff = sample(c(0.5,0.6,0.7)), UES_const= 0)
  output_designs_Var_s[,,i]<- Var_S_plus_opt$design
  Var_s_values [i]<-Var_S_plus_opt$P_0
  
}




best_UES_sq <-output_designs_UES_sq[,,sort(UES_sq_values, index.return = TRUE)$ix[1:50]]

best_var_s <-output_designs_Var_s[,,sort(Var_s_values, index.return = TRUE)$ix[1:50]]

best_huer <-abind(best_UES_sq, best_var_s, along = 3)

#best_huer = best_var_s

candidate_designs = best_huer

# save these designs in the output folder
for (i in 1:50){ write.csv(best_huer[,,i], file = paste("./n9_p10_k3_allsigns_designs_HILS/best_huer_UES_sq_design_", paste(i, ".csv", sep=''), sep=""), row.names = FALSE, col.names=FALSE)}

for(i in 51:100){write.csv(best_huer[,,i], file = paste("./n9_p10_k3_allsigns_designs_HILS/best_huer_Var_s_design_", paste(i-50, ".csv", sep=''), sep=""), row.names = FALSE, col.names=FALSE)}






#Read in PEDs and add them to the candidate set

d1<-as.matrix(read_in_design(paste0(root_dir, "/Design_Catalog/n9_p10_k3_allsigns/d1.txt")))
candidate_designs<- abind(candidate_designs, d1, along=3)



d2<-as.matrix(read_in_design(paste0(root_dir, "/Design_Catalog/n9_p10_k3_allsigns/d1.txt")))
candidate_designs<- abind(candidate_designs, d2, along=3)
#### Get the lasso recovery measure of these  designs



sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)
all_signs <- gen_all_sign_vects(seq(1,k_star,1))

# only using All positive signs here
start_time <- Sys.time()
#parallelizing
cl <- makeCluster(16)
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
                                          sign_vects=all_signs,
                                          log_lambda = 1,
                                          submodels = NULL,
                                          log_lambda_strategy="integrate",
                                          sigma = 1)$results)$mean
  }
  prob
  
}

stopCluster(cl)
end_time = Sys.time()
print(end_time - start_time)

# write the measure results for the heuristic designs to a file. 

write.csv(prob_results, file ="./n9_p10_k3_allsigns_designs_HILS/HILS_eval_results.csv", row.names =FALSE, col.names = FALSE)


print(prob_results)

#### Take the top 4 and paralellize the construction from here


print("Index of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1:4])

print("Measure of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$x[1:4])

#starts_for_pairwise = best_huer[,,sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1:4]]

HILS_opt <-  best_huer[,,sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1]]


# Save the design
write.csv(HILS_opt, file="./n9_p10_k3_allsigns_designs_HILS/HILS_opt_n9_p10_k3_allsigns.csv", row.names = F, col.names=F)






