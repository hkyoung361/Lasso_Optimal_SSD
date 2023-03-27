#source('/Users/hkyoung/Desktop/SIC_Paper_Code/SIC_Paper_Function_Library.R')


root_dir = "PATH TO DIRECTORY CONTAINING CODE"

source(paste0(root_dir,'/SIC_Paper_Function_Library.R'))
### The above path should point to whatever the path is in your local files


library(abind)
library(readr)
# Setting design parameters
n=14
p=20
k_star = 5
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)


perturb_design<-function(D_mat,num_elts){
  n = dim(D_mat)[1]
  p = dim(D_mat)[2]
  
  for (elt in c(1:num_elts)){
    i = sample(seq(1:n), size =1)
    j = sample(seq(1:p), size =1)
    D_mat[i,j] = -D_mat[i,j]
  }
  return(D_mat)
}


# Creating a folder for the output
output_folder = "./n14_k20_kstar5_allpos_designs_HILS"

if(file.exists(output_folder)){
  print("Output folder already exists, overwritting output")
}else{dir.create(output_folder)}



# Generating heuristic designs


# start_design_list_UES_sq = replicate(25, matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p))
# output_designs_UES_sq = start_design_list_UES_sq
# UES_sq_values = rep(NA, 25)
# for ( i in 1:25){
#   UES_squared_optimal <- coord_exchange_UES_sq(start_design_list_UES_sq[,,i], max_iter=10000)
#   output_designs_UES_sq[,,i] <-UES_squared_optimal$design
#   UES_sq_values[i] = UES_squared_optimal$P_0
# }

start_design_list_Var_s= replicate(50, matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p))
#Generating Var(s+) designs may take some time given Byran's proposed changes have been implemented
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

# best_huer <-abind(best_UES_sq, best_var_s, along = 3)

best_huer = best_var_s

candidate_designs = best_huer

# save these designs in the output folder
# for (i in 1:25){ write.csv(best_huer[,,i], file = paste("./n9_k10_kstar3_allpos_designs_HILS/best_huer_UES_sq_design_", paste(i, ".csv", sep=''), sep=""), row.names = FALSE, col.names=FALSE)}

for(i in 1:50){write.csv(best_huer[,,i], file = paste("./n14_k20_kstar5_allpos_designs_HILS/best_huer_Var_s_design_", paste(i, ".csv", sep=''), sep=""), row.names = FALSE, col.names=FALSE)}

# for (i in 1:50){
#   perturb_design_1 <- perturb_design(best_huer[,,i],1)
#   candidate_designs<- abind(candidate_designs, perturb_design_1, along=3)
#   perturb_design_2<-perturb_design(best_huer[,,i],2)
#   candidate_designs<- abind(candidate_designs, perturb_design_2, along=3)
#   # perturb_design_3<-perturb_design(best_huer[,,i],3)
#   # candidate_designs<- abind(candidate_designs, perturb_design_3, along=3)
#   write.csv(perturb_design_1, file = paste("./n14_k20_kstar5_allpos_designs_HILS/best_huer_Var_s_design_", paste(i, "perturb_1.csv", sep=''), sep=""), row.names = FALSE, col.names=FALSE)
#   write.csv(perturb_design_2, file = paste("./n14_k20_kstar5_allpos_designs_HILS/best_huer_Var_s_design_", paste(i, "perturb_2.csv", sep=''), sep=""), row.names = FALSE, col.names=FALSE)
#   #write.csv(perturb_design_3, file = paste("./n14_k20_kstar5_allpos_designs_HILS/best_huer_Var_s_design_", paste(i, "perturb_3.csv", sep=''), sep=""), row.names = FALSE, col.names=FALSE)
# }

#Read in Rahki's design and add them to the candidate set

d1<-as.matrix(read_in_design("./n14_k20_kstar5_allpos_designs_HILS/d1.txt"))
candidate_designs<- abind(candidate_designs, d1, along=3)



d2<-as.matrix(read_in_design("./n14_k20_kstar5_allpos_designs_HILS/d2.txt"))
candidate_designs<- abind(candidate_designs, d2, along=3)

sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)

submodel_samples = stack_NBIBD(p=p, k=k_star, s_1=64, num_stacks = 15)

A_val = submodel_samples$A_final
A_64 = submodel_samples$A_NBIBD


# only using All positive signs here
start_time <- Sys.time()
cl <- makeCluster(8)
registerDoParallel(cl)
F_0 = candidate_designs[,,1]
V_0_half= diag(lapply(as.data.frame(F_0),fn))
F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))

prob <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half, 
                                          B_mag=rep(3,k_star),k=k_star,
                                          sign_vects=sign_all_pos,
                                          log_lambda = 1,
                                          submodels = A_val,
                                          log_lambda_strategy="integrate",
                                          sigma = 1)$results)$mean
prob




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

write.csv(prob_results, file ="./n14_k20_kstar5_allpos_designs_HILS/HILS_eval_results.csv", row.names =FALSE, col.names = FALSE)


print(prob_results)

#### Take the top 4 and paralellize the construction from here


print("Index of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1:4])

print("Measure of top 4 designs")

print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$x[1:4])

#starts_for_pairwise = best_huer[,,sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1:4]]

HILS_opt <-  candidate_designs[,,sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1]]


# Save the design
write.csv(HILS_opt, file="./n14_k20_kstar5_allpos_designs_HILS/HILS_opt_n14_k20_kstar5_allpos.csv", row.names = F, col.names=F)


