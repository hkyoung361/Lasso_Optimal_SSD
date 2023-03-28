# This script produces the plots in the n=14, p=20 (scenario 2 in section 5) case 

root_dir = "PATH TO DIRECTORY CONTAINING CODE AND DESIGN CATALOG"

source(paste0(root_dir,'/Lasso_optimal_SSD_function_library.R'))
### The above path should point to whatever the path is in your local files



library(readr)
library(doParallel)
library(foreach)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(purrr)



# Read in the PED

d1<-as.matrix(read_in_design(paste0(root_dir, "/Design_Catalog/n14_p20_k5_allpos/d1.txt")))

# Read in the next best non-PED selected by HILS, in this case it is one of the Var(s+) designs in the calalog


HILS_opt <- as.matrix(read_csv(paste0(root_dir, "/Design_Catalog/n14_p20_k5_allpos/best_huer_Var_s_design_18.csv")))


# Creating a folder for the output
output_folder = "./n14_p20_k5_allpos_designs_HILS_lam_opt"

if(file.exists(output_folder)){
  print("Output folder already exists, overwritting output")
}else{dir.create(output_folder)}



#Function to permute columns
permute_cols <- function(X){
  X_out= X[,sample(ncol(X))]
  return(X_out)
}



n=14
k=20
k_star=5
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)
log_lam = seq(from = -4.5, to= 2, by = 0.05)
submodel_samples = stack_NBIBD(p=k, k=k_star, s_1=64, num_stacks = 10)
A_val = submodel_samples$A_final
A_64 = submodel_samples$A_NBIBD


lambda = exp(log_lam)

I_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0))
S_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0))
joint_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0))

I_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0))
S_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0))
joint_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0))


I_results_d1$design = factor(I_results_d1$design, levels=c("d1","HILS"))
S_results_d1$design = factor(S_results_d1$design, levels=c("d1","HILS"))
joint_results_d1$design = factor(S_results_d1$design, levels=c("d1","HILS"))


I_results_HILS$design = factor(I_results_HILS$design, levels=c("d1","HILS"))
S_results_HILS$design = factor(S_results_HILS$design, levels=c("d1","HILS"))
joint_results_HILS$design = factor(S_results_HILS$design, levels=c("d1","HILS"))




F_0=d1

V_0_half= diag(lapply(as.data.frame(F_0),fn))
F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
for (i in c(1:length(lambda))){
  print(i)
  lam = lambda[i]
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))

  prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                            B_mag=rep(3,k_star),k=k_star,
                                            sign_vects=sign_all_pos,
                                            log_lambda = log(lam),
                                            submodels = A_val,
                                            log_lambda_strategy="fixed",
                                            sigma = 1, output_option="I")$results)$mean
  prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                            B_mag=rep(3,k_star),k=k_star,
                                            sign_vects=sign_all_pos,
                                            log_lambda = log(lam),
                                            submodels = A_val,
                                            log_lambda_strategy="fixed",
                                            sigma = 1, output_option="S")$results)$mean
  prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                            B_mag=rep(3,k_star),k=k_star,
                                            sign_vects=sign_all_pos,
                                            log_lambda = log(lam),
                                            submodels = A_val,
                                            log_lambda_strategy="fixed",
                                            sigma = 1)$results)$mean

  results_I<- c(design= "d1", lambda = lam, prob_I)
  results_S <- c(design="d1", lambda = lam, prob_S)
  results_joint <- c(design="d1", lambda = lam, prob_joint)

  I_results_d1[i,]=results_I
  S_results_d1[i,]=results_S
  joint_results_d1[i,]= results_joint
}



overall_results_d1 = I_results_d1
overall_results_d1$P_S=as.numeric(S_results_d1$P_S)

overall_results_d1$joint = as.numeric(joint_results_d1$P_joint)
write.csv(overall_results_d1, file ="./n14_p20_k5_allpos_designs_HILS_lam_opt/prob_comparison_d1_log_lam_opt.csv", row.names =FALSE, col.names = FALSE)



F_0=HILS_opt

V_0_half= diag(lapply(as.data.frame(F_0),fn))
F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
for (i in c(1:length(lambda))){
  print(i)
  lam = lambda[i]
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))

  prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                            B_mag=rep(3,k_star),k=k_star,
                                            sign_vects=sign_all_pos,
                                            log_lambda = log(lam),
                                            submodels = A_val,
                                            log_lambda_strategy="fixed",
                                            sigma = 1, output_option="I")$results)$mean
  prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                            B_mag=rep(3,k_star),k=k_star,
                                            sign_vects=sign_all_pos,
                                            log_lambda = log(lam),
                                            submodels = A_val,
                                            log_lambda_strategy="fixed",
                                            sigma = 1, output_option="S")$results)$mean
  prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                                B_mag=rep(3,k_star),k=k_star,
                                                sign_vects=sign_all_pos,
                                                log_lambda = log(lam),
                                                submodels = A_val,
                                                log_lambda_strategy="fixed",
                                                sigma = 1)$results)$mean

  results_I<- c(design= "HILS", lambda = lam, prob_I)
  results_S <- c(design="HILS", lambda = lam, prob_S)
  results_joint <- c(design="HILS", lambda = lam, prob_joint)

  I_results_HILS[i,]=results_I
  S_results_HILS[i,]=results_S
  joint_results_HILS[i,]= results_joint
  
}


overall_results_HILS = I_results_HILS
overall_results_HILS$P_S=as.numeric(S_results_HILS$P_S)
overall_results_HILS$joint = as.numeric(joint_results_HILS$P_joint)

write.csv(overall_results_HILS, file ="./n14_p20_k5_allpos_designs_HILS_lam_opt/prob_comparison_HILS_log_lam_opt.csv", row.names =FALSE, col.names = FALSE)



I_results=rbind(I_results_d1, I_results_HILS)

S_results=rbind(S_results_d1, S_results_HILS)
joint_results = rbind(joint_results_d1, joint_results_HILS)
overall_results = I_results
overall_results$P_S=as.numeric(S_results$P_S)
overall_results$joint=as.numeric(joint_results$P_joint)


write.csv(overall_results, file ="./n14_p20_k5_allpos_designs_HILS_lam_opt/prob_comparison_log_lam_opt.csv", row.names =FALSE, col.names = FALSE)



prob_comparison_d1 <- read_csv("prob_comparison_d1_log_lam_opt.csv")

prob_comparison_HILS <- read_csv("prob_comparison_HILS_log_lam_opt.csv")

overall_results= rbind(prob_comparison_d1, prob_comparison_HILS)

trimmed_overall=overall_results
trimmed_overall[trimmed_overall$design=="d1",]$design="PED"
trimmed_overall[trimmed_overall$design=="HILS",]$design="Var(s+)_Eff80"



# Plotting

ggjoint = ggplot(data = trimmed_overall,  aes(x=log(as.numeric(lambda)), y= as.numeric(joint), color= design, linetype = design))+
  geom_line()+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\phi_l(\\textbf{z}_A=\\textbf{1}_5, \\lambda)"))+
  scale_color_manual(values =c("blue","red"))+
  theme_bw()+
  scale_x_continuous(limits = c(-2,2))+
  scale_linetype_manual(values =c("solid", "dashed"))+
  theme( legend.title=element_blank(), text=element_text(size=16), legend.position = c(0.2,.75))

 ggsave("./Section5_n14_p20_joint_figure.png", dpi=600, device="png", ggjoint) ## this saves a 600 dpi image.
# 
# 
