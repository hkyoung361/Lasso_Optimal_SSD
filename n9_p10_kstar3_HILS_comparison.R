
# This script produces the plots in the n=9, p=20 (scenario 2 in section 5) case 

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

d1<-as.matrix(read_in_design(paste0(root_dir, "/Design_Catalog/n9_p10_k3_allsigns/d1.txt")))

# Read in the HILS optimal design selected by integrating over lambda

HILS_opt <- as.matrix(read_csv(paste0(root_dir, "/Design_Catalog/n9_p10_k3_allsigns/HILS_opt_n9_k10_kstar3_allsigns_int.csv")))


# Read in the HILS optimal design selected by optimizing over lambda

HILS_lam_opt <- as.matrix(read_csv(paste0(root_dir, "/Design_Catalog/n9_p10_k3_allsigns/HILS_opt_n9_k10_kstar3_allsigns_lam_opt.csv")))


# Creating a folder for the output
output_folder = "./n9_p10_k3_allsigns_designs_HILS"

if(file.exists(output_folder)){
  print("Output folder already exists, overwritting output")
}else{dir.create(output_folder)}




n=9
k=10
k_star=3
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)
all_signs = gen_all_sign_vects(c(1,2,3))

log_lam = seq(-4.5, 2, by = 0.05)
lambda = exp(log_lam)

I_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0), sign=character(0))
S_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0),sign=character(0))
joint_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0),sign=character(0))

I_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0),sign=character(0))
S_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0),sign=character(0))
joint_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0),sign=character(0))
I_results_HILS_2 = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0),sign=character(0))
S_results_HILS_2 = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0),sign=character(0))
joint_results_HILS_2 = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0),sign=character(0))

I_results_d1$design = factor(I_results_d1$design, levels=c("d1","HILS"))
S_results_d1$design = factor(S_results_d1$design, levels=c("d1","HILS"))
joint_results_d1$design = factor(S_results_d1$design, levels=c("d1","HILS"))
I_results_d1$sign = factor(I_results_d1$sign, levels=c("s_1","s_2","s_3","s_4"))
S_results_d1$sign = factor(S_results_d1$sign, levels=c("s_1","s_2","s_3","s_4"))
joint_results_d1$sign = factor(S_results_d1$sign, levels=c("s_1","s_2","s_3","s_4"))

I_results_HILS$design = factor(I_results_HILS$design, levels=c("d1","HILS"))
S_results_HILS$design = factor(S_results_HILS$design, levels=c("d1","HILS"))
joint_results_HILS$design = factor(S_results_HILS$design, levels=c("d1","HILS"))
I_results_HILS$sign = factor(I_results_HILS$sign, levels=c("s_1","s_2","s_3","s_4"))
S_results_HILS$sign = factor(S_results_HILS$sign, levels=c("s_1","s_2","s_3","s_4"))
joint_results_HILS$sign = factor(S_results_HILS$sign, levels=c("s_1","s_2","s_3","s_4"))


I_results_HILS_2$design = factor(I_results_HILS_2$design, levels=c("d1","HILS","HILS2"))
S_results_HILS_2$design = factor(S_results_HILS_2$design, levels=c("d1","HILS", "HILS2"))
joint_results_HILS_2$design = factor(S_results_HILS_2$design, levels=c("d1","HILS", "HILS2"))
I_results_HILS_2$sign = factor(I_results_HILS_2$sign, levels=c("s_1","s_2","s_3","s_4"))
S_results_HILS_2$sign = factor(S_results_HILS_2$sign, levels=c("s_1","s_2","s_3","s_4"))
joint_results_HILS_2$sign = factor(S_results_HILS_2$sign, levels=c("s_1","s_2","s_3","s_4"))


for (s in c(1:4)){
  sign_vec =all_signs[s,]
  print(sign_vec)
  print(paste0("s_",s))
  F_0=d1
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  for (i in c(1:length(log_lam))){
    print(i)
    log_lambda = log_lam[i]
    V_0_half= diag(lapply(as.data.frame(F_0),fn))
    F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))

    prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                              B_mag=rep(3,k_star),k=k_star,
                                              sign_vects=all_signs[s,],
                                              log_lambda = log_lambda,
                                              submodels = NULL,
                                              log_lambda_strategy="fixed",
                                              sigma = 1, output_option="I")$results)$mean
    prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                              B_mag=rep(3,k_star),k=k_star,
                                              sign_vects=all_signs[s,],
                                              log_lambda = log_lambda,
                                              submodels = NULL,
                                              log_lambda_strategy="fixed",
                                              sigma = 1, output_option="S")$results)$mean
    prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                               B_mag=rep(3,k_star),k=k_star,
                                               sign_vects=all_signs[s,],
                                               log_lambda = log_lambda,
                                               submodels = NULL,
                                               log_lambda_strategy="fixed",
                                               sigma = 1)$results)$mean

    results_I<- c(design= "d1", lambda = exp(log_lambda), prob_I, paste0("s_",s))
    results_S <- c(design="d1", lambda = exp(log_lambda), prob_S, paste0("s_",s))
    results_joint <- c(design="d1", lambda = exp(log_lambda), prob_joint, paste0("s_",s))
    index = ((s-1)*length(log_lam))+i
    print(paste0("Index=",index))

    I_results_d1[index,]=results_I
    S_results_d1[index,]=results_S
    joint_results_d1[index,]= results_joint
  }
  F_0=HILS_opt

  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  for (i in c(1:length(log_lam))){
    print(i)
    #lam = lambda[i]
    log_lambda = log_lam[i]
    V_0_half= diag(lapply(as.data.frame(F_0),fn))
    F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))

   

    prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                              B_mag=rep(3,k_star),k=k_star,
                                              sign_vects=all_signs[s,],
                                              log_lambda = log_lambda,
                                              submodels = NULL,
                                              log_lambda_strategy="fixed",
                                              sigma = 1, output_option="I")$results)$mean
    prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                              B_mag=rep(3,k_star),k=k_star,
                                              sign_vects=all_signs[s,],
                                              log_lambda = log_lambda,
                                              submodels = NULL,
                                              log_lambda_strategy="fixed",
                                              sigma = 1, output_option="S")$results)$mean
    prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                                  B_mag=rep(3,k_star),k=k_star,
                                                  sign_vects=all_signs[s,],
                                                  log_lambda = log_lambda,
                                                  submodels = NULL,
                                                  log_lambda_strategy="fixed",
                                                  sigma = 1)$results)$mean


    results_I<- c(design= "HILS", lambda = exp(log_lambda), prob_I, paste0("s_",s))
    results_S <- c(design="HILS", lambda = exp(log_lambda), prob_S, paste0("s_",s))
    results_joint <- c(design="HILS", lambda = exp(log_lambda), prob_joint, paste0("s_",s))
    index = ((s-1)*length(log_lam))+i
    print(paste0("Index=",index))

    I_results_HILS[index,]=results_I
    S_results_HILS[index,]=results_S
    joint_results_HILS[index,]= results_joint
  }

}


for (s in c(1:4)){
  sign_vec =all_signs[s,]
  print(sign_vec)
  print(paste0("s_",s))
  F_0=HILS_lam_opt
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  for (i in c(1:length(log_lam))){
    print(i)
    log_lambda = log_lam[i]
    V_0_half= diag(lapply(as.data.frame(F_0),fn))
    F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
    
    prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                              B_mag=rep(3,k_star),k=k_star,
                                              sign_vects=all_signs[s,],
                                              log_lambda = log_lambda,
                                              submodels = NULL,
                                              log_lambda_strategy="fixed",
                                              sigma = 1, output_option="I")$results)$mean
    prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                              B_mag=rep(3,k_star),k=k_star,
                                              sign_vects=all_signs[s,],
                                              log_lambda = log_lambda,
                                              submodels = NULL,
                                              log_lambda_strategy="fixed",
                                              sigma = 1, output_option="S")$results)$mean
    prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                                  B_mag=rep(3,k_star),k=k_star,
                                                  sign_vects=all_signs[s,],
                                                  log_lambda = log_lambda,
                                                  submodels = NULL,
                                                  log_lambda_strategy="fixed",
                                                  sigma = 1)$results)$mean
    
    results_I<- c(design= "HILS2", lambda = exp(log_lambda), prob_I, paste0("s_",s))
    results_S <- c(design="HILS2", lambda = exp(log_lambda), prob_S, paste0("s_",s))
    results_joint <- c(design="HILS2", lambda = exp(log_lambda), prob_joint, paste0("s_",s))
    index = ((s-1)*length(log_lam))+i
    print(paste0("Index=",index))
    
    I_results_HILS_2[index,]=results_I
    S_results_HILS_2[index,]=results_S
    joint_results_HILS_2[index,]= results_joint
  }
 
  
}











I_results=rbind(I_results_d1, I_results_HILS)
# #
S_results=rbind(S_results_d1, S_results_HILS)
# #
joint_results = rbind(joint_results_d1, joint_results_HILS)
#
overall_results = I_results
overall_results$P_S=as.numeric(S_results$P_S)
overall_results$joint=as.numeric(joint_results$P_joint)
#
#
write.csv(overall_results, file ="./n9_k10_k3_allsigns_designs_HILS/prob_comparison_n9_by_sign.csv", row.names =FALSE, col.names = FALSE)


I_results = I_results_HILS_2
S_results = S_results_HILS_2
joint_results = joint_results_HILS_2
overall_results = I_results
overall_results$P_S=as.numeric(S_results$P_S)
overall_results$joint=as.numeric(joint_results$P_joint)
write.csv(overall_results, file ="./n9_p10_k3_allsigns_designs_HILS/prob_comparison_lam_opt_n9_by_sign.csv", row.names =FALSE, col.names = FALSE)


print("file written successfully")


library(readr)


overall_results_by_sign1<- read_csv("prob_comparison_n9_by_sign.csv")
overall_results_by_sign2<- read_csv("prob_comparison_lam_opt_n9_by_sign.csv")
overall_results_by_sign2$design="HILS_max_lam"

overall_results_by_sign= rbind(overall_results_by_sign1, overall_results_by_sign2)

overall_results_by_sign[overall_results_by_sign$design=="d1",]$design="PED"





overall_results_mixed_sign= overall_results_by_sign[overall_results_by_sign$sign!="s_4",]
overall_results_mixed_sign_mean= overall_results_mixed_sign%>%group_by(design,lambda)%>%summarise(across(joint, ~ mean(.x, na.rm = TRUE)), .groups="keep")

overall_results=overall_results_by_sign%>%group_by(design,lambda)%>%summarise(across(joint, ~ mean(.x, na.rm = TRUE)), .groups="keep")

overall_results[overall_results$design=="d1",]$design="PED"

overall_results = overall_results[overall_results$design!="d2",]

# All signs graph, HILS_max_lam is omitted because it overlaps with the other HILS (integrated) design
ggjoint = ggplot(data = overall_results[overall_results$design!="HILS_max_lam",], aes(x=log(as.numeric(lambda)), y= as.numeric(joint), color= design, linetype=design))+
  geom_line()+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  scale_color_manual(values =c("red","blue"), labels = c(latex2exp::TeX("$HILS_{\\Lambda} / HILS_{max}$"), "PED"))+
  scale_x_continuous(limits = c(-2,2))+
  scale_linetype_manual(values =c("dashed", "solid"), labels = c(latex2exp::TeX("$HILS_{\\Lambda} / HILS_{max}$"), "PED"))+
  theme_bw()+
  theme(axis.title.y=element_text( vjust=0.55), legend.title = element_blank(), text=element_text(size=16), legend.position = c(0.2,.75))


# All positive sign graph
ggjoint_s4 = ggplot(data = overall_results_by_sign[overall_results_by_sign$sign=="s_4" & overall_results_by_sign$lambda >= .12,], aes(x=log(as.numeric(lambda)), y= as.numeric(joint), color= design, linetype=design))+
  geom_line(show.legend = TRUE)+
  xlab(latex2exp::TeX("$\\lambda"))+
  theme_bw()+
  scale_color_manual(values =c("red","red", 'blue', "orange"), labels = c(latex2exp::TeX("$HILS_{\\Lambda}$"),latex2exp::TeX("$HILS_{max}$"), "PED"))+
  scale_linetype_manual(values =c("solid", "dashed", "solid"), labels = c(latex2exp::TeX("$HILS_{\\Lambda}$"),latex2exp::TeX("$HILS_{max}$"), "PED"))+
  theme(axis.title.y=element_text( vjust=0.55), legend.title = element_blank(), text=element_text(size=16), legend.position = c(0.2,.75))



# ggall= ggarrange(ggjoint,ggarrange(ggjoint_s4, ggjoint_mixed, ncol=2, labels=c("B","C")), nrow=2,labels="A")
# # # ggall
#  ggsave("./Section4_n9_k10_joint_figure_pre_PP2.png", dpi=600, device="png", ggall) ## this saves a 600 dpi image.
# # # 
# # 













# F_0=d2
# 
# V_0_half= diag(lapply(as.data.frame(F_0),fn))
# F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
# for (i in c(1:length(lambda))){
#   print(i)
#   lam = lambda[i]
#   V_0_half= diag(lapply(as.data.frame(F_0),fn))
#   F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
# 
#   prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                             B_mag=rep(3,k_star),k=k_star,
#                                             sign_vects=all_signs,
#                                             log_lambda = log(lam),
#                                             submodels = NULL,
#                                             log_lambda_strategy="fixed",
#                                             sigma = 1, output_option="I")$results)$mean
#   prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                             B_mag=rep(3,k_star),k=k_star,
#                                             sign_vects=all_signs,
#                                             log_lambda = log(lam),
#                                             submodels = NULL,
#                                             log_lambda_strategy="fixed",
#                                             sigma = 1, output_option="S")$results)$mean
# 
#   results_I<- c(design= "d2", lambda = lam, prob_I)
#   results_S <- c(design="d2", lambda = lam, prob_S)
# 
#   I_results_d2[i,]=results_I
#   S_results_d2[i,]=results_S
# }
# 
# overall_results_d2 = I_results_d2
# overall_results_d2$P_S=as.numeric(S_results_d2$P_S)
# overall_results_d2$joint=as.numeric(I_results_d2$P_I)*as.numeric(S_results_d2$P_S)
# write.csv(overall_results_d2, file ="./n9_k10_kstar3_allsigns_designs_HILS/prob_comparison_d2_n9.csv", row.names =FALSE, col.names = FALSE)









# F_0=HILS_2
# 
# V_0_half= diag(lapply(as.data.frame(F_0),fn))
# F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
# for (i in c(1:length(lambda))){
#   print(i)
#   lam = lambda[i]
#   V_0_half= diag(lapply(as.data.frame(F_0),fn))
#   F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
# 
#   prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                             B_mag=rep(3,k_star),k=k_star,
#                                             sign_vects=all_signs,
#                                             log_lambda = log(lam),
#                                             submodels = NULL,
#                                             log_lambda_strategy="fixed",
#                                             sigma = 1, output_option="I")$results)$mean
#   prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                             B_mag=rep(3,k_star),k=k_star,
#                                             sign_vects=all_signs,
#                                             log_lambda = log(lam),
#                                             submodels = NULL,
#                                             log_lambda_strategy="fixed",
#                                             sigma = 1, output_option="S")$results)$mean
# 
#   results_I<- c(design= "HILS 2", lambda = lam, prob_I)
#   results_S <- c(design="HILS 2", lambda = lam, prob_S)
# 
#   I_results_HILS_2[i,]=results_I
#   S_results_HILS_2[i,]=results_S
# }
# 
# 
# overall_results_HILS_2 = I_results_HILS_2
# overall_results_HILS_2$P_S=as.numeric(S_results_HILS_2$P_S)
# overall_results_HILS_2$joint=as.numeric(I_results_HILS_2$P_I)*as.numeric(S_results_HILS_2$P_S)
# write.csv(overall_results_HILS_2, file ="./n9_k10_kstar3_allsigns_designs_HILS/prob_comparison_HILS_2_n9.csv", row.names =FALSE, col.names = FALSE)

