library(readr)
library(doParallel)
library(foreach)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(purrr)
source('~/SIC_Paper/SIC_Paper_Function_Library.R')

# d1<-as.matrix(read_in_design("./n9_k10_kstar3_allpos_designs_HILS/d1.txt"))
d1_n9<-as.matrix(read_in_design("./n9_k10_designs/d1.txt"))
# #d2<-as.matrix(read_in_design("./n9_k10_kstar3_allpos_designs_HILS/d2.txt"))
# 
HILS_opt<- as.matrix(read_csv("./HILS_opt_n9_k10_kstar3_allsigns.csv"))
#HILS_opt<- as.matrix(read_csv("./n9_k10_kstar3_allsigns_designs_HILS/HILS_opt_n9_k10_kstar3_allsigns.csv"))
HILS_lam_opt <- as.matrix(read_csv("./n9_k10_kstar3_allsigns_designs_HILS/best_huer_Var_s_design_14.csv"))
HILS_lam_opt <- as.matrix(read_csv("./best_huer_Var_s_design_14.csv"))

# true_y<-function(X, k=10, sn=3, error=1,active=c(1,2,3),R=c(1,1,1)){
#   effects=1:k
#   n = dim(X)[1]
#   a = length(active)
#   #if(size=="med") { a=floor(n/3) } else if (size=="hi") {a=3} else {a=floor(0.75*n)}
#   # active=sample(k, a, replace = FALSE, prob = NULL)
#   true=effects%in%active
#   inactive=which(true %in% c(0))
#   #R=(rbinom(a,1,p)*-2)+1 #creates a vector of 1 and -1 of signs
#   #acoeff=(rexp(a, rate=1)+sn)*R
#   acoeff=(rep(sn,a))*R
#   inacoeff=rep(0,(k-a))
#  
#   Y=X[,active]%*%t(as.matrix(acoeff))+X[,inactive]%*%inacoeff+rnorm(n,mean=0, sd=error)
#   
#   return(list(Y=Y, active=active, inactive=inactive))
#   
# }
# 
# 
# ### Fit lasso for one value of lambda
# ### X=design matrix, linear effects
# ### Y=true response vector
# ### lambda=value of lambda for the fit
# 
# lasso.fit<-function(X, Y){
#   require(broom)
#   require(glmnet)
#   n = dim(X)[1]
#   Yc=Y-mean(Y)
#   Xc<-apply(X, 2, function(X) X - mean(X))
#   vlength<-sqrt((1/n)*(diag(t(Xc)%*%Xc)))
#   Xs<-t(t(Xc)*(1/(vlength)))
#   #log_lam = seq(from = -4.5, to= 2, by = 0.1)
#   log_lam = log(seq(exp(-4.5), exp(2), by=exp(log(0.1))))
#   lambda = exp(log_lam)
#   fit<-tidy(glmnet(Xs, Yc, standardize = FALSE, intercept = FALSE, 
#                    lambda=lambda), return_zeros=TRUE)
#   
#   return(fit)
# }
# 
# 
# #Function to permute columns
# permute_cols <- function(X){
#   X_out= X[,sample(ncol(X))]
#   return(X_out)
# }
# take_median<- function(df){
#   median(as.matrix(df))
# }
# take_mean<- function(df){
#   mean(as.matrix(df))
# }
# get_quantile_25<- function(df){
#   quantile(df, na.rm = T, probs=0.25)
# }
# get_quantile_75<- function(df){
#   quantile(df, na.rm = T, probs=0.75)}
# 
# construct_beta_target<-function(nlambdas,k=10, sn=3,active=c(1,2,3),R=c(1,1,1)){
#   effects=1:k
#   a = length(active)
#   #if(size=="med") { a=floor(n/3) } else if (size=="hi") {a=3} else {a=floor(0.75*n)}
#   # active=sample(k, a, replace = FALSE, prob = NULL)
#   true=effects%in%active
#   inactive=which(true %in% c(0))
#   #R=(rbinom(a,1,p)*-2)+1 #creates a vector of 1 and -1 of signs
#   #acoeff=(rexp(a, rate=1)+sn)*R
#   acoeff=(rep(sn,a))*as.numeric(as.vector(R))
#   inacoeff=rep(0,(k-a))
#   target_beta= c(acoeff, inacoeff)
#   #target_beta_rep = as.matrix(rep(target_beta, each = nlambdas), nrow=nlambdas, ncol=1)
#   target_beta_rep = rep(target_beta, each=nlambdas)
#   
#   return(target_beta_rep)
# }
# 
# all_signs = gen_all_sign_vects(c(1,2,3))
# resp_d1 <-true_y(d1, R=all_signs[1,])
# # resp_d2<- true_y(d2)
# resp_HILS<-true_y(HILS_opt, R=all_signs[1,])
# # resp_HILS_2<- true_y(HILS_2)
# #run glmnet for each design
# fit_d1<-lasso.fit(d1, Y=resp_d1$Y)
# # fit_d2<-lasso.fit(d2, Y=resp_d2$Y)
# fit_HILS <- lasso.fit(HILS_opt, Y=resp_HILS$Y)
# # fit_HILS_2 <- lasso.fit(HILS_2, Y=resp_HILS_2$Y)
# 
# fit_d1<- fit_d1[fit_d1$term!="(Intercept)",]
# # fit_d2<- fit_d2[fit_d2$term!="(Intercept)",]
# fit_HILS<- fit_HILS[fit_HILS$term!="(Intercept)",]
# # fit_HILS_2<- fit_HILS_2[fit_HILS_2$term!="(Intercept)",]
# 
# # fit_UE<- fit_UE[order(fit_UE$term),]
# # fit_VarS<- fit_VarS[order(fit_VarS$term),]
# 
# ### Initialize coefficient paths df for each design, these are p by number of lambda dfs
# df_fit_d1_all <-fit_d1%>%dplyr::select(term, lambda)
# # df_fit_d2_all <-fit_d2%>%dplyr::select(term, lambda)
# df_fit_HILS_all <-fit_HILS%>%dplyr::select(term, lambda)
# # df_fit_HILS_2_all <-fit_HILS_2%>%dplyr::select(term, lambda)
# 
# count=1
# for (s in c(1:nrow(all_signs))){
#   target_effects= construct_beta_target(nlambdas = 74, R= all_signs[s,])
#   #browser()
#   
#   for(i in c(1:50)){
#     #print(i)
#     #permute the cols for each design
#     d1_perm <- permute_cols(d1)
#     # d2_perm <- permute_cols(d2)
#     HILS_perm<-permute_cols(HILS_opt)
#     # HILS_2_perm<-permute_cols(HILS_2)
#     
#     # make the i specific coefficient path dfs for each design matrix of all zeros
#     for (j in c(1:30)){
#       #generate y for each design
#       resp_d1 <-true_y(d1_perm, R=all_signs[s,])
#       # resp_d2 <-true_y(d2_perm)
#       resp_HILS<-true_y(HILS_perm, R=all_signs[s,])
#       # resp_HILS_2<-true_y(HILS_2_perm)
#       #run glmnet for each design
#       fit_d1<-lasso.fit(d1_perm, Y=resp_d1$Y)
#       # fit_d2<-lasso.fit(d2_perm, Y=resp_d2$Y)
#       fit_HILS<-lasso.fit(HILS_perm, Y=resp_HILS$Y)
#       # fit_HILS_2<-lasso.fit(HILS_2_perm, Y=resp_HILS_2$Y)
#       
#       fit_d1<- fit_d1[fit_d1$term!="(Intercept)",]
#       # fit_d2<- fit_d2[fit_d2$term!="(Intercept)",]
#       fit_HILS<- fit_HILS[fit_HILS$term!="(Intercept)",]
#       # fit_HILS_2<- fit_HILS_2[fit_HILS_2$term!="(Intercept)",]
#       
#       # fit_UE<- fit_UE[order(fit_UE$term),]
#       # fit_VarS<- fit_VarS[order(fit_VarS$term),]
#       #add coeff path to the i specific coeff path dfs.
#       #browser()
#       df_fit_d1_all= cbind(df_fit_d1_all, abs(fit_d1$estimate-target_effects))
#       names(df_fit_d1_all)[length(names(df_fit_d1_all))]<-paste("estimation_error",count)
#       
#       df_fit_HILS_all= cbind(df_fit_HILS_all, abs(fit_HILS$estimate-target_effects))
#       names(df_fit_HILS_all)[length(names(df_fit_HILS_all))]<-paste("estimation_error",count)
#       count=count+1
#     }
#   }
# }
# df_fit_d1_all = df_fit_d1_all%>%group_by(term,lambda)%>%nest()
# df_fit_HILS_all= df_fit_HILS_all%>%group_by(term, lambda)%>%nest()
# 
# 
# # df_fit_d1_all =df_fit_d1_all%>% 
# #   mutate(median_est=map_dbl(.data[["data"]], take_median))%>%
# #   mutate(mean_est=map_dbl(.data[["data"]], take_mean))%>%
# #   mutate(percentile_25=map_dbl(.data[["data"]], get_quantile_25))%>%
# #   mutate(percentile_75=map_dbl(.data[["data"]], get_quantile_75))
# # 
# # df_fit_HILS_all =df_fit_HILS_all%>% 
# #   mutate(median_est=map_dbl(.data[["data"]], take_median))%>%
# #   mutate(mean_est=map_dbl(.data[["data"]], take_mean))%>%
# #   mutate(percentile_25=map_dbl(.data[["data"]], get_quantile_25))%>%
# #   mutate(percentile_75=map_dbl(.data[["data"]], get_quantile_75))
# 
# avg_coeff_path_d1 = df_fit_d1_all
# 
# 
# avg_coeff_path_HILS = df_fit_HILS_all
# 
# avg_coeff_path_d1$design="d1"
# subset_d1=avg_coeff_path_d1[1:222,]
# subset_I_d1=avg_coeff_path_d1[223:740,]
# subset_I_d1$design="Inactive Effects d1"
# 
# 
# avg_coeff_path_HILS$design="HILS"
# avg_coeff_path_HILS$term = paste0(avg_coeff_path_HILS$term, "HILS")
# subset_HILS=avg_coeff_path_HILS[1:222,]
# subset_I_HILS=avg_coeff_path_HILS[223:740,]
# subset_I_HILS$design="Inactive Effects HILS"
# 
# 
# 
# 
# 
# 
# 
# 
# #combined_subset = rbind(subset_d1, subset_d2, subset_HILS, subset_HILS_2, subset_I_d1, subset_I_d2, subset_I_HILS, subset_I_HILS_2)
# combined_subset = rbind(subset_d1, subset_HILS, subset_I_d1, subset_I_HILS)
# #combined_subset$design= factor(combined_subset$design, levels= c("HILS", "d1", "d2","HILS 2", "Inactive Effects HILS", "Inactive Effects d1", "Inactive Effects d2","Inactive Effects HILS 2" ))
# combined_subset$design= factor(combined_subset$design, levels= c("HILS", "d1",  "Inactive Effects HILS", "Inactive Effects d1"))
# combined_subset_raw= combined_subset%>%unnest(c("data"))
# active_subset= ungroup(combined_subset_raw[1:444,])%>%group_by(lambda, design)%>%
#     summarise(across('estimation_error 1':"estimation_error 6000", ~ mean(.x, na.rm = TRUE)), .groups="keep")%>%
#     nest()%>%
#     mutate(median_est=map_dbl(.data[["data"]], take_median))%>%
#     mutate(percentile_25=map_dbl(.data[["data"]], get_quantile_25))%>%
#     mutate(percentile_75=map_dbl(.data[["data"]], get_quantile_75))
# #   
# inactive_subset= combined_subset_raw[445:1480,]%>%group_by(lambda, design)%>%
#   summarise(across('estimation_error 1':"estimation_error 6000", ~ mean(.x, na.rm = TRUE)), .groups="keep")%>%
#   nest()%>%
#   mutate(median_est=map_dbl(.data[["data"]], take_median))%>%
#   mutate(percentile_25=map_dbl(.data[["data"]], get_quantile_25))%>%
#   mutate(percentile_75=map_dbl(.data[["data"]], get_quantile_75))
# 
# 
# save(combined_subset, file ="./n9_k10_kstar3_allsigns_designs_HILS/coeff_path_comparison_raw_n9_k10.RData")
# 
# 
# 
# # trimmed_subset= combined_subset[combined_subset$lambda<6.0,]
# # 
# # 
# ggactive=ggplot(active_subset[ active_subset$lambda<5,], aes(lambda, median_est, group = design, colour=design)) +
#   scale_color_manual(values =c("red", "blue"))+
#   geom_ribbon(aes(x= lambda, ymin= percentile_25, ymax= percentile_75, group = design, fill=design),inherit.aes = FALSE,  alpha=0.50, linetype=3, color=NA, show.legend = FALSE)+
#   geom_line(show.legend = TRUE)+
#   geom_line(aes(x=lambda, y=percentile_25,group=design, colour=design), linetype=2, show.legend=FALSE)+
#   geom_line(aes(x=lambda, y=percentile_75,group=design, colour=design), linetype=2, show.legend=FALSE)+
#   ylim(-.1, 3)+
#   ggtitle("Simulated Solution Path IQR Active Effects")+
#   #ylab(expression(paste("median(",hat(beta),")")))+
#   ylab(expression(hat(beta)))+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   theme_bw()+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55), legend.title = element_blank())
# 
# # 
# gginactive=ggplot(inactive_subset[inactive_subset$lambda<5,], aes(lambda, median_est, group = design, colour=design)) +
#   scale_color_manual(labels=c("HILS", "d1"),values =c("red" ,"blue"))+
#   geom_ribbon(aes(x= lambda, ymin= percentile_25, ymax= percentile_75, group = design, fill=design),inherit.aes = FALSE, alpha=0.50, linetype=3, color=NA, show.legend = FALSE)+
#   #geom_line(show.legend = FALSE)+
#   geom_line()+
#   geom_line(aes(x=lambda, y=percentile_25,group=design, colour=design), linetype=2)+
#   geom_line(aes(x=lambda, y=percentile_75,group=design, colour=design), linetype=2)+
#   #ylim(-.1/12, .25)+
#   #ylim(-.1, 3)+
#   ggtitle("Simulated Solution Path IQR Inactive Effects")+
#   #ylab(expression(paste("median(",hat(beta),")")))+
#   ylab(expression(hat(beta)))+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   theme_bw()+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55), legend.title = element_blank())+
#   scale_y_continuous(limits=c(-.1/12, .5), breaks= c(0,0.1,0.2,0.3,0.4,0.5))
# 
# 
# ggSP=ggarrange(ggactive, gginactive, ncol = 1, nrow=2)



#ggsave("./Section4_n14_k20_soln_path_figure.png", dpi=600, device="png", ggSP) ## this saves a 600 dpi image.



n=9
k=10
k_star=3
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)
all_signs = gen_all_sign_vects(c(1,2,3))
#log_lam = seq(from = -4.5, to= 2, by = 0.1)
# submodel_samples = stack_NBIBD(p=k, k=k_star, s_1=64, num_stacks = 10)
# A_val = submodel_samples$A_final
# A_64 = submodel_samples$A_NBIBD

#log_lam = log(seq(exp(-4.5), exp(2), by=exp(log(0.1))))
log_lam = seq(-4.5, 2, by = 0.05)
#lambda = exp(log_lam)

# I_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0), sign=character(0))
# S_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0),sign=character(0))
# joint_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0),sign=character(0))
# 
# I_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0),sign=character(0))
# S_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0),sign=character(0))
# joint_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0),sign=character(0))
I_results_HILS_2 = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0),sign=character(0))
S_results_HILS_2 = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0),sign=character(0))
joint_results_HILS_2 = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0),sign=character(0))

# I_results_d1$design = factor(I_results_d1$design, levels=c("d1","HILS"))
# S_results_d1$design = factor(S_results_d1$design, levels=c("d1","HILS"))
# joint_results_d1$design = factor(S_results_d1$design, levels=c("d1","HILS"))
# I_results_d1$sign = factor(I_results_d1$sign, levels=c("s_1","s_2","s_3","s_4"))
# S_results_d1$sign = factor(S_results_d1$sign, levels=c("s_1","s_2","s_3","s_4"))
# joint_results_d1$sign = factor(S_results_d1$sign, levels=c("s_1","s_2","s_3","s_4"))
# 
# I_results_HILS$design = factor(I_results_HILS$design, levels=c("d1","HILS"))
# S_results_HILS$design = factor(S_results_HILS$design, levels=c("d1","HILS"))
# joint_results_HILS$design = factor(S_results_HILS$design, levels=c("d1","HILS"))
# I_results_HILS$sign = factor(I_results_HILS$sign, levels=c("s_1","s_2","s_3","s_4"))
# S_results_HILS$sign = factor(S_results_HILS$sign, levels=c("s_1","s_2","s_3","s_4"))
# joint_results_HILS$sign = factor(S_results_HILS$sign, levels=c("s_1","s_2","s_3","s_4"))


I_results_HILS_2$design = factor(I_results_HILS_2$design, levels=c("d1","HILS","HILS2"))
S_results_HILS_2$design = factor(S_results_HILS_2$design, levels=c("d1","HILS", "HILS2"))
joint_results_HILS_2$design = factor(S_results_HILS_2$design, levels=c("d1","HILS", "HILS2"))
I_results_HILS_2$sign = factor(I_results_HILS_2$sign, levels=c("s_1","s_2","s_3","s_4"))
S_results_HILS_2$sign = factor(S_results_HILS_2$sign, levels=c("s_1","s_2","s_3","s_4"))
joint_results_HILS_2$sign = factor(S_results_HILS_2$sign, levels=c("s_1","s_2","s_3","s_4"))


# for (s in c(1:4)){
#   sign_vec =all_signs[s,]
#   print(sign_vec)
#   print(paste0("s_",s))
#   F_0=d1
#   V_0_half= diag(lapply(as.data.frame(F_0),fn))
#   F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
#   for (i in c(1:length(log_lam))){
#     print(i)
#     log_lambda = log_lam[i]
#     V_0_half= diag(lapply(as.data.frame(F_0),fn))
#     F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
#     
#     prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                               B_mag=rep(3,k_star),k=k_star,
#                                               sign_vects=all_signs[s,],
#                                               log_lambda = log_lambda,
#                                               submodels = NULL,
#                                               log_lambda_strategy="fixed",
#                                               sigma = 1, output_option="I")$results)$mean
#     prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                               B_mag=rep(3,k_star),k=k_star,
#                                               sign_vects=all_signs[s,],
#                                               log_lambda = log_lambda,
#                                               submodels = NULL,
#                                               log_lambda_strategy="fixed",
#                                               sigma = 1, output_option="S")$results)$mean
#     prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                                B_mag=rep(3,k_star),k=k_star,
#                                                sign_vects=all_signs[s,],
#                                                log_lambda = log_lambda,
#                                                submodels = NULL,
#                                                log_lambda_strategy="fixed",
#                                                sigma = 1)$results)$mean
#     
#     results_I<- c(design= "d1", lambda = exp(log_lambda), prob_I, paste0("s_",s))
#     results_S <- c(design="d1", lambda = exp(log_lambda), prob_S, paste0("s_",s))
#     results_joint <- c(design="d1", lambda = exp(log_lambda), prob_joint, paste0("s_",s))
#     index = ((s-1)*length(log_lam))+i
#     print(paste0("Index=",index))
#    
#     I_results_d1[index,]=results_I
#     S_results_d1[index,]=results_S
#     joint_results_d1[index,]= results_joint
#   }
#   F_0=HILS_opt
#   
#   V_0_half= diag(lapply(as.data.frame(F_0),fn))
#   F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
#   for (i in c(1:length(log_lam))){
#     print(i)
#     #lam = lambda[i]
#     log_lambda = log_lam[i]
#     V_0_half= diag(lapply(as.data.frame(F_0),fn))
#     F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
#     
#     # prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#     #                                           B_mag=rep(3,k_star),k=k_star,
#     #                                           sign_vects=all_signs[s,],
#     #                                           log_lambda = log(lam),
#     #                                           submodels =NULL,
#     #                                           log_lambda_strategy="fixed",
#     #                                           sigma = 1, output_option="I")$results)$mean
#     # prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#     #                                           B_mag=rep(3,k_star),k=k_star,
#     #                                           sign_vects=all_signs[s,],
#     #                                           log_lambda = log(lam),
#     #                                           submodels = NULL,
#     #                                           log_lambda_strategy="fixed",
#     #                                           sigma = 1, output_option="S")$results)$mean
#     
#     prob_I <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                               B_mag=rep(3,k_star),k=k_star,
#                                               sign_vects=all_signs[s,],
#                                               log_lambda = log_lambda,
#                                               submodels = NULL,
#                                               log_lambda_strategy="fixed",
#                                               sigma = 1, output_option="I")$results)$mean
#     prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                               B_mag=rep(3,k_star),k=k_star,
#                                               sign_vects=all_signs[s,],
#                                               log_lambda = log_lambda,
#                                               submodels = NULL,
#                                               log_lambda_strategy="fixed",
#                                               sigma = 1, output_option="S")$results)$mean
#     prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                                   B_mag=rep(3,k_star),k=k_star,
#                                                   sign_vects=all_signs[s,],
#                                                   log_lambda = log_lambda,
#                                                   submodels = NULL,
#                                                   log_lambda_strategy="fixed",
#                                                   sigma = 1)$results)$mean
#     
#     
#     results_I<- c(design= "HILS", lambda = exp(log_lambda), prob_I, paste0("s_",s))
#     results_S <- c(design="HILS", lambda = exp(log_lambda), prob_S, paste0("s_",s))
#     results_joint <- c(design="HILS", lambda = exp(log_lambda), prob_joint, paste0("s_",s))
#     index = ((s-1)*length(log_lam))+i
#     print(paste0("Index=",index))
#     
#     I_results_HILS[index,]=results_I
#     S_results_HILS[index,]=results_S
#     joint_results_HILS[index,]= results_joint
#   }
#   
# }


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
    
    results_I<- c(design= "d1", lambda = exp(log_lambda), prob_I, paste0("s_",s))
    results_S <- c(design="d1", lambda = exp(log_lambda), prob_S, paste0("s_",s))
    results_joint <- c(design="d1", lambda = exp(log_lambda), prob_joint, paste0("s_",s))
    index = ((s-1)*length(log_lam))+i
    print(paste0("Index=",index))
    
    I_results_HILS_2[index,]=results_I
    S_results_HILS_2[index,]=results_S
    joint_results_HILS_2[index,]= results_joint
  }
 
  
}



# overall_results_d1 = I_results_d1
# overall_results_d1$P_S=as.numeric(S_results_d1$P_S)
# overall_results_d1$joint=as.numeric(I_results_d1$P_I)*as.numeric(S_results_d1$P_S)
# write.csv(overall_results_d1, file ="./n9_k10_kstar3_allsigns_designs_HILS/prob_comparison_d1_n9.csv", row.names =FALSE, col.names = FALSE)
# 
# 
# 
# overall_results_HILS = I_results_HILS
# overall_results_HILS$P_S=as.numeric(S_results_HILS$P_S)
# overall_results_HILS$joint=as.numeric(I_results_HILS$P_I)*as.numeric(S_results_HILS$P_S)
# write.csv(overall_results_HILS, file ="./n9_k10_kstar3_allsigns_designs_HILS/prob_comparison_HILS_n9.csv", row.names =FALSE, col.names = FALSE)









# # I_results=rbind(I_results_d1, I_results_HILS)
# # 
# # S_results=rbind(S_results_d1, S_results_HILS)
# # 
# # joint_results = rbind(joint_results_d1, joint_results_HILS)
# 
# overall_results = I_results
# overall_results$P_S=as.numeric(S_results$P_S)
# overall_results$joint=as.numeric(joint_results$P_joint)
# 
# 
# write.csv(overall_results, file ="./n9_k10_kstar3_allsigns_designs_HILS/prob_comparison_n9_by_sign.csv", row.names =FALSE, col.names = FALSE)


I_results = I_results_HILS_2
S_results = S_results_HILS_2
joint_results = joint_results_HILS_2
overall_results = I_results
overall_results$P_S=as.numeric(S_results$P_S)
overall_results$joint=as.numeric(joint_results$P_joint)
write.csv(overall_results, file ="./n9_k10_kstar3_allsigns_designs_HILS/prob_comparison_lam_opt_n9_by_sign.csv", row.names =FALSE, col.names = FALSE)


print("file written successfully")


library(readr)


overall_results_by_sign1<- read_csv("prob_comparison_n9_by_sign.csv")
overall_results_by_sign2<- read_csv("prob_comparison_lam_opt_n9_by_sign.csv")
overall_results_by_sign2$design="HILS_max_lam"

overall_results_by_sign= rbind(overall_results_by_sign1, overall_results_by_sign2)

overall_results_by_sign[overall_results_by_sign$design=="d1",]$design="PED"




# prob_comparison_d1 <- read_csv("prob_comparison_d1_n9.csv")
# prob_comparison_d2 <- read_csv("prob_comparison_d2_n9.csv")
# prob_comparison_HILS <- read_csv("prob_comparison_HILS_n9.csv")
# prob_comparison_HILS_2 <- read_csv("prob_comparison_HILS_2_n9.csv")
# 
# overall_results= rbind(prob_comparison_d1, prob_comparison_d2, prob_comparison_HILS, prob_comparison_HILS_2)
# 



# 
# trimmed_overall= overall_results[overall_results$lambda<4.5,]
# 
# ggI = ggplot(data = overall_results[overall_results$design!="d2" & overall_results$design!="HILS 2",], aes(x=as.numeric(lambda), y= as.numeric(P_I), color= design))+
#   geom_line(show.legend = FALSE)+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   ylab(expression(paste("P(",I,")")))+
#   scale_color_manual(values =c("blue","red", 'green', "orange"))+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55))
# 
# ggS = ggplot(data = overall_results[overall_results$design!="d2" & overall_results$design!="HILS 2",], aes(x=as.numeric(lambda), y= as.numeric(P_S), color= design))+
#   geom_line(show.legend = FALSE)+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   ylab(expression(paste("P(",italic(S),")")))+
#   #ylab(expression("P(\uD835\uDCAE)"))+
#   scale_color_manual(values =c("blue","red", 'green', "orange"))+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55))
# 
overall_results_mixed_sign= overall_results_by_sign[overall_results_by_sign$sign!="s_4",]
overall_results_mixed_sign_mean= overall_results_mixed_sign%>%group_by(design,lambda)%>%summarise(across(joint, ~ mean(.x, na.rm = TRUE)), .groups="keep")

overall_results=overall_results_by_sign%>%group_by(design,lambda)%>%summarise(across(joint, ~ mean(.x, na.rm = TRUE)), .groups="keep")

overall_results[overall_results$design=="d1",]$design="PED"

overall_results = overall_results[overall_results$design!="d2",]

ggjoint = ggplot(data = overall_results[overall_results$design!="HILS_max_lam",], aes(x=log(as.numeric(lambda)), y= as.numeric(joint), color= design, linetype=design))+
  geom_line()+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  scale_color_manual(values =c("red","blue"), labels = c(latex2exp::TeX("$HILS_{\\Lambda} / HILS_{max}$"), "PED"))+
  scale_x_continuous(limits = c(-2,2))+
  scale_linetype_manual(values =c("dashed", "solid"), labels = c(latex2exp::TeX("$HILS_{\\Lambda} / HILS_{max}$"), "PED"))+
  theme_bw()+
  theme(axis.title.y=element_text( vjust=0.55), legend.title = element_blank(), text=element_text(size=16), legend.position = c(0.2,.75))

# ggjoint_s1 = ggplot(data = overall_results_by_sign[overall_results_by_sign$sign=="s_1",], aes(x=as.numeric(lambda), y= as.numeric(joint), color= design))+
#   geom_line()+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   ylab(latex2exp::TeX("$\\phi_l(\\textbf{z}_A=(-1,-1,1)^T, \\lambda)"))+
#   scale_color_manual(values =c("blue","red", 'green', "orange"))+
#   theme_bw()+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55), legend.title = element_blank())
# 
#   
# ggjoint_s2 = ggplot(data = overall_results_by_sign[overall_results_by_sign$sign=="s_2",], aes(x=as.numeric(lambda), y= as.numeric(joint), color= design))+
#     geom_line()+
#     xlab(latex2exp::TeX("$\\lambda"))+
#     ylab(latex2exp::TeX("$\\phi_l(\\textbf{z}_A=(1,-1,1)^T, \\lambda)"))+
#     scale_color_manual(values =c("blue","red", 'green', "orange"))+
#     theme_bw()+
#     theme(axis.title.y=element_text(angle=0, vjust=0.55), legend.title = element_blank())
# 
# ggjoint_s3 = ggplot(data = overall_results_by_sign[overall_results_by_sign$sign=="s_3",], aes(x=as.numeric(lambda), y= as.numeric(joint), color= design))+
#   geom_line()+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   ylab(latex2exp::TeX("$\\phi_l(\\textbf{z}_A=(-1,1,1)^T, \\lambda)"))+
#   scale_color_manual(values =c("blue","red", 'green', "orange"))+
#   theme_bw()+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55), legend.title = element_blank())


ggjoint_s4 = ggplot(data = overall_results_by_sign[overall_results_by_sign$sign=="s_4" & overall_results_by_sign$lambda >= .12,], aes(x=log(as.numeric(lambda)), y= as.numeric(joint), color= design, linetype=design))+
  geom_line(show.legend = TRUE)+
  xlab(latex2exp::TeX("$\\lambda"))+
  theme_bw()+
  #ylab(latex2exp::TeX("$\\phi_l(\\textbf{z}_A=(1,1,1)^T, \\lambda)"))+
  #ylab(" ")+
  scale_color_manual(values =c("red","red", 'blue', "orange"), labels = c(latex2exp::TeX("$HILS_{\\Lambda}$"),latex2exp::TeX("$HILS_{max}$"), "PED"))+
  scale_linetype_manual(values =c("solid", "dashed", "solid"), labels = c(latex2exp::TeX("$HILS_{\\Lambda}$"),latex2exp::TeX("$HILS_{max}$"), "PED"))+
  # scale_fill_discrete(labels = c(latex2exp::TeX("$HILS_{\\Lambda}$"),latex2exp::TeX("$HILS_{max}$"), "PED"))+
  theme(axis.title.y=element_text( vjust=0.55), legend.title = element_blank(), text=element_text(size=16), legend.position = c(0.2,.75))

ggjoint_mixed = ggplot(data = overall_results_mixed_sign_mean[overall_results_mixed_sign_mean$lambda<4.5,], aes(x=log(as.numeric(lambda)), y= as.numeric(joint), color= design))+
  geom_line()+
  xlab(latex2exp::TeX("$\\lambda"))+
  #ylab(latex2exp::TeX("$\\phi_l(\\textbf{z}_A=(-1,1,1)^T, \\lambda)"))+
  #ylab(" ")+
  scale_color_manual(values =c("red","blue", 'green', "orange"))+
  theme_bw()+
  theme(axis.title.y=element_text( vjust=0.55), legend.title = element_blank(), text=element_text(size=14))

ggall= ggarrange(ggjoint,ggarrange(ggjoint_s4, ggjoint_mixed, ncol=2, labels=c("B","C")), nrow=2,labels="A")
# # ggall
 ggsave("./Section4_n9_k10_joint_figure_pre_PP2.png", dpi=600, device="png", ggall) ## this saves a 600 dpi image.
# # 
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

