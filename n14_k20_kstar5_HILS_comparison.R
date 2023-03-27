library(readr)
library(doParallel)
library(foreach)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(purrr)
source('~/SIC_Paper/SIC_Paper_Function_Library.R')

d1<-as.matrix(read_in_design("./n14_k20_kstar5_allpos_designs_HILS/d1.txt"))
d2<-as.matrix(read_in_design("./n14_k20_kstar5_allpos_designs_HILS/d2.txt"))

HILS_2<- as.matrix(read_csv("./n14_k20_kstar5_allpos_designs_HILS_lam_opt/best_huer_Var_s_design_15.csv"))
HILS_opt <- as.matrix(read_csv("./n14_k20_kstar5_allpos_designs_HILS_lam_opt/best_huer_Var_s_design_18.csv"))
HILS_opt <- as.matrix(read_csv("./best_huer_Var_s_design_18.csv"))
#HILS_opt <- as.matrix(read_csv("./best_huer_Var_s_design_18.csv"))
true_y<-function(X, k=14, p=0, sn=3, error=1,active=c(1,2,3,4,5)){
  effects=1:k
  n = dim(X)[1]
  a = length(active)
  #if(size=="med") { a=floor(n/3) } else if (size=="hi") {a=3} else {a=floor(0.75*n)}
  # active=sample(k, a, replace = FALSE, prob = NULL)
  true=effects%in%active
  inactive=which(true %in% c(0))
  R=(rbinom(a,1,p)*-2)+1 #creates a vector of 1 and -1 of signs
  #acoeff=(rexp(a, rate=1)+sn)*R
  acoeff=(rep(sn,a))*R
  inacoeff=rep(0,(k-a))
  Y=X[,active]%*%acoeff+X[,inactive]%*%inacoeff+rnorm(n,mean=0, sd=error)
  
  return(list(Y=Y, active=active, inactive=inactive))
  
}


### Fit lasso for one value of lambda
### X=design matrix, linear effects
### Y=true response vector
### lambda=value of lambda for the fit

lasso.fit<-function(X, Y){
  require(broom)
  require(glmnet)
  n = dim(X)[1]
  Yc=Y-mean(Y)
  Xc<-apply(X, 2, function(X) X - mean(X))
  vlength<-sqrt((1/n)*(diag(t(Xc)%*%Xc)))
  Xs<-t(t(Xc)*(1/(vlength)))
  log_lam = seq(from = -4.5, to= 2, length.out=74)
  #log_lam = log(seq(exp(-4.5), exp(2), by=exp(log(0.1))))
  lambda = exp(log_lam)
  fit<-tidy(glmnet(Xs, Yc, standardize = FALSE, intercept = FALSE, 
                   lambda=lambda), return_zeros=TRUE)
  
  return(fit)
}


#Function to permute columns
permute_cols <- function(X){
  X_out= X[,sample(ncol(X))]
  return(X_out)
}



resp_d1 <-true_y(d1)
# resp_d2<- true_y(d2)
resp_HILS<-true_y(HILS_opt)
# resp_HILS_2<- true_y(HILS_2)
#run glmnet for each design
fit_d1<-lasso.fit(d1, Y=resp_d1$Y)
# fit_d2<-lasso.fit(d2, Y=resp_d2$Y)
fit_HILS <- lasso.fit(HILS_opt, Y=resp_HILS$Y)
# fit_HILS_2 <- lasso.fit(HILS_2, Y=resp_HILS_2$Y)

fit_d1<- fit_d1[fit_d1$term!="(Intercept)",]
# fit_d2<- fit_d2[fit_d2$term!="(Intercept)",]
fit_HILS<- fit_HILS[fit_HILS$term!="(Intercept)",]
# fit_HILS_2<- fit_HILS_2[fit_HILS_2$term!="(Intercept)",]

# fit_UE<- fit_UE[order(fit_UE$term),]
# fit_VarS<- fit_VarS[order(fit_VarS$term),]

### Initialize coefficient paths df for each design, these are p by number of lambda dfs
# df_fit_d1_all <-fit_d1%>%dplyr::select(term, lambda)
# # df_fit_d2_all <-fit_d2%>%dplyr::select(term, lambda)
# df_fit_HILS_all <-fit_HILS%>%dplyr::select(term, lambda)
# # df_fit_HILS_2_all <-fit_HILS_2%>%dplyr::select(term, lambda)
# count=1
# for(i in c(1:50)){
#   #print(i)
#   #permute the cols for each design
#   d1_perm <- permute_cols(d1)
#   # d2_perm <- permute_cols(d2)
#   HILS_perm<-permute_cols(HILS_opt)
#   # HILS_2_perm<-permute_cols(HILS_2)
# 
#   # make the i specific coefficient path dfs for each design matrix of all zeros
#   for (j in c(1:30)){
#     #generate y for each design
#     resp_d1 <-true_y(d1_perm)
#     # resp_d2 <-true_y(d2_perm)
#     resp_HILS<-true_y(HILS_perm)
#     # resp_HILS_2<-true_y(HILS_2_perm)
#     #run glmnet for each design
#     fit_d1<-lasso.fit(d1_perm, Y=resp_d1$Y)
#     # fit_d2<-lasso.fit(d2_perm, Y=resp_d2$Y)
#     fit_HILS<-lasso.fit(HILS_perm, Y=resp_HILS$Y)
#     # fit_HILS_2<-lasso.fit(HILS_2_perm, Y=resp_HILS_2$Y)
# 
#     fit_d1<- fit_d1[fit_d1$term!="(Intercept)",]
#     # fit_d2<- fit_d2[fit_d2$term!="(Intercept)",]
#     fit_HILS<- fit_HILS[fit_HILS$term!="(Intercept)",]
#     # fit_HILS_2<- fit_HILS_2[fit_HILS_2$term!="(Intercept)",]
# 
#     # fit_UE<- fit_UE[order(fit_UE$term),]
#     # fit_VarS<- fit_VarS[order(fit_VarS$term),]
#     #add coeff path to the i specific coeff path dfs.
#     df_fit_d1_all= cbind(df_fit_d1_all, fit_d1$estimate)
#     names(df_fit_d1_all)[length(names(df_fit_d1_all))]<-paste("estimate",count)
#     
#     df_fit_HILS_all= cbind(df_fit_HILS_all, fit_HILS$estimate)
#     names(df_fit_HILS_all)[length(names(df_fit_HILS_all))]<-paste("estimate",count)
#     count=count+1
#   }
# }
# 
# df_fit_d1_all = df_fit_d1_all%>%group_by(term,lambda)%>%nest()
# df_fit_HILS_all= df_fit_HILS_all%>%group_by(term, lambda)%>%nest()
# 
# take_median<- function(df){
#   median(as.matrix(df))
# }
# take_mean<- function(df){
#   mean(as.matrix(df))
# }
# 
# 
# 
# df_fit_d1_all =df_fit_d1_all%>% 
#   mutate(median_est=map_dbl(.data[["data"]], take_median))%>%
#   mutate(mean_est=map_dbl(.data[["data"]], take_mean))
# 
# df_fit_HILS_all =df_fit_HILS_all%>% 
#   mutate(median_est=map_dbl(.data[["data"]], take_median))%>%
#   mutate(mean_est=map_dbl(.data[["data"]], take_mean))
# 
# avg_coeff_path_d1 = df_fit_d1_all
# 
# 
# avg_coeff_path_HILS = df_fit_HILS_all


# avg_coeff_path_d1$estimate = avg_coeff_path_d1$estimate/(50*100)
# # avg_coeff_path_d2$estimate = avg_coeff_path_d2$estimate/(50*100)
# 
# avg_coeff_path_HILS$estimate = avg_coeff_path_HILS$estimate/(50*100)
# avg_coeff_path_HILS_2$estimate = avg_coeff_path_HILS_2$estimate/(50*100)



# ggplot(avg_coeff_path_UE, aes(lambda, estimate, group = term, colour=term)) +
#   geom_line()+
#   ylim(-.1, 3)+
#   scale_color_manual(values=(c("red","black","black","black","black","black","black","black", "blue", "green", "orange","black","black","black","black","black")))+
#   ggtitle("UE(s^2) design")
#
#
# ggplot(avg_coeff_path_VarS, aes(lambda, estimate, group = term, colour=term)) +
#   geom_line()+
#   ylim(-.1, 3)+
#   scale_color_manual(values=(c("red","black","black","black","black","black","black","black", "blue", "green", "orange","black","black","black","black","black")))+
#   ggtitle("Best design")+
#   ylab("Estimate")+
#   xlab("Lambda")
#
# trimmed_df_UE=avg_coeff_path_UE[avg_coeff_path_UE$lambda>2.99,]
#
# ggplot(trimmed_df_UE, aes(lambda, estimate, group = term, colour=term)) +
#   geom_line()+
#   ylim(-.1, 0.8)+
#   scale_color_manual(values=(c("red","black","black","black","black","black","black","black", "blue", "green", "orange","black","black","black","black","black")))+
#   ggtitle("UE(s^2) design")+
#   ylab("Estimate")+
#   xlab("Lambda")
#
# trimmed_df=avg_coeff_path_VarS[avg_coeff_path_VarS$lambda>2.99,]
#
# ggplot(trimmed_df, aes(lambda, estimate, group = term, colour=term)) +
#   geom_line()+
#   ylim(-.1, 0.8)+
#   scale_color_manual(values=(c("red","black","black","black","black","black","black","black", "blue", "green", "orange","black","black","black","black","black")))+
#   ggtitle("Best design")+ylab("Estimate")+
#   xlab("Lambda")



# avg_coeff_path_d1$design="d1"
# subset_d1=avg_coeff_path_d1[1:370,]
# subset_I_d1=avg_coeff_path_d1[371:1480,]
# subset_I_d1$design="Inactive Effects d1"
# 
# # avg_coeff_path_d2$design="d2"
# # avg_coeff_path_d2$term = paste0(avg_coeff_path_d2$term, "d2")
# # subset_d2=avg_coeff_path_d2[1:370,]
# # subset_I_d2=avg_coeff_path_d2[371:1480,]
# # subset_I_d2$design="Inactive Effects d2"
# 
# avg_coeff_path_HILS$design="HILS"
# avg_coeff_path_HILS$term = paste0(avg_coeff_path_HILS$term, "HILS")
# subset_HILS=avg_coeff_path_HILS[1:370,]
# subset_I_HILS=avg_coeff_path_HILS[371:1480,]
# subset_I_HILS$design="Inactive Effects HILS"

# avg_coeff_path_HILS_2$design="HILS 2"
# avg_coeff_path_HILS_2$term = paste0(avg_coeff_path_HILS_2$term, "HILS_2")
# subset_HILS_2=avg_coeff_path_HILS_2[1:370,]
# subset_I_HILS_2=avg_coeff_path_HILS_2[371:1480,]
# subset_I_HILS_2$design="Inactive Effects HILS 2"




#combined_subset = rbind(subset_d1, subset_d2, subset_HILS, subset_HILS_2, subset_I_d1, subset_I_d2, subset_I_HILS, subset_I_HILS_2)
#combined_subset = rbind(subset_d1, subset_HILS, subset_I_d1, subset_I_HILS)
#combined_subset$design= factor(combined_subset$design, levels= c("HILS", "d1", "d2","HILS 2", "Inactive Effects HILS", "Inactive Effects d1", "Inactive Effects d2","Inactive Effects HILS 2" ))
#combined_subset$design= factor(combined_subset$design, levels= c("HILS", "d1",  "Inactive Effects HILS", "Inactive Effects d1"))
#head(combined_subset)
#save(combined_subset, file ="./n14_k20_kstar5_allpos_designs_HILS/coeff_path_comparison_raw_log.RData")
# combined_subset = dplyr::select(combined_subset, -data)
# write.csv(combined_subset, file ="./n14_k20_kstar5_allpos_designs_HILS/coeff_path_comparison.csv", row.names =FALSE, col.names = FALSE)
# 

#Start plotting here


# combined_subset_raw= combined_subset%>%unnest(c("data"))
# library(purrr)
# get_quantile_25<- function(df){
#   quantile(df, na.rm = T, probs=0.25)
# }
# get_quantile_75<- function(df){
#   quantile(df, na.rm = T, probs=0.75)
# }
# 
# active_subset= ungroup(combined_subset_raw[1:740,])%>%group_by(lambda, design)%>%
#   summarise(across('estimate 1':"estimate 2500", ~ mean(abs(.x-3), na.rm = TRUE)), .groups="keep")%>%
#   nest()%>%
#   mutate(median_est=map_dbl(.data[["data"]], take_median))%>%
#   mutate(percentile_25=map_dbl(.data[["data"]], get_quantile_25))%>%
#   mutate(percentile_75=map_dbl(.data[["data"]], get_quantile_75))
# 
# inactive_subset= combined_subset_raw[741:2960,]%>%group_by(lambda, design)%>%
#   summarise(across('estimate 1':"estimate 2500", ~ mean(.x, na.rm = TRUE)), .groups="keep")%>%
#   nest()%>%
#   mutate(median_est=map_dbl(.data[["data"]], take_median))%>%
#   mutate(percentile_25=map_dbl(.data[["data"]], get_quantile_25))%>%
#   mutate(percentile_75=map_dbl(.data[["data"]], get_quantile_75))
# #
# trimmed_subset= combined_subset[combined_subset$lambda<6.0,]
# 
# active_subset$design=as.character(active_subset$design)
# active_subset[active_subset$design=="d1",]$design="PED"
# inactive_subset$design=as.character(inactive_subset$design)
# inactive_subset[inactive_subset$design=="Inactive Effects d1",]$design="PED"
# inactive_subset[inactive_subset$design=="Inactive Effects HILS",]$design="HILS"




# ggactive=ggplot(active_subset[ active_subset$lambda<5,], aes(lambda, median_est, group = design, colour=design)) +
#   scale_color_manual(labels=c("HILS", "PED"),values =c("red" ,"blue"))+
#   geom_ribbon(aes(x= lambda, ymin= percentile_25, ymax= percentile_75, group = design, fill=design),inherit.aes = FALSE, alpha=0.40, linetype=3, color=NA, show.legend = FALSE)+
#   geom_line(size=1.25, show.legend = TRUE)+
#   geom_line(aes(x=lambda, y=percentile_25,group=design, colour=design), linetype=2, show.legend=FALSE)+
#   geom_line(aes(x=lambda, y=percentile_75,group=design, colour=design), linetype=2, show.legend=FALSE)+
#   ylim(-.1, 3)+
#   ggtitle("Simulated Estimation Error IQR Active Effects")+
#   #ylab(expression(paste("median(",hat(beta),")")))+
#   #ylab(expression(hat(beta)))+
#   ylab(expression(paste("|",hat(beta)-beta, "|")))+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   theme_bw()+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55), text=element_text(size=16),legend.title = element_blank())
# 
# 
# gginactive=ggplot(inactive_subset[inactive_subset$lambda<5,], aes(lambda, median_est, group = design, colour=design)) +
#   scale_color_manual(labels=c("HILS", "PED"),values =c("red" ,"blue"))+
#   geom_ribbon(aes(x= lambda, ymin= percentile_25, ymax= percentile_75, group = design, fill=design),inherit.aes = FALSE, alpha=0.4, linetype=3, color=NA, show.legend = FALSE)+
#   #geom_line(show.legend = FALSE)+
#   geom_line(size=1.25)+
#   geom_line(aes(x=lambda, y=percentile_25,group=design, colour=design), linetype=2)+
#   geom_line(aes(x=lambda, y=percentile_75,group=design, colour=design), linetype=2)+
#   ylim(-.1/20, 1)+
#   ggtitle("Simulated Estimation Error IQR Inactive Effects")+
#   #ylab(expression(paste("median(",hat(beta),")")))+
#   ylab(expression(paste("|",hat(beta)-beta, "|")))+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   theme_bw()+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55), legend.title = element_blank(), text=element_text(size=16))
#   # geom_line()+
#   # ylim(-.1/6, 0.5)+
#   # ggtitle("Solution Path Inactive Effects")+ylab(expression(paste("median(",hat(beta),")")))+
#   # xlab(latex2exp::TeX("$\\lambda"))+
#   # theme(axis.title.y=element_text(angle=0, vjust=0.55))+
#   # scale_fill_discrete(name="Design", labels=c("d1", "HILS", "HILS 2"))+
#   # scale_color_manual( labels=c("d1", "HILS", "HILS 2"),values =c("red", "blue", "green","orange","black", "gray", "slategray","lightblue"))
# 
# ggSP=ggarrange(ggactive, gginactive, nrow=2)

# gg1 = ggplot(trimmed_subset, aes(lambda, estimate, group = term, colour=design)) +
#   geom_line()+
#   ylim(-.1, 3)+
#   ggtitle("Solution Path")+ylab(expression(hat(beta)))+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55))+
#   scale_color_manual(values =c("red", "blue", "green","orange","black", "gray", "slategray","lightblue"))

# ggsave("./Section4_n14_k20_soln_path_figure.png", dpi=600, device="png", ggSP) ## this saves a 600 dpi image.



n=14
k=20
k_star=5
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)
log_lam = seq(from = -4.5, to= 2, by = 0.05)
submodel_samples = stack_NBIBD(p=k, k=k_star, s_1=64, num_stacks = 10)
A_val = submodel_samples$A_final
A_64 = submodel_samples$A_NBIBD

#log_lam = log(seq(exp(-4.5), exp(2), by=exp(log(0.1))))
lambda = exp(log_lam)

I_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0))
S_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0))
joint_results_d1 = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0))
# I_results_d2 = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0))
# S_results_d2 = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0))
I_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0))
S_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0))
joint_results_HILS = data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0))
# I_results_HILS_2 = data.frame(design=character(0), lambda=numeric(0), P_I=numeric(0))
# S_results_HILS_2 = data.frame(design=character(0), lambda=numeric(0), P_S=numeric(0))

I_results_d1$design = factor(I_results_d1$design, levels=c("d1","HILS"))
S_results_d1$design = factor(S_results_d1$design, levels=c("d1","HILS"))
joint_results_d1$design = factor(S_results_d1$design, levels=c("d1","HILS"))
# I_results_d2$design = factor(I_results_d2$design, levels=c("d1","d2","HILS","HILS 2"))
# S_results_d2$design = factor(S_results_d2$design, levels=c("d1","d2","HILS","HILS 2"))
I_results_HILS$design = factor(I_results_HILS$design, levels=c("d1","HILS"))
S_results_HILS$design = factor(S_results_HILS$design, levels=c("d1","HILS"))
joint_results_HILS$design = factor(S_results_HILS$design, levels=c("d1","HILS"))
# I_results_HILS_2$design = factor(I_results_HILS_2$design, levels=c("d1","d2","HILS","HILS 2"))
# S_results_HILS_2$design = factor(S_results_HILS_2$design, levels=c("d1","d2","HILS","HILS 2"))



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
#overall_results_d1$joint=as.numeric(I_results_d1$P_I)*as.numeric(S_results_d1$P_S)
overall_results_d1$joint = as.numeric(joint_results_d1$P_joint)
write.csv(overall_results_d1, file ="./n14_k20_kstar5_allpos_designs_HILS_lam_opt/prob_comparison_d1_log_lam_opt.csv", row.names =FALSE, col.names = FALSE)

F_0=d2

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
#                                             sign_vects=sign_all_pos,
#                                             log_lambda = log(lam),
#                                             submodels = A_val,
#                                             log_lambda_strategy="fixed",
#                                             sigma = 1, output_option="I")$results)$mean
#   prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                             B_mag=rep(3,k_star),k=k_star,
#                                             sign_vects=sign_all_pos,
#                                             log_lambda = log(lam),
#                                             submodels = A_val,
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
# write.csv(overall_results_d2, file ="./n14_k20_kstar5_allpos_designs_HILS/prob_comparison_d2.csv", row.names =FALSE, col.names = FALSE)



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
#overall_results_HILS$joint=as.numeric(I_results_HILS$P_I)*as.numeric(S_results_HILS$P_S)
write.csv(overall_results_HILS, file ="./n14_k20_kstar5_allpos_designs_HILS_lam_opt/prob_comparison_HILS_log_lam_opt.csv", row.names =FALSE, col.names = FALSE)



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
#                                             sign_vects=sign_all_pos,
#                                             log_lambda = log(lam),
#                                             submodels = A_val,
#                                             log_lambda_strategy="fixed",
#                                             sigma = 1, output_option="I")$results)$mean
#   prob_S <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
#                                             B_mag=rep(3,k_star),k=k_star,
#                                             sign_vects=sign_all_pos,
#                                             log_lambda = log(lam),
#                                             submodels = A_val,
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
# write.csv(overall_results_HILS_2, file ="./n14_k20_kstar5_allpos_designs_HILS/prob_comparison_HILS_2.csv", row.names =FALSE, col.names = FALSE)

# I_results=rbind(I_results_d1,I_results_d2, I_results_HILS, I_results_HILS_2)
I_results=rbind(I_results_d1, I_results_HILS)
# S_results=rbind(S_results_d1,S_results_d2, S_results_HILS, S_results_HILS_2)
S_results=rbind(S_results_d1, S_results_HILS)
joint_results = rbind(joint_results_d1, joint_results_HILS)
overall_results = I_results
overall_results$P_S=as.numeric(S_results$P_S)
overall_results$joint=as.numeric(joint_results$P_joint)
#overall_results$joint=as.numeric(I_results$P_I)*as.numeric(S_results$P_S)

write.csv(overall_results, file ="./n14_k20_kstar5_allpos_designs_HILS_lam_opt/prob_comparison_log_lam_opt.csv", row.names =FALSE, col.names = FALSE)


library(readr)
prob_comparison_d1 <- read_csv("prob_comparison_d1_log_lam_opt.csv")
# prob_comparison_d2 <- read_csv("prob_comparison_d2.csv")
prob_comparison_HILS <- read_csv("prob_comparison_HILS_log_lam_opt.csv")
# prob_comparison_HILS_2 <- read_csv("prob_comparison_HILS_2.csv")

# overall_results= rbind(prob_comparison_d1, prob_comparison_d2, prob_comparison_HILS, prob_comparison_HILS_2)
overall_results= rbind(prob_comparison_d1, prob_comparison_HILS)
#trimmed_overall= overall_results[overall_results$lambda<4.5,]
trimmed_overall=overall_results
trimmed_overall[trimmed_overall$design=="d1",]$design="PED"
trimmed_overall[trimmed_overall$design=="HILS",]$design="Var(s+)_Eff80"

trimmed_overall$transformed_joint = abs(1/(trimmed_overall$lambda))*trimmed_overall$joint
# 
# ggI = ggplot(data = overall_results[overall_results$design!="d2",], aes(x=as.numeric(lambda), y= as.numeric(P_I), color= design))+
#   geom_line(show.legend = FALSE)+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   ylab(expression(paste("P(",I,")")))+
#   scale_color_manual(values =c("blue","red", 'green', "orange"))+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55))
# 
# ggS = ggplot(data = overall_results[overall_results$design!="d2",], aes(x=as.numeric(lambda), y= as.numeric(P_S), color= design))+
#   geom_line(show.legend = FALSE)+
#   xlab(latex2exp::TeX("$\\lambda"))+
#   ylab(expression(paste("P(",italic(S),")")))+
#   #ylab(expression("P(\uD835\uDCAE)"))+
#   scale_color_manual(values =c("blue","red", 'green', "orange"))+
#   theme(axis.title.y=element_text(angle=0, vjust=0.55))
# 
ggjoint = ggplot(data = trimmed_overall[trimmed_overall$design!="d2" & trimmed_overall$design!="HILS 2" ,], aes(x=log(as.numeric(lambda)), y= as.numeric(joint), color= design, linetype = design))+
  geom_line()+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\phi_l(\\textbf{z}_A=\\textbf{1}_5, \\lambda)"))+
  scale_color_manual(values =c("blue","red"))+
  theme_bw()+
  scale_x_continuous(limits = c(-2,2))+
  scale_linetype_manual(values =c("solid", "dashed"))+
  theme( legend.title=element_blank(), text=element_text(size=16), legend.position = c(0.2,.75))
# 
# 
# ggall= ggarrange(ggjoint,ggarrange(ggI, ggS, ncol=2, labels=c("B","C")), nrow=2,labels="A")
# ggall
 ggsave("./Section4_n14_k20_joint_figure.png", dpi=600, device="png", ggjoint) ## this saves a 600 dpi image.
# 
# 
