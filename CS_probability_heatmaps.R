
# This script produces the contour plots for CS information matrices given in section 3 of the paper. 
# A simple change in parameters will give the contour plots in section 5




library(glmnet)
library(dplyr)
library(tidyr)
library(eigeninv)
library(Matrix)
library(pracma)
library(MASS)
library(ggplot2)
library(mvtnorm)
library(MixMatrix)
library(plotly)
library(ggpubr)
library(gridExtra)


root_dir = "PATH TO DIRECTORY CONTAINING CODE AND DESIGN CATALOG"

source(paste0(root_dir,'/Lasso_optimal_SSD_function_library.R'))

### The above path should point to whatever the path is in your local files

make_cs_matrix <-function(p,c){
  # This generates a completely symmetric matrix with off diagonals c
  return(CSgenerate(p,c))
}
prob_conditon_CS <- function(c, sigma_sq, n, kstar,q,lambda,B1, sign_vec){
  # This generates the probabilites of inclusion and exclusion for a CS X1TX1 and X2TX2
  
  
  X1TX1= n*make_cs_matrix(p=kstar,c=c)
  X2TX2 =n* make_cs_matrix(p=q,c=c)
  
  X2TX1 = n*c*matrix(rep(1,kstar*q), nrow=q, ncol=kstar)
  
  X1TX2=t(X2TX1)
  
  Sig1<-diag(1/(1-c),kstar)-(c/(1-c))*(1/(1+c*(kstar-1)))*matrix(1,nrow=kstar,ncol=kstar)
  Var_inc<-sigma_sq*diag(sign_vec)%*%Sig1%*%diag(sign_vec)
  #Var_inc = sigma_sq * n *(diag(sign_vec)%*% solve(X1TX1)%*% diag(sign_vec))
  
  
  one_vec_inc = rep(1,kstar)
  
  z = t(one_vec_inc)%*%sign_vec
  #browser()
  
  
  mean_inc = sqrt(n)*lambda*( ( (1/(1-c))*one_vec_inc ) - (z*(c/((1-c)*(1+c*(kstar-1))))*sign_vec) )
  
  
  #rhs_inc = n*lambda* solve(X1TX1)%*%one_vec_inc-B1
  
  #mean_vec = rhs_inc
  
  #browser()
  rhs_vec = sqrt(n) * one_vec_inc*B1
  prob_inc =pmvnorm(upper = c(rhs_vec), mean = c(mean_inc), sigma=Var_inc)
  
  #browser()
  
  I_q = diag(rep(1,q))
  
  J_q = matrix(rep(1,q*q), nrow=q, ncol=q)
  
  Var_ex = sigma_sq* ( ((1-c)*I_q) + ((c*(1-c)/(1+c*(kstar-1)))*J_q) )
  one_vec_ex =rep(1,q)
  
  mean_ex = sqrt(n)*lambda*z* (c/(1+c*(kstar-1)))*one_vec_ex
  
  #bound = n*lambda*(one_vec_ex-abs( ((c*pstar)/(1+c*(pstar-1))) *one_vec_ex))
  bound = sqrt(n)*lambda*one_vec_ex
  if (any(bound < 0)){
    prob_ex =0
  }
  else{prob_ex=pmvnorm(lower=-bound, upper=bound,mean = c(mean_ex), sigma=Var_ex)}
  return(list('prob_inclusion'=prob_inc[1], 'prob_exclusion' = prob_ex[1], "joint"= prob_inc[1]*prob_ex[1]))
}

get_c_list<- function(p){
  # this returns a list of evenly spaced c values based on the number of columns p
  c_start=(-0.99)/(p-1)
  c_end=0.99/(p-1)
  c_int = (c_end-c_start)/98
  c_list = seq(c_start, c_end,c_int)
  return(c_list)
}

c_list = seq(-1/10,0.98,length.out = 109)
#c_list[100]=0
h = 0.0001
kstar = 4
q=8
beta = 3
lambda_list= seq(0.501,3.501, by=0.03)
log_lambda_list = seq(-3,4, by = 0.06)
#log_lambda_list = log(seq(exp(-4.5), exp(2), by=exp(log(0.1))))
result_mat_inc = matrix(NA, nrow=109, ncol=length(log_lambda_list))
result_mat_ex=matrix(NA, nrow=109, ncol=length(log_lambda_list))
# result_mat_inc_der = matrix(NA, nrow=109, ncol=length(log_lambda_list))
# result_mat_ex_der = matrix(NA, nrow=109, ncol=length(log_lambda_list))
# result_mat_overall_der = matrix(NA, nrow=109, ncol=length(log_lambda_list))


result_df = data.frame(c=numeric(0), log_lambda=numeric(0), P_I=numeric(0), P_S= numeric(0), Joint=numeric(0))
n=10

sigma_sq=1

# for(c in c_list){
#   lambda = 0.1
#   condition_prob = prob_conditon_CS(c=c, sigma_sq=1, n=n, kstar=kstar,q=q,lambda=lambda,B1=3, sign_vec= rep(1,kstar))
#   #condition_prob = screen_prob_all_fixed(A= seq(1,kstar, by =1), n, C= make_cs_matrix(p=4, c=c),3, Z_A= rep(1,kstar), lambda)
#   #test_out[nrow(test_out)+1,]= c(c,condition_prob$prob_exclusion, condition_prob$prob_inclusion)
#   test_out[nrow(test_out)+1,]= c(c,condition_prob$P_I, condition_prob$P_S)}
k=1
for (j in 1:length(log_lambda_list)) {
  lambda= exp(log_lambda_list[j])
  #need to get valid c s.t. XTX is a PD matrix, not just X1TX1
  #c_list = get_c_list(kstar+q)
  print(j)
  #browser()
  for(i in 1:length(c_list)){
    c = c_list[i]
    condition_prob = prob_conditon_CS(c=c, sigma_sq=1, n=n, kstar=kstar,q=q,lambda=lambda,B1=beta, sign_vec= rep(1,kstar))
    #condition_prob_plus_h = prob_conditon_CS(c=c+h, sigma_sq=1, n=n, kstar=kstar,q=q,lambda=lambda,B1=beta, sign_vec= rep(1,kstar))
    #condition_prob_minus_h = prob_conditon_CS(c=c-h, sigma_sq=1, n=n, kstar=kstar,q=q,lambda=lambda,B1=beta, sign_vec= rep(1,kstar))
    
    result_mat_inc[i,j] <- condition_prob$prob_inclusion
    result_mat_ex[i,j]<-condition_prob$prob_exclusion
    # if (i==11 & j==25){
    #   print(condition_prob_plus_h$prob_exclusion)
    #   print(condition_prob_minus_h$prob_exclusion)
    #   print(condition_prob_plus_h$prob_exclusion-condition_prob$prob_exclusion)
    # }
    #result_mat_inc_der[i,j] = (condition_prob_plus_h$prob_inclusion - condition_prob_minus_h$prob_inclusion)/(2*h)
    #result_mat_ex_der[i,j] = (condition_prob_plus_h$prob_exclusion - condition_prob_minus_h$prob_exclusion)/(2*h)
    #result_mat_overall_der[i,j] = (condition_prob_plus_h$joint - condition_prob_minus_h$joint)/(2*h)
    
    # results_row_ex= c(c_list[i],log_lambda_list[j], condition_prob$prob_exclusion, condition_prob$prob_inclusion, condition_prob$joint)
    # result_df[k,] =results_row_ex
    # k = k+1
  }
  
}
combine_mat=result_mat_ex*result_mat_inc
c_opt = c_list[which(rowSums(combine_mat)==max(rowSums(combine_mat), na.rm=TRUE), arr.ind = TRUE)[1]]
print(c_opt)



# Plotting P(S)
fig1<-plot_ly(x=log_lambda_list, y=get_c_list(kstar+q), z=result_mat_inc, type='contour')


fig1<-plot_ly(x=log_lambda_list, y=c_list, z=result_mat_inc, type='contour')
fig1%>%layout(title= TeX("$P(S_{\\lambda}|F,\\beta)"),
              yaxis=list(title = "c", titlefont=list(size=25), tickfont= list(size =18), range=c(-1/9,0.99)),
              xaxis=list(title = TeX("log(\\lambda)"),tickfont= list(size =18), range=c(-2, 2))) %>% config(mathjax = "cdn")



# Plotting P(I)


fig2<-plot_ly(x=log_lambda_list, y=c_list, z=result_mat_ex, type='contour')
fig2%>%layout(title= TeX("$P(I_{\\lambda}|F,\\beta)"),
              yaxis=list(title = "c", titlefont=list(size=25), tickfont= list(size =18), range=c(-1/9,0.99)),
              xaxis=list(title = TeX("log(\\lambda)"),tickfont= list(size =18),tick0=-2, dtick=2, range=c(-2, 4))) %>% config(mathjax = "cdn")





# Plotting Joint

fig3<-plot_ly(x=log_lambda_list, y=c_list, z=result_mat_inc*result_mat_ex, type='contour',autocontour=T)%>%
  layout(title= TeX("$P(z|F,\\beta)"),
         yaxis=list(title = "c", titlefont=list(size=25), tickfont= list(size =22), range=c(-1/9,0.99)),
         xaxis=list(title = TeX("log(\\lambda)"),tickfont= list(size =22), dtick=2,range=c(-2, 4))) %>% config(mathjax = "cdn")
 


fig3%>%layout(yaxis=list(title = "c", font=list(size=15), tickfont=list(size =15),range=c(-1/12,0.99)),
              xaxis=list(title = expression(lambda), tickfont=list(size =15),range=c(-2, 4)))%>%
        layout(shapes=list(type="line", x0=-2, x1=4, y0=c_opt, y1=c_opt, line=list(color="red", width=7)))
              


combine_mat=result_mat_ex*result_mat_inc

# Get optimal value of c and lambda

c_opt = c_list[which(rowSums(combine_mat)==max(rowSums(combine_mat), na.rm=TRUE), arr.ind = TRUE)[1]]
print(c_opt)

lam_opt = lambda_list[which(combine_mat==max(combine_mat, na.rm=TRUE), arr.ind = TRUE)[2]]


#Now calculating across all signs

Z= as.matrix(gen_all_sign_vects(c(1,2,3,4)))

all_signs_inc = matrix(0, nrow=109, ncol=length(log_lambda_list))
all_signs_ex = matrix(0, nrow=109, ncol=length(log_lambda_list))
all_signs_joint = matrix(0, nrow=109, ncol=length(log_lambda_list))
result_df_all_signs = data.frame(c=numeric(0), log_lambda=numeric(0), Joint=numeric(0))

for(l in 1:8){
  sign_vect = Z[l,]
  result_mat_inc = matrix(NA, nrow=109, ncol=length(log_lambda_list))
  result_mat_ex=matrix(NA, nrow=109, ncol=length(log_lambda_list))
  result_mat_joint=matrix(NA, nrow=109, ncol=length(log_lambda_list))
  result_df_all_signs_fixed = data.frame(c=numeric(0), log_lambda=numeric(0), Joint=numeric(0))
  k=1
  for (j in 1:length(log_lambda_list)) {
    lambda= exp(log_lambda_list[j])
    #need to get valid c s.t. XTX is a PD matrix, not just X1TX1
    #c_list = get_c_list(kstar+q)
    print(j)
    #browser()
    for(i in 1:length(c_list)){
      c = c_list[i]
      condition_prob = prob_conditon_CS(c=c, sigma_sq=1, n=n, kstar=kstar,q=q,lambda=lambda,B1=2, sign_vec= sign_vect)
      result_mat_inc[i,j] <- condition_prob$prob_inclusion
      result_mat_ex[i,j]<-condition_prob$prob_exclusion
      result_mat_joint[i,j]<-condition_prob$joint
      results_row_ex= c(c_list[i],log_lambda_list[j], condition_prob$joint)
      result_df_all_signs_fixed[k,] =results_row_ex
      k = k+1
    }
    
  }
  all_signs_inc = all_signs_inc+result_mat_inc
  all_signs_ex=all_signs_ex+result_mat_ex
  all_signs_joint=all_signs_joint+result_mat_joint
  if (l==1){
    result_df_all_signs= result_df_all_signs_fixed
  }
  else{ result_df_all_signs$Joint = result_df_all_signs$Joint + result_df_all_signs_fixed$Joint}
  
}

mean_all_signs_inc = all_signs_inc/8
mean_all_signs_ex = all_signs_ex/8
combine_mat_all_signs = all_signs_joint/8
c_opt_all =  c_list[which(rowSums(combine_mat_all_signs)==max(rowSums(combine_mat_all_signs), na.rm=TRUE), arr.ind = TRUE)[1]]

print(c_opt_all)
result_df_all_signs_mean = result_df_all_signs
result_df_all_signs_mean$Joint = result_df_all_signs$joint/8




fig4<-plot_ly(x=log_lambda_list, y=c_list, z=combine_mat_all_signs, type='contour',autocontour=T)
fig4%>%layout(title= TeX("$P(Z|F,\\beta)"),
              yaxis=list(title = "c", titlefont=list(size=18), tickfont= list(size =22), range=c(-1/9,0.99)),
              xaxis=list(title = TeX("log(\\lambda)"),tickfont= list(size =22), tick0=-2, dtick=2,range=c(-2, 4))) %>% config(mathjax = "cdn")




# plot_ly(x=log_lambda_list, y=c_list, z=mean_all_signs_, type='contour',autocontour=T)%>%layout(title= TeX("$P(Z|F,\\beta)"),
#                                                                                                      yaxis=list(title = "c", titlefont=list(size=25), tickfont= list(size =18), range=c(-1/9,0.99)),
#                                                                                                     xaxis=list(title = TeX("log(\\lambda)"),tickfont= list(size =18), tick0=-2, dtick=2,range=c(-2, 4))) %>% config(mathjax = "cdn")
# 
