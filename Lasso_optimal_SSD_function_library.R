### This script houses all functions used for the Lasso Optimal SSD paper

library(mvtnorm)
library(dplyr)
library(tidyr)
library(MASS)
library(ibd)
library(doParallel)
library(foreach)
library(HadamardR)
library(readr)
library(ggpubr)

read_in_design<-function(path){
  df <- read_table2(path,
                    col_names = FALSE)
  
  df <- df[, 1:ncol(df)-1]
  return(df)
}


### Function for creating the diagonals of V_half, use lapply to find V
fn<-function(x){
  sqrt(sum((x-mean(x))^2)/length(x))
}
### Function for centering and scaling, use lapply to center and scale F
cent_scale <- function(x){
  (x - mean(x))/ (sqrt((sum((x-mean(x))^2))/length(x)))
}

get_cs_and_V <- function(mat){
  F_cs = as.data.frame(lapply(mat, cent_scale))
  V_half=diag(lapply(mat,fn))
  return(list("F_cs"= F_cs, "V_half"=V_half))
}



check_inverse <- function(M){
  class(try(solve(M), silent = T))=="matrix"
}

fix_neg_diag <- function(M){
  if (any(diag(M)<=0)){
    # When we have negative diags in covariance matrix for I, they are very close to zero.
    # So, using the absolute value will flip them to slightly positive
    diag(M)=abs(diag(M))+ 1e-5
  }
  return(M)
}


#'Joint_prob_all_fixed
#'
#'Calculates the joint probability of sign recovery for a given design when beta, A, Z_A and lambda are fixed
#'
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param Z_A Sign vector of active effects
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda scalar, log of tuning parameter
#'
#' @return P_joint, scalar probability of sign recovery
#' @export
#'
#' @examples
Joint_prob_all_fixed <-function(A, F_cs, V_half, B_mag, Z_A, log_lambda, sigma = 1){
  # Subsetting design into active and inactive columns
  print(A)
  F_A=as.matrix(F_cs[,A])
  F_I = as.matrix(F_cs[,-A])
  V_half=diag(diag(V_half)[A])

  lambda = exp(log_lambda)

  n = dim(F_cs)[1]
  k = length(Z_A)
  one_k = rep(1,k)
  one_q = rep(1, dim(F_I)[2])
  B_min_vec = B_mag*Z_A
  C_IA = (1/n)*t(F_I)%*%F_A
  C_II = (1/n)*t(F_I)%*%F_I
  C_AA = (1/n)*t(F_A)%*%F_A



  # Check if C_AA has an inverse, if not, use moore-penrose. 
    C_AA_inv = n*ginv(t(F_A)%*%F_A)
    Mean_I = lambda*sqrt(n)*C_IA%*%C_AA_inv%*%Z_A
    Cov_I = (sigma^2)*(C_II - (C_IA%*%C_AA_inv%*%t(C_IA)))

    # Fix the diagonals of Cov_I so that they are slightly positive, if they are slightly negative
    Cov_I = fix_neg_diag(Cov_I)
    P_I = mvtnorm::pmvnorm(lower = -sqrt(n)*lambda*one_q, upper = sqrt(n)*lambda*one_q,
                           mean = as.vector(Mean_I), sigma = as.matrix(Cov_I), algorithm = GenzBretz() )[1]
    Cov_S = (sigma^2)* diag(Z_A)%*%C_AA_inv%*%diag(Z_A)
    upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%B_min_vec -
      lambda*sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
    
    P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S), mean= rep(0,k), sigma = Cov_S)[1]
    P_joint = P_S*P_I

 return(P_joint)

}

## This function is a deconstructed version of Joint_prob_all_fixed for optimization/ integration over lambda
Joint_prob_all_fixed_deconstructed<-function(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S,mean_S_WO_lam, log_lambda){

  lambda = exp(log_lambda)
  k = dim(Cov_S)[1]


  Mean_I = mean_I_WO_lam*lambda
  

  P_I = mvtnorm::pmvnorm(lower = lambda*as.vector(lower_I_WO_lam), upper = lambda * as.vector(upper_I_WO_lam),
                         mean = as.vector(Mean_I), sigma = as.matrix(Cov_I), algorithm = GenzBretz() )[1]
  

  P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S)-lambda*as.vector(mean_S_WO_lam), mean= rep(0,k), sigma = Cov_S)[1]


  P_joint = P_S*P_I
  return(P_joint)

}


Joint_prob_all_fixed_deconstructed_more_output<-function(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S,mean_S_WO_lam, log_lambda, output_option="joint"){
  
  lambda = exp(log_lambda)
  k = dim(Cov_S)[1]
  
  
  Mean_I = mean_I_WO_lam*lambda
  
  P_I = mvtnorm::pmvnorm(lower = lambda*as.vector(lower_I_WO_lam), upper = lambda * as.vector(upper_I_WO_lam),
                         mean = as.vector(Mean_I), sigma = as.matrix(Cov_I), algorithm = GenzBretz() )[1]
  #browser()
  
  P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S)-lambda*as.vector(mean_S_WO_lam), mean= rep(0,k), sigma = Cov_S)[1]
  
  
  P_joint = P_S*P_I
  if (output_option =="joint"){ return(P_joint)}
  if (output_option =="I"){return(P_I)}
  else{ return(P_S)}
  
}

### Strategies to handle lambda


#' opt_log_lambda
#' 
#' 
#'
#' @param log_lam_min Scalar, minimum for grid of log lambdas for a warm start
#' @param log_lam_max Scalar, maximum for grid of log lambdas for a warm start
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param Z_A Sign vector of active effects
#' @param sigma Scalar, variance parameter, default is 1
#' @param stepsize Scalar, stepsize between grid points for warm start, default is 0.05
#'
#' @return list with one item denoting the optimal log lambda value, and another giving the P_joint value at the optimal log lambda
#' @export
#'
#' @examples
opt_log_lambda <- function( log_lam_min, log_lam_max, F_cs, V_half, B_mag, k, submodels, sign_vects, sigma =1, stepsize= 0.05){
  # Get warm start for lamnda values
  #browser()
  log_lambda_grid = seq(log_lam_min,log_lam_max,stepsize)
  grid_probs <-sapply(log_lambda_grid, Joint_prob_general_fixed_lambda, F_cs=F_cs, V_half = V_half, B_mag=B_mag, k=k, submodels=submodels, sign_vects=sign_vects,
                      log_lambda_strategy="fixed", sigma=sigma, output_option="joint")
  warm_start <- cbind(log_lambda_grid,grid_probs)
  #print(max(warm_start[,2]))
  #browser()

  if (is.na(max(warm_start[,2]))|| max(warm_start[,2])<0.0001){
    start_lam <- runif(1, 0, 3)
  }
  else{
    start_lam <- warm_start[which.max(warm_start[,2]),1]
  }
  # maybe consider doing box-constaints, but since log-lam can take any real number, I thought other solvers could be faster
  #browser()
  opt_lam <-optim(par=start_lam, fn=Joint_prob_general_fixed_lambda, F_cs=F_cs, V_half = V_half, B_mag=B_mag, k=k, submodels=submodels, sign_vects=sign_vects,
                           log_lambda_strategy="fixed", sigma=sigma, output_option="joint",
                           method="Brent",lower=-5, upper =5,
                           control = list("fnscale"=-1,warn.1d.NelderMead=FALSE))
  # opt_lam <-tryCatch(optim(par=start_lam, fn=Joint_prob_all_fixed_deconstructed,lower_I_WO_lam=lower_I_WO_lam,
  #                 upper_I_WO_lam=upper_I_WO_lam, mean_I_WO_lam=mean_I_WO_lam,
  #                 Cov_I=Cov_I, Cov_S = Cov_S, upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam,
  #                 method="L-BFGS-B", upper = 10, lower = -Inf,
  #                 control = list("fnscale"=-1,warn.1d.NelderMead=FALSE)), error = function(c){return(list("par"=start_lam, "value"=0))})
  return(list("opt_log_lam"= opt_lam$par, "opt_val" = opt_lam$value))


}



#' This gives the area under the joint probability curve as it depends on log lambda
#'
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param Z_A Sign vector of active effects
#' @param sigma Scalar, variance parameter, default is 1
#' @param int_lower_bound Scalar, lower bound on region of integration, default is infinity
#' @param int_upper_bound Scalar, upper bound on region of integration, default is 3
#' @param method String,  Quadrature or Sum, Sum indicates Riemann sum
#'
#'
#' @return resulting value of integral
#' @export
#'
#' @examples
integrate_log_lambda <- function(lower_I_WO_lam, upper_I_WO_lam,
                                 mean_I_WO_lam, Cov_I, Cov_S, upper_S, mean_S_WO_lam, int_lower_bound=-5, int_upper_bound = 1.5, method = "Quadrature", step_size = 0.02){
 

  
  int_methods = c("Quadrature", "Sum")
  if(!method%in% int_methods){
    stop("The input integration method must be either  'Quadrature', or 'Sum'")
  }
  if (method == "Quadrature"){
    Joint_prob_all_fixed_vec = Vectorize(Joint_prob_all_fixed_deconstructed, "log_lambda")
    #browser()
    #print("Integrating via quadrature")
    
    integrate_results = integrate(Joint_prob_all_fixed_vec, lower= int_lower_bound, upper= int_upper_bound,
                                  lower_I_WO_lam=lower_I_WO_lam,
                                  upper_I_WO_lam=upper_I_WO_lam,
                                  mean_I_WO_lam=mean_I_WO_lam,
                                  Cov_I=Cov_I, Cov_S = Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam)
    return(integrate_results$value)
  }
  else{
    #print("Integrating via Riemann Sum")
    log_lambda_grid <- seq(int_lower_bound, int_upper_bound, by=step_size)
  
    
    above_0_yet = FALSE
    sum_prob = Joint_prob_all_fixed_deconstructed(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S, mean_S_WO_lam, log_lambda = int_lower_bound)
    for( j in c(1:1000)){
      log_lambda = int_lower_bound + step_size*j
      #print(log_lambda)
      prob = tryCatch({Joint_prob_all_fixed_deconstructed(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S, mean_S_WO_lam, log_lambda)}, error =function(c){return(0)})
      sum_prob = sum_prob + prob
      if (is.na(prob)){
        #something with MVN algorithm givens NaNs sometimes
        prob =0
      }
     
      #print(paste0("prob: ", prob))
      if( prob > 0.0001 && above_0_yet == FALSE){ 
        #print("prob is above 0 for first time")
        above_0_yet= TRUE}
      if( prob < 0.0001 && above_0_yet == TRUE){ 
        #print("breaking loop")
        break}
      
    }
    return(step_size* sum_prob)
  }
  


}


#'This is a function for a fixed submodel and sign vector, with an option for handling lambda
#'
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param Z_A Sign vector of active effects
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda Scalar, fixed value of log_lambda, this is set to NULL, if the log_lambda strategy is fixed, this must be defined
#' @param log_lambda_strategy string, strategy to handle log_lambdas, takes in either "fixed", "optimal", "integrate"
#' @param int_lower_bound
#' @param int_upper_bound
#'
#' @return P_joint, scalar probability of sign recovery
#' @export
#'
#' @examples
Joint_prob_submodel_sign_fixed <-function(A, F_cs, V_half, B_mag, Z_A, log_lambda=NULL, log_lambda_strategy="integrate",
                                          sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method = "Quadrature", output_option = "joint"){
  strategies = c("fixed", "optimal", "integrate")
  if(!log_lambda_strategy%in% strategies){
    stop("The input log_lam strategy must be either 'fixed, 'optimal', or 'integrate'")
  }
  # Subsetting design into active and inactive columns
  F_A=as.matrix(F_cs[,A])
  F_I = as.matrix(F_cs[,-A])
  V_half=diag(diag(V_half)[A])

  #lambda = exp(log_lambda)

  n = dim(F_cs)[1]
  k = length(Z_A)
  one_k = rep(1,k)
  one_q = rep(1, dim(F_I)[2])
  B_min_vec = B_mag*Z_A
  C_IA = (1/n)*t(F_I)%*%F_A
  C_II = (1/n)*t(F_I)%*%F_I
  C_AA = (1/n)*t(F_A)%*%F_A



  # Check if C_AA has an inverse, if not, set the joint prob to zero. (Discuss with group)
  
  C_AA_inv = n*ginv(t(F_A)%*%F_A)
  
  mean_I_WO_lam = sqrt(n)*C_IA%*%C_AA_inv%*%Z_A
  Cov_I = (sigma^2)*(C_II - (C_IA%*%C_AA_inv%*%t(C_IA)))
  upper_I_WO_lam = sqrt(n)*one_q
  lower_I_WO_lam = -upper_I_WO_lam
  
  # Fix the diagonals of Cov_I so that they are slightly positive, if they are slightly negative
  
  Cov_I = fix_neg_diag(Cov_I)
  Cov_S = (sigma^2)* diag(Z_A)%*%C_AA_inv%*%diag(Z_A)
  #browser()
  upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%B_min_vec
  mean_S_WO_lam = sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
  if( log_lambda_strategy=="fixed"){
    if(is.null(log_lambda)){
      stop("For a fixed lambda strategy, you must specify what the fixed log_lambda is. It is currently NULL ")
    }
    else{
      return(Joint_prob_all_fixed_deconstructed_more_output(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S, mean_S_WO_lam, log_lambda, output_option))
    }
  }
  if( log_lambda_strategy=="integrate"){
    if(any(is.null(c(int_lower_bound, int_lower_bound)))){
      warning(" Integral upper and/or lower bounds are not specified, using defualt of lower bound at -5.")
      return(integrate_log_lambda(lower_I_WO_lam, upper_I_WO_lam,
                                  mean_I_WO_lam, Cov_I, Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam, method =int_method ))
    }
    else{ 
      return(integrate_log_lambda(lower_I_WO_lam, upper_I_WO_lam,
                                  mean_I_WO_lam, Cov_I, Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam,
                                  int_lower_bound=int_lower_bound ,
                                  int_upper_bound = int_upper_bound, method= int_method))
    }
    
    }

    if(log_lambda_strategy=="optimal"){
      #print(A)
      #browser()
      return(opt_log_lambda( log_lam_min=-2, log_lam_max=3, lower_I_WO_lam, upper_I_WO_lam,
                                         mean_I_WO_lam, Cov_I, Cov_S, upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam, stepsize= 0.05)$opt_val)
    }

  

}

### Now we move to a fixed submodel, but the ability to loop over sign vectors


gen_all_sign_vects<-function(A){
  k = length(A)
  l=rep(list(c(-1,1)), as.integer(k))
  # expand.grid will give a df where each row is a unique sign vector. The first half of the rows are the negartives of the second half
  # a is the index of the half way point of the rows
  #b is the largest row index
  a = ((2^(k-1))+1)
  b = 2^k
  sign_vec_grid = expand.grid(l)[a:b,]
  return(sign_vec_grid)
}

#'This is a function for a fixed submodel with the option to specify sign vectors and lambda strategy
#'
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param sign_vects Matrix, with |A| columns where each row is a sign vector, if NULL then it loops over all sign vectors. Defualt is NULL.
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda Scalar, fixed value of log_lambda, this is set to NULL, if the log_lambda strategy is fixed, this must be defined
#' @param log_lambda_strategy string, strategy to handle log_lambdas, takes in either "fixed", "optimal", "integrate"
#' @param int_lower_bound Scalar, lower bound in log scale for integral, default is NULL. 
#' @param int_upper_bound Scalar, lower bound in log scale for integral, default is NULL. 
#'
#' @return Dataframe where the first |A| columns represent the sign vector and the last column is the joint probablity
#' @export
#'
#' @examples
Joint_prob_submodel_fixed <-function(A, F_cs, V_half, B_mag, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                          sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Quadrature", output_option = "joint"){
  #print(A)
  #browser()
  if(is.null(sign_vects)){
    sign_vects= gen_all_sign_vects(A)
  }
  colnames(sign_vects) <- paste("s", 1:length(A),sep="")
  sign_vects = as.data.frame(sign_vects)
  #browser()
  results <- apply(as.data.frame(sign_vects), 1, Joint_prob_submodel_sign_fixed, A=A, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
                   log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
                   sigma = sigma, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound, int_method=int_method, output_option=output_option)
  #browser()
  output = cbind(sign_vects, results)
  return(output)
}


### Now we move to the completely general case, looping over all (or a given subset of) submodels and sign vectors


#This function gives a set of all p choose k submodels of size k. 
gen_all_submodels<-function(p,k){
  return(t(combn(p,k)))
}

# This function samples the submodels using the NBIBD approach. s_1 is the sample size. 
sample_submodels <-function(p, k, s_1){
  
  # this Uses NBIBD approach
  return(ibd(v=p, b= s_1, k=k)$design)
}

#' Stack_NBIBD
#' 
#' Function that permutes the labels and stacks NBIBD submodel sampling sets to make a larger submodel sampling set.
#'
#' @param p number of main effect factors
#' @param k number of factors in active set
#' @param s_1 The number for one NBIBD design, recommend no more than 64
#' @param num_stacks Intieger, number of times the NBIBD design of size s_1 will be permuted and stacked to give the final sample submodel set
#'
#' @return list: A_final is the final set of sampled submodels of size s_1*num_stacks, A_NBIBD is the original NBIBD samples submodels of size s_1. 
#' @export
#'
#' @examples
stack_NBIBD <-function(p,k,s_1,num_stacks){
  A_NBIBD = sample_submodels(p=p,k=k,s_1=s_1)
  A_final = A_NBIBD
  for(s in 1:(num_stacks-1)){
    A_perm = A_NBIBD
    mapping = sample(1:p,size =p)
    for (i in 1:s_1){
      for (j in 1:k){
        A_perm[i,j]= mapping[A_NBIBD[i,j]]
        
      }
    }
    A_final = rbind(A_final, A_perm)
  }
  
  return(list("A_final"= A_final, "A_NBIBD"=A_NBIBD))
}

#'This is a function for that gives the summary measure (mean only) over a set of submodels and sign vectors for a fixed lambda for lambda optimization
#'
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param k Integer, number of active effects. If NULL, it is the floor of n/3
#' @param sign_vects Matrix, with |A| columns where each row is a sign vector, if NULL then it loops over all sign vectors. Defualt is NULL.
#' @param sigma Scalar, variance parameter, default is 1
#' @param submodels Matrix, with k columns where each row represents a particular submodel column index, if NULL then it loops over all p choose k submodels. Defualt is NULL.
#' @param log_lambda Scalar, fixed value of log_lambda, this is set to NULL, this must be defined
#'
#' @return scalar, mean across submodels and sign vectors for a given lambda
#' @export
#'
#' @examples

Joint_prob_general_fixed_lambda <-function( F_cs, V_half, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=1,
                               sigma = 1,  output_option = "joint"){
  log_lambda_strategy = "fixed"
  p=dim(F_cs)[2]
  n = dim(F_cs)[1]
  if(any(is.na(F_cs))){
    return(c(0,0,0)) 
  }
  if(is.null(submodels)){
    if(is.null(k)){
      k = floor(n/3)
    }
    submodels= gen_all_submodels(p,k)
  }
  else{ k = dim(submodels)[2]}
  if(is.null(sign_vects)){
    num_sign_vects = 2^(k-1)
  }
  else{
    num_sign_vects = dim(sign_vects)[1]
  }
  
  colnames(submodels) <- paste("A", 1:k,sep="")
  submodels.expanded = submodels[rep(seq_len(nrow(submodels)), rep(num_sign_vects,nrow(submodels))), 1:k]
  results <- do.call(rbind,apply(submodels,1, Joint_prob_submodel_fixed,sign_vects=sign_vects, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
                                 log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
                                 sigma = sigma, output_option=output_option))
  
  output = cbind(submodels.expanded, results)
  return(mean(output$results))
}


# This function wraps Joint_prob_general_fixed_lambda for optimizing over lambda

Joint_prob_general_optimal_lambda <-function( F_cs, V_half, B_mag,k=NULL, submodels=NULL, sign_vects=NULL,
                                            sigma = 1,  output_option = "joint", log_lam_min=-2, log_lam_max= 2, stepsize=0.1){
  p=dim(F_cs)[2]
  n = dim(F_cs)[1]
  if(any(is.na(F_cs))){
    return(c(0,0,0)) 
  }
  if(is.null(submodels)){
    if(is.null(k)){
      k = floor(n/3)
    }
    submodels= gen_all_submodels(p,k)
  }
  else{ k = dim(submodels)[2]}
  if(is.null(sign_vects)){
    num_sign_vects = 2^(k-1)
    sign_vects = gen_all_sign_vects(seq(1,k, by =1))
  }
  else{
    num_sign_vects = dim(sign_vects)[1]
  }
 
  results = opt_log_lambda(log_lam_min = log_lam_min, log_lam_max = log_lam_max, F_cs, V_half, k, B_mag=B_mag, submodels, sign_vects,sigma,stepsize)
  return(results$opt_val)
}


#'This is a function that gives the joint probability, or the intergral across log lambda, for all submodels and sign vectors supplied by the user. 
#'
#'
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param k Integer, size of submodels, only used if submodels=NULL, if k=NULL, it will be set to floor of n/3.
#' @param submodels Matrix, with k columns where each row represents a particular submodel column index, if NULL then it loops over all p choose k submodels. Defualt is NULL.
#' @param sign_vects Matrix, with k columns where each row is a sign vector, if NULL then it loops over all sign vectors. Defualt is NULL.
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda Scalar, fixed value of log_lambda, this is set to NULL, if the log_lambda strategy is fixed, this must be defined
#' @param log_lambda_strategy string, strategy to handle log_lambdas, takes in either "fixed", "optimal", "integrate". Note, do not use "optimal" for this arguement, use Joint_prob_general_optimal_lambda instead
#' @param int_lower_bound Scalar, log lambda lower bound of integration
#' @param int_upper_bound Scalar, log lambda upper bound of integration
#' @param int_method string, "Sum" or "Quadrature", Quadrature takes a long time
#' @param output_option string, must be one of "I","S",or "joint", specifies which event, or joint event the output is in terms of. 
#'
#' @return Dataframe where the first |A| columns represent the submodel index, the next |A| column represent the sign vector and the last column is the output_option probablity
#' @export
#'
#' @examples
Joint_prob_general <-function( F_cs, V_half, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                     sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", output_option = "joint"){
  p=dim(F_cs)[2]
  n = dim(F_cs)[1]
  if(any(is.na(F_cs))){
    return(c(0,0,0)) 
  }
  if(is.null(submodels)){
    if(is.null(k)){
      k = floor(n/3)
    }
    submodels= gen_all_submodels(p,k)
  }
  else{ k = dim(submodels)[2]}
  if(is.null(sign_vects)){
    num_sign_vects = 2^(k-1)
  }
  else{
    num_sign_vects = dim(sign_vects)[1]
  }
  
  colnames(submodels) <- paste("A", 1:k,sep="")
  submodels.expanded = submodels[rep(seq_len(nrow(submodels)), rep(num_sign_vects,nrow(submodels))), 1:k]
  results <- do.call(rbind,apply(submodels,1, Joint_prob_submodel_fixed,sign_vects=sign_vects, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
                   log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
                   sigma = sigma, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound, int_method=int_method, output_option=output_option))
  
  output = cbind(submodels.expanded, results)
  return(output)
}

### Function comparing the Joint_prob_general to Var(s) and UE(S^2) measures for a set of designs 



### Measures of output of Joint_prob_general

#' Title
#'
#' @param results vector of the results for the set of submodels and sign vectors. This is the last column of the output from Joint_prob_general.
#'
#' @return List of several summary measures.
#' @export
#'
#' @examples
measure_list <-function(results){
  percent_0 = sum(results ==0)/length(results)
  percent_LT_0.1 = sum(results <0.1)/length(results)
  percent_GT_0.8 = sum(results >0.8)/length(results)
  summary_results = summary(results)
  percentile_25=as.numeric(summary_results[2])
  percentile_75=as.numeric(summary_results[5])
  mean_results = as.numeric(summary_results[4])
  median_results = as.numeric(summary_results[3])
  return(list("percent_0"= percent_0, "percent_LT_0.1"= percent_LT_0.1,
              "percent_GT_0.8"=percent_GT_0.8, "percentile_25"=percentile_25, 
              "percentile_75"=percentile_75, "mean"=mean_results,
              "median"=median_results))
}





### Hueristics and Hueristic design construction

UES <- function(F_0){
  p = dim(F_0)[2]
  n = dim(F_0)[1]
  F_0_w_int =as.data.frame(F_0)
  F_0_w_int$int = 1
  
  FTF= t(as.matrix(F_0_w_int))%*%as.matrix(F_0_w_int)
  FTF_upper <- FTF[upper.tri(FTF, diag = FALSE)]
  UES <- (2/(p*(p+1))) * sum(FTF_upper)
  return(UES)
}

UES_squared <- function(F_0){
  p = dim(F_0)[2]
  n = dim(F_0)[1]
  F_0_w_int =as.data.frame(F_0)
  F_0_w_int$int = 1
  
  FTF= t(as.matrix(F_0_w_int))%*%as.matrix(F_0_w_int)
  FTF_upper <- FTF[upper.tri(FTF, diag = FALSE)]
  UES_sq = (2/(p*(p+1))) * sum(FTF_upper^2)
  return(UES_sq)
}

#' Title
#'
#' @param F_0 matrix: matrix to evaluate
#' @param UES_sq_star scalar, optimal UES_sq value for a design with dimensions matching F_0 
#' @param eff scalar, between 0 and 1. Efficiency constraint
#' @param UES_const scalar, lower bound on UES value. Default is 0
#'
#' @return list of var_s value, UES value, and "valid": true or false value on if the design meets the efficiency and UES constraints
#' @export 
#'
#' @examples
Var_S_Plus <- function(F_0, UES_sq_star, eff = 0.8, UES_const = 0){
  p = dim(F_0)[2]
  n = dim(F_0)[1]
  F_0_w_int =as.data.frame(F_0)
  F_0_w_int$int = 1
  
  FTF= t(as.matrix(F_0_w_int))%*%as.matrix(F_0_w_int)
  FTF_upper <- FTF[upper.tri(FTF, diag = FALSE)]
  UES_sq = (2/(p*(p+1))) * sum(FTF_upper^2)
  UES <- (2/(p*(p+1))) * sum(FTF_upper)
  #browser()
  if(UES >UES_const & (UES_sq_star/UES_sq) > eff){
    valid_result <- TRUE
  }
  else{ valid_result = FALSE}
  Var_s = UES_sq - (UES^2)
  return(list("Var_s"=Var_s, "valid"= valid_result, "UES"=UES))
}


## Based on the Jones Paper

UES_sq_construction<-function(n,main_effects){
  p = main_effects+1
  if(p%%4 ==0){
    H_p = Normcol(Hadamard_Matrix(p))
    X_0 =  H_p[sample(p, size = n, replace = FALSE),]
  }
  if (p%%4 ==1){
    H = Normcol(Hadamard_Matrix(p-1))
    V = H[sample(p-1, size = n, replace = FALSE),]
    phi = sample(c(-1,1), n, replace=TRUE)
    X_0 = as.matrix(cbind(V, phi))
    
  }
  if (p%%4 ==2){
    if(n%%2==0){
      m = n/2
      H = Normcol(Hadamard_Matrix(p-2))
      X_star = H[sample(p-2, size = n, replace = FALSE),]
      G = matrix(c(1,1,-1,-1,-1,1,1,-1), nrow = 4, ncol=2, byrow = TRUE)
      U_11 = G[sample(2,m, replace = TRUE),]
      U_12 = G[sample(3:4,m, replace = TRUE),]
      U_1 = rbind(U_11, U_12)
      X_0 = as.matrix(cbind(X_star, U_1))
    }
    if(n%%2!=0){
      m = (n-1)/2
      H = Normcol(Hadamard_Matrix(p-2))
      X_star = H[sample(p-2, size = n, replace = FALSE),]
      G = matrix(c(1,1,-1,-1,-1,1,1,-1), nrow = 4, ncol=2, byrow = TRUE)
      U_21 = G[sample(2,m, replace = TRUE),]
      U_22 = G[sample(3:4,m+1, replace = TRUE),]
      U_2 = rbind(U_21, U_22)
      X_0 = as.matrix(cbind(X_star, U_2))
       }
    
  }
  if (p%%4 ==3){
    H = t(Hadamard_Matrix(p+1))
    X_star = H[sample(p+1, size = n, replace = FALSE),]
    X_0 = X_star[,1:p]
  }
  return(X_0)
}


# Function to find minimum UES_sq value for a given design dimensions
min_UES_sq_val <- function(n, main_effects){
  p = main_effects+1
  if(p%%2 ==1){
    min_UES_sq = ((n*(n-1)) + (n*p*(p-n)))/(p*(p-1))
  }
  if (p%%4 ==0){
    min_UES_sq = (n*(p-n))/(p-1)
    
  }
  if (p%%4 ==2){
    if(n%%2==0){
      min_UES_sq = ((2*n*(n-2)) + (n*p*(p-n)))/(p*(p-1))
    }
    if(n%%2!=0){
      min_UES_sq = ((2*((n-1)**2)) + (n*p*(p-n)))/(p*(p-1))
    }
    
  }
  return(min_UES_sq)
}

coord_exchange_UES_sq<-function( F_init, max_iter=10000){
  
  F_0 = F_init
  n = dim(F_0)[1]
  p = dim(F_0)[2]
  
  
  P_0 =  as.numeric(UES_squared(F_0))
  
  
  iterations = 0
  # Intitialize the count between finding a design with a higher measure.
  count_between_flips=0
  # exchange_data = data.frame( exchange_ind = numeric(0), sum_correl=numeric(0), max_correl=numeric(0), 
  #                             min_correl= numeric(0), improvement=numeric(0))
  stop = FALSE
  while ( iterations <= max_iter ){
    #print(paste0("Iteration: ", iterations))
    for (i in c(1:n)){
      for (j in c(1:p)){
        # print(i)
        # print(j)
        #print(paste0("Count Between Flips: ", count_between_flips))
        F_1 = F_0
        # Change Coordinate and update scaling
        F_1[i,j]= -F_0[i,j]
        
        P_1 = as.numeric(UES_squared(F_1))
        if( P_1 < P_0) {
          #Capture the sum and max absolute column wise correlations of the factor that is flipped
          #browser()
          # abs_correl_mat = (abs(t(as.matrix(F_0_cs))%*%as.matrix(F_0_cs))/n)- diag(1, nrow=p, ncol=p)
          # sum_correl = rowSums(abs_correl_mat)[j]
          # max_correl = max(abs_correl_mat[j,])
          # min_correl = min(abs_correl_mat[j,])
          # imp= P_1-P_0
          # exchange_data[nrow(exchange_data)+1,] <- c(1, sum_correl,max_correl , min_correl, imp)
          #browser()
          F_0 = F_1
          P_0 = P_1
          # We just flipped designs so now we reset count
          count_between_flips = 0
          
        }
        else{
          #Capture the sum and max absolute column wise correlations of the factor that is NOT flipped
          count_between_flips = count_between_flips+1
        }
        if(count_between_flips> n*p){
          # If we have looped through the whole matrix without changing designs break. 
          stop = TRUE
          break
        }
        #Breaking out of outter loops
        
      }
      if(stop){break}
      
    }
    if(stop){break}
    iterations=iterations+1
  }
  
  return(list("design"= F_0, "P_0"=P_0))
}


# Var(s+) coordinate exchange code
coord_exchange_Var_s<-function( n, p, eff = 0.8,UES_const = 0, max_iter=10000){
  UES_sq_standard = coord_exchange_UES_sq(matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p), max_iter = 1000)
  UES_sq_star = UES_sq_standard$P_0
  #print(UES_sq_star)
  
  
  invalid = TRUE
  # First exchange to get a valid var(s +) design
  while(invalid){
    # Indicator for negative average off diagonals
    neg_off <- TRUE
    #browser()
    # Start with random matrices until we get one with UES >0
    while( neg_off){
      
      F_init = matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p)
      
      if (UES(F_init)>UES_const){
        neg_off = FALSE
        #print("found positive matrix")
      }
      else{ neg_off = TRUE}
    }
    F_0 = F_init
    # Now we exchange to get up to the efficiency limit
    eff_0 = UES_sq_star / UES_squared(F_0)
    row = 1
    col = 1
    count = 0
    overall_count = 1000
    #browser()
    while( eff_0 < eff){
      F_1= F_0
      #exchange coordinates
      F_1[row,col] = -F_1[row,col]
      eff_1 = UES_sq_star / UES_squared(F_1)
      if(eff_1>eff_0){
        if( UES(F_1)>UES_const){
          F_0 = F_1
          eff_0 = eff_1
          count = 0
        }
      }
      row = (row + 1)%%n
      col = (col+ 1)%%p
      count = count +1
      if(count >n*p){
        #print("Not efficient enough, trying again.")
        #browser()
        break}
    }
    if(eff_0 >= eff){
      invalid = FALSE
    }
    
  }
  
  P_0 = as.numeric(Var_S_Plus(F_0, UES_sq_star, eff,UES_const = UES_const)$Var_s)
  
  #browser()
  
  
  iterations = 0
  # Intitialize the count between finding a design with a higher measure.
  count_between_flips=0
  # exchange_data = data.frame( exchange_ind = numeric(0), sum_correl=numeric(0), max_correl=numeric(0), 
  #                             min_correl= numeric(0), improvement=numeric(0))
  stop = FALSE
  while ( iterations <= max_iter ){
    #print(paste0("Iteration: ", iterations))
    for (i in c(1:n)){
      for (j in c(1:p)){
        # print(i)
        # print(j)
        #print(paste0("Count Between Flips: ", count_between_flips))
        F_1 = F_0
        # Change Coordinate and update scaling
        F_1[i,j]= -F_0[i,j]
        
        Var_s_results = Var_S_Plus(F_1, UES_sq_star, eff)
        if(Var_s_results$valid){
          P_1 = Var_s_results$Var_s
        }
        else{P_1 = P_0 +0.5}
        if( P_1 <= P_0) {
          #Capture the sum and max absolute column wise correlations of the factor that is flipped
          #browser()
          # abs_correl_mat = (abs(t(as.matrix(F_0_cs))%*%as.matrix(F_0_cs))/n)- diag(1, nrow=p, ncol=p)
          # sum_correl = rowSums(abs_correl_mat)[j]
          # max_correl = max(abs_correl_mat[j,])
          # min_correl = min(abs_correl_mat[j,])
          # imp= P_1-P_0
          # exchange_data[nrow(exchange_data)+1,] <- c(1, sum_correl,max_correl , min_correl, imp)
          #browser()
          F_0 = F_1
          P_0 = P_1
          # We just flipped designs so now we reset count
          count_between_flips = 0
          
        }
        else{
          #Capture the sum and max absolute column wise correlations of the factor that is NOT flipped
          count_between_flips = count_between_flips+1
        }
        if(count_between_flips> n*p){
          # If we have looped through the whole matrix without changing designs break. 
          stop = TRUE
          break
        }
        #Breaking out of outter loops
        
      }
      if(stop){break}
      
    }
    if(stop){break}
    iterations=iterations+1
  }
  #print(F_init - F_0)
  return(list("design"= F_0, "P_0"=P_0, "UES"= UES(F_0)))
}







