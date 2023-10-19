########## 已经废弃

source("utils.R")
fun_linear_regression = function(D,B,level,e,r,BS_true,k_par){
  n = dim(D)[1]
  k = dim(D)[2]
  Dt = t(D)
  DtD = Dt%*%D
  n_total = n
  n_r = length(r)
  sensi_mean = max(D)-min(D)
  gamma = max(abs(D))/2
  
  theta_hat = apply(D,2,mean)
  beta_hat = theta_hat[1:k_par]
  cov_matrix = 1/(n-1)*DtD - n/(n-1)*theta_hat%*%t(theta_hat)
  BS_hat = max(beta_hat)
  s_max = which.max(beta_hat)
  
  Lap_noise = Lap_noise_Gaussian(k,e,k_par,theta_hat,sensi_mean,gamma)
  scl = Lap_noise$scale
  
  theta_hat_pri_whl = theta_hat + Lap_noise$whl_mean
  theta_hat_pri_par = theta_hat + Lap_noise$par_mean
  
  beta_hat_pri_whl = theta_hat_pri_whl[1:k_par]
  beta_hat_pri_par = theta_hat_pri_par[1:k_par]
  
  BS_hat_pri_whl = max(beta_hat_pri_whl)
  BS_hat_pri_par = max(beta_hat_pri_par)
  
  cov_matrix_whl = 1/(n-1)*(DtD+Lap_noise$whl_cov) - n/(n-1)*theta_hat_pri_whl%*%t(theta_hat_pri_whl)
  cov_matrix_par = 1/(n-1)*(DtD+Lap_noise$par_cov) - n/(n-1)*theta_hat_pri_par%*%t(theta_hat_pri_par)
  
  lb_naive = beta_hat[s_max] - qt(1-level,df = n_total-1)*sqrt(cov_matrix[s_max,s_max]/n)
  lb_naive_whl = beta_hat_pri_whl[s_max] - qt(1-level,df = n_total-1)*sqrt(cov_matrix_whl[s_max,s_max]/n)
  lb_naive_par = beta_hat_pri_par[s_max] - qt(1-level,df = n_total-1)*sqrt(cov_matrix_par[s_max,s_max]/n)
  
  d = matrix(0,n_r,k_par)
  d_whl = matrix(0,n_r,k_par)
  d_par = matrix(0,n_r,k_par)
  n_r = length(r)
  for (i in 1:n_r){
    d[i,] = (1-n_total^(r[i]-0.5))*(BS_hat-beta_hat)
    d_whl[i,] = (1-n_total^(r[i]-0.5))*(BS_hat_pri_whl-beta_hat_pri_whl)
    d_par[i,] = (1-n_total^(r[i]-0.5))*(BS_hat_pri_par-beta_hat_pri_par)
  }
  
  # parametric bootstrap
  T_b = matrix(0,n_r,B)
  T_b_whl = matrix(0,n_r,B)
  T_b_par = matrix(0,n_r,B)
  
  for (b in 1:B){
    D_boot = mvrnorm(n,rep(0,k), cov_matrix)
    D_boot_whl = mvrnorm(n,rep(0,k), cov_matrix_whl)
    D_boot_par = mvrnorm(n,rep(0,k), cov_matrix_par)
    
    theta_boot = apply(D_boot,2,mean)
    Lap_noise_b = Lap_noise_Gaussian(k,e,k_par,theta_boot,sensi_mean,gamma)
    
    theta_boot_whl = theta_boot + Lap_noise_b$whl_mean
    theta_boot_par = theta_boot + Lap_noise_b$par_mean

    beta_boot = theta_boot[1:k_par]
    beta_boot_whl = theta_boot_whl[1:k_par]
    beta_boot_par = theta_boot_par[1:k_par]
    for (i in 1:n_r){
      T_b[i,b] = sqrt(n_total)*(max(beta_boot+d[i,])-BS_hat)
      T_b_whl[i,b] = sqrt(n_total)*(max(beta_boot_whl+d_whl[i,])-BS_hat_pri_whl)
      T_b_par[i,b] = sqrt(n_total)*(max(beta_boot_par+d_par[i,])-BS_hat_pri_par)
    }
  }
  T_b = na.omit(T_b)
  T_b_whl = na.omit(T_b_whl)
  T_b_par = na.omit(T_b_par)
  lb = 0
  lb_whl = 0
  lb_par = 0
  for (i in 1:n_r){
    lb[i] = unname(BS_hat - quantile(T_b[i,],1-level)/sqrt(n_total)) # lower bound
    lb_whl[i] = unname(BS_hat_pri_whl - quantile(T_b_whl[i,],1-level)/sqrt(n_total)) # lower bound
    lb_par[i] = unname(BS_hat_pri_par - quantile(T_b_par[i,],1-level)/sqrt(n_total)) # lower bound
  }
  result <- list(lb,lb_whl,lb_par,lb_naive,lb_naive_whl,lb_naive_par,s_max,scl)
  names(result) <- c("nonpri_bstp","priwhole_bstp","pripartial_bstp","nonpri_naive","priwhole_naive","pripartial_naive","best_subgroup","scale_of_noise")
  return(result)
}

n = 400
k = 2
k_par = 2
B = 200
e = 500
level = 0.05
r = c(-Inf,1/30,1/15,1/10,1/5,0.5)
D = mvrnorm(n,rep(0,k), diag(k))





fun_mv_Gassian = function(D,B,level,e,r,BS_true,k_par){
  n = dim(D)[1]
  k = dim(D)[2]
  Dt = t(D)
  DtD = Dt%*%D
  n_total = n
  n_r = length(r)
  sensi_mean = max(D)-min(D)
  gamma = max(abs(D))/2
  
  theta_hat = apply(D,2,mean)
  beta_hat = theta_hat[1:k_par]
  cov_matrix = 1/(n-1)*DtD - n/(n-1)*theta_hat%*%t(theta_hat)
  BS_hat = max(beta_hat)
  s_max = which.max(beta_hat)
  
  Lap_noise = Lap_noise_mv_Gaussian(k,e,k_par,theta_hat,sensi_mean,gamma)
  scl = Lap_noise$scale
  
  theta_hat_pri_whl = theta_hat + Lap_noise$whl_mean
  theta_hat_pri_par = theta_hat + Lap_noise$par_mean
  
  beta_hat_pri_whl = theta_hat_pri_whl[1:k_par]
  beta_hat_pri_par = theta_hat_pri_par[1:k_par]
  
  BS_hat_pri_whl = max(beta_hat_pri_whl)
  BS_hat_pri_par = max(beta_hat_pri_par)
  
  cov_matrix_whl = 1/(n-1)*(DtD+Lap_noise$whl_cov) - n/(n-1)*theta_hat_pri_whl%*%t(theta_hat_pri_whl)
  cov_matrix_par = 1/(n-1)*(DtD+Lap_noise$par_cov) - n/(n-1)*theta_hat_pri_par%*%t(theta_hat_pri_par)
  
  lb_naive = beta_hat[s_max] - qt(1-level,df = n_total-1)*sqrt(cov_matrix[s_max,s_max]/n)
  lb_naive_whl = beta_hat_pri_whl[s_max] - qt(1-level,df = n_total-1)*sqrt(cov_matrix_whl[s_max,s_max]/n)
  lb_naive_par = beta_hat_pri_par[s_max] - qt(1-level,df = n_total-1)*sqrt(cov_matrix_par[s_max,s_max]/n)
  
  d = matrix(0,n_r,k_par)
  d_whl = matrix(0,n_r,k_par)
  d_par = matrix(0,n_r,k_par)
  n_r = length(r)
  for (i in 1:n_r){
    d[i,] = (1-n_total^(r[i]-0.5))*(BS_hat-beta_hat)
    d_whl[i,] = (1-n_total^(r[i]-0.5))*(BS_hat_pri_whl-beta_hat_pri_whl)
    d_par[i,] = (1-n_total^(r[i]-0.5))*(BS_hat_pri_par-beta_hat_pri_par)
  }
  
  # parametric bootstrap
  T_b = matrix(0,n_r,B)
  T_b_whl = matrix(0,n_r,B)
  T_b_par = matrix(0,n_r,B)
  
  for (b in 1:B){
    D_boot = mvrnorm(n,rep(0,k), cov_matrix)
    D_boot_whl = mvrnorm(n,rep(0,k), cov_matrix_whl)
    D_boot_par = mvrnorm(n,rep(0,k), cov_matrix_par)
    
    theta_boot = apply(D_boot,2,mean)
    Lap_noise_b = Lap_noise_mv_Gaussian(k,e,k_par,theta_boot,sensi_mean,gamma)
    
    theta_boot_whl = theta_boot + Lap_noise_b$whl_mean
    theta_boot_par = theta_boot + Lap_noise_b$par_mean
    
    beta_boot = theta_boot[1:k_par]
    beta_boot_whl = theta_boot_whl[1:k_par]
    beta_boot_par = theta_boot_par[1:k_par]
    for (i in 1:n_r){
      T_b[i,b] = sqrt(n_total)*(max(beta_boot+d[i,])-BS_hat)
      T_b_whl[i,b] = sqrt(n_total)*(max(beta_boot_whl+d_whl[i,])-BS_hat_pri_whl)
      T_b_par[i,b] = sqrt(n_total)*(max(beta_boot_par+d_par[i,])-BS_hat_pri_par)
    }
  }
  T_b = na.omit(T_b)
  T_b_whl = na.omit(T_b_whl)
  T_b_par = na.omit(T_b_par)
  lb = 0
  lb_whl = 0
  lb_par = 0
  for (i in 1:n_r){
    lb[i] = unname(BS_hat - quantile(T_b[i,],1-level)/sqrt(n_total)) # lower bound
    lb_whl[i] = unname(BS_hat_pri_whl - quantile(T_b_whl[i,],1-level)/sqrt(n_total)) # lower bound
    lb_par[i] = unname(BS_hat_pri_par - quantile(T_b_par[i,],1-level)/sqrt(n_total)) # lower bound
  }
  result <- list(lb,lb_whl,lb_par,lb_naive,lb_naive_whl,lb_naive_par,s_max,scl)
  names(result) <- c("nonpri_bstp","priwhole_bstp","pripartial_bstp","nonpri_naive","priwhole_naive","pripartial_naive","best_subgroup","scale_of_noise")
  return(result)
}
