source("utils.R")
############################ Gaussian ############################
fun_Gassian = function(D,B,level,e,r,BS_true,k_par){
  n = dim(D)[1]
  k = dim(D)[2]
  n_total = n
  n_d = n*k
  n_r = length(r)
  
  theta_hat = apply(D,2,mean)
  #beta_hat = theta_hat[1:k_par]
  beta_hat = theta_hat
  
  sensi_mean = 0
  sensi_sd = 0
  sd_hat = 0
  for (m in 1:k){
    sensi_mean[m] = 1/n*abs(max(D[,m])-min(D[,m]))
    sensi_sd[m] = 2/n*abs(max(D[,m])-min(D[,m]))
    sd_hat[m] = sqrt(1/(n-1)*sum((D[,m]-theta_hat[m])^2))
  }
  BS_hat = max(beta_hat)
  s_max = which.max(beta_hat)
  
  Lap_noise = Lap_noise_Gaussian(k,e,k_par,theta_hat,sensi_mean,sensi_sd)
  scl = Lap_noise$scale
  
  theta_hat_pri_whl = theta_hat + Lap_noise$whl_mean
  theta_hat_pri_par = theta_hat + Lap_noise$par_mean
  
  # beta_hat_pri_whl = theta_hat_pri_whl[1:k_par]
  # beta_hat_pri_par = theta_hat_pri_par[1:k_par]
  beta_hat_pri_whl = theta_hat_pri_whl
  beta_hat_pri_par = theta_hat_pri_par
  
  BS_hat_pri_whl = max(beta_hat_pri_whl)
  BS_hat_pri_par = max(beta_hat_pri_par)
  
  sd_hat_pri_whl = sd_hat + Lap_noise$whl_cov
  sd_hat_pri_par = sd_hat + Lap_noise$par_cov
  
  lb_naive = beta_hat[s_max] - qt(1-level,df = n_total-1)*sd_hat[s_max]/sqrt(n)
  lb_naive_whl = beta_hat_pri_whl[s_max] - qt(1-level,df = n_total-1)*sd_hat_pri_whl[s_max]/sqrt(n)
  lb_naive_par = beta_hat_pri_par[s_max] - qt(1-level,df = n_total-1)*sd_hat_pri_par[s_max]/sqrt(n)

  # d = matrix(0,n_r,k_par)
  # d_whl = matrix(0,n_r,k_par)
  # d_par = matrix(0,n_r,k_par)
  d = matrix(0,n_r,k)
  d_whl = matrix(0,n_r,k)
  d_par = matrix(0,n_r,k)
  
  for (i in 1:n_r){
    d[i,] = (1-n_d^(r[i]-0.5))*(BS_hat-beta_hat)
    d_whl[i,] = (1-n_d^(r[i]-0.5))*(BS_hat_pri_whl-beta_hat_pri_whl)
    d_par[i,] = (1-n_d^(r[i]-0.5))*(BS_hat_pri_par-beta_hat_pri_par)
  }
  
  T_b = matrix(0,n_r,B)
  T_b_whl = matrix(0,n_r,B)
  T_b_par = matrix(0,n_r,B)
  
  for (b in 1:B){
    D_boot = mvrnorm(n,theta_hat, diag(sd_hat^2))
    D_boot_whl = mvrnorm(n,theta_hat_pri_whl, diag(sd_hat_pri_whl^2))
    D_boot_par = mvrnorm(n,theta_hat_pri_par, diag(sd_hat_pri_par^2))
    
    theta_boot = apply(D_boot,2,mean)
    Lap_noise_b = Lap_noise_Gaussian(k,e,k_par,theta_hat,sensi_mean,sensi_sd)
    
    theta_boot_whl = theta_boot + Lap_noise_b$whl_mean
    theta_boot_par = theta_boot + Lap_noise_b$par_mean
    
    # beta_boot = theta_boot[1:k_par]
    # beta_boot_whl = theta_boot_whl[1:k_par]
    # beta_boot_par = theta_boot_par[1:k_par]
    beta_boot = theta_boot
    beta_boot_whl = theta_boot_whl
    beta_boot_par = theta_boot_par
    
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

fun_mv_Gassian = function(D,B,level,e,r,BS_true,k_par){
  n = dim(D)[1]
  k = dim(D)[2]
  Dt = t(D)
  DtD = Dt%*%D
  n_total = n
  n_d = n*k
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
    d[i,] = (1-n_d^(r[i]-0.5))*(BS_hat-beta_hat)
    d_whl[i,] = (1-n_d^(r[i]-0.5))*(BS_hat_pri_whl-beta_hat_pri_whl)
    d_par[i,] = (1-n_d^(r[i]-0.5))*(BS_hat_pri_par-beta_hat_pri_par)
  }
  
  # parametric bootstrap
  T_b = matrix(0,n_r,B)
  T_b_whl = matrix(0,n_r,B)
  T_b_par = matrix(0,n_r,B)
  
  for (b in 1:B){
    D_boot = mvrnorm(n,theta_hat, cov_matrix)
    D_boot_whl = mvrnorm(n,theta_hat_pri_whl, cov_matrix_whl)
    D_boot_par = mvrnorm(n,theta_hat_pri_par, cov_matrix_par)
    
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

############################ linear regression ############################
fun_linear_regression = function(D,Y,B,level,e,r,BS_true,k_par){
  ### k_par is the dim of vector of interest
  n = dim(D)[1]
  k = dim(D)[2]
  n_total = n
  Dt = t(D)
  DtD = Dt%*%D
  DtD_inv = solve(DtD)
  DtY = Dt%*%Y
  #ZtY = DtY[1:k_par]
  
  gamma = max(abs(D))/2
  zeta = max(abs(Y))/2
  
  ###################### parameter estimate
  theta_hat = DtD_inv %*% DtY
  beta_hat = theta_hat[1:k_par]
  BS_hat = max(beta_hat)
  s_max = which.max(beta_hat)
  
  Q_hat = 1/n * DtD
  
  Lap_noise = Lap_noise_OLS2(k,gamma,zeta,e,k_par,DtY)
  scl = Lap_noise$scale
  
  Q_hat_pri_whl = Q_hat + 1/n * Lap_noise$whl_DtD
  Q_hat_pri_par = Q_hat + 1/n * Lap_noise$par_DtD
  
  cov_matrix_whl = sigma_sq_hat_pri * Q_hat_pri_whl
  cov_matrix_par = sigma_sq_hat_pri * Q_hat_pri_par
  if (is.positive.definite(cov_matrix_whl) == FALSE){
    cov_matrix_whl = make.positive.definite(cov_matrix_whl, tol=1e-3)
  }
  if (is.positive.definite(cov_matrix_par) == FALSE){
    cov_matrix_par = make.positive.definite(cov_matrix_par, tol=1e-3)
  }
  
  theta_hat_pri_whl = solve(DtD+Lap_noise$whl_DtD)%*%(DtY+Lap_noise$whl_DtY)
  theta_hat_pri_par = solve(DtD+Lap_noise$par_DtD)%*%(DtY+Lap_noise$par_DtY)
  
  beta_hat_pri_whl = theta_hat_pri_whl[1:k_par]
  beta_hat_pri_par = theta_hat_pri_par[1:k_par]
  
  BS_hat_pri_whl = max(beta_hat_pri_whl)
  BS_hat_pri_par = max(beta_hat_pri_par)
  
  width_term = max((zeta - sum(-gamma * abs(theta_hat))) ** 2,(-zeta - sum(gamma * abs(theta_hat))) ** 2)
  Delta_sigma_sq = 1/(n-k) * width_term
  sigma_sq_hat = 1/(n-k) * sum((Y - D%*%theta_hat)^2)
  sigma_sq_hat_pri = sigma_sq_hat + rlaplace(1,0, Delta_sigma_sq/(e/3))
  if (sigma_sq_hat_pri<0){
    sigma_sq_hat_pri = 0.1
  }
  
  lb_naive = beta_hat[s_max] - qt(1-level,df = n_total-1)*sqrt(1.0 * sigma_sq_hat * DtD_inv[s_max,s_max])
  lb_naive_whl = beta_hat_pri_whl[s_max] - qt(1-level,df = n_total-1)*sqrt(1.0 * sigma_sq_hat_pri * solve(DtD+Lap_noise$whl_DtD)[s_max,s_max])
  lb_naive_par = beta_hat_pri_par[s_max] - qt(1-level,df = n_total-1)*sqrt(1.0 * sigma_sq_hat_pri * solve(DtD+Lap_noise$par_DtD)[s_max,s_max])
  
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
    Y_boot = D%*%theta_hat + rnorm(n,0,sqrt(sigma_sq_hat))
    Lap_noise_b = Lap_noise_OLS2(k,gamma,zeta,e,k_par,DtY)
    Z_b_whl = mvrnorm(1,rep(0,k), cov_matrix_whl)
    Z_b_par = mvrnorm(1,rep(0,k), cov_matrix_par)
    Q_hat_b_whl= Q_hat_pri_whl + 1 /n * Lap_noise_b$whl_DtD
    Q_hat_b_par= Q_hat_pri_par + 1 /n * Lap_noise_b$par_DtD
    
    theta_boot = DtD_inv%*% Dt %*% Y_boot
    theta_boot_whl = solve(Q_hat_b_whl)%*%Q_hat_pri_whl%*%theta_hat_pri_whl+solve(Q_hat_b_whl)%*%(1/sqrt(n)*Z_b_whl+1/n * Lap_noise_b$whl_DtY)
    theta_boot_par = solve(Q_hat_b_par)%*%Q_hat_pri_par%*%theta_hat_pri_par+solve(Q_hat_b_par)%*%(1/sqrt(n)*Z_b_par+1/n * Lap_noise_b$par_DtY)
    
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
  #dis_lb = BS_true - lb # distance form the true best subgroup effect
  #result <- list(lb,dis_lb,s_max,lb_naive)
  result <- list(lb,lb_whl,lb_par,lb_naive,lb_naive_whl,lb_naive_par,s_max,scl)
  names(result) <- c("nonpri_bstp","priwhole_bstp","pripartial_bstp","nonpri_naive","priwhole_naive","pripartial_naive","best_subgroup","scale_of_noise")
  return(result)
}

############################ simulation ############################
simulation = function(rep,B,level,e,r,theta_true,k_par,cov_matrix=diag(length(theta_true)),model){
  n_r = length(r)
  count_PB = matrix(0,rep,n_r)
  dis_PB = matrix(0,rep,n_r)
  count_whl = matrix(0,rep,n_r)
  dis_whl = matrix(0,rep,n_r)
  count_par = matrix(0,rep,n_r)
  dis_par = matrix(0,rep,n_r)
  count_naive = 0
  dis_naive = 0
  count_whl_naive = 0
  dis_whl_naive = 0
  count_par_naive = 0
  dis_par_naive = 0
  scl = 0
  
  if (model == 'Gaussian' | model == 'mv_Gaussian'){
    BS_true = max(theta_true)
  }
  else{  
    beta_true = theta_true[1:k_par]
    BS_true = max(beta_true)
  }
  
  for (j in 1:rep){
    if (model == 'Gaussian'){
      D = Gaussian_data(n,theta_true,cov_matrix)
      result = fun_Gassian(D,B,level,e,r,BS_true,k_par)
    }
    
    if (model == 'mv_Gaussian'){
      D = Gaussian_data(n,theta_true,cov_matrix)
      result = fun_mv_Gassian(D,B,level,e,r,BS_true,k_par)
    }
    if (model == 'OLS_data_1'){
      data = OLS_data_1(n,theta_true)
      D = data$D
      Y = data$Y
      result = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
    }
    
    if (model == 'OLS_data_2'){
      data = OLS_data_2(n,k_par,theta_true)
      D = data$D
      Y = data$Y
      result = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
    }
    
    if (model == 'OLS_data_real'){
      data = OLS_data_real()
      D = data$D
      Y = data$Y
      result = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
    }
    
    scl[j] = abs(result$scale_of_noise)
    count_PB[j,] = (result$nonpri_bstp <= BS_true)
    count_whl[j,] = (result$priwhole_bstp <= BS_true)
    count_par[j,] = (result$pripartial_bstp <= BS_true)
    count_naive[j] = (result$nonpri_naive <= BS_true)
    count_whl_naive[j] = (result$priwhole_naive <= BS_true)
    count_par_naive[j] = (result$pripartial_naive <= BS_true)

    dis_PB[j,] = BS_true - result$nonpri_bstp
    dis_whl[j,] = BS_true - result$priwhole_bstp
    dis_par[j,] = BS_true - result$pripartial_bstp
    dis_naive[j] = BS_true - result$nonpri_naive
    dis_whl_naive[j] = BS_true - result$priwhole_naive  
    dis_par_naive[j] = BS_true - result$pripartial_naive
  }
  c_PB = colSums(count_PB)/rep
  avedis_PB = apply(dis_PB,2,mean)
  c_whl = colSums(count_whl)/rep
  c_par= colSums(count_par)/rep
  avedis_whl = apply(dis_whl,2,mean)
  avedis_par = apply(dis_par,2,mean)
  c_naive = mean(count_naive)
  avedis_naive = mean(dis_naive)
  c_whl_naive = mean(count_whl_naive)
  avedis_whl_naive = mean(dis_whl_naive)
  c_par_naive = mean(count_par_naive)
  avedis_par_naive = mean(dis_par_naive)
  scl = mean(scl)
  result <- list(c_PB,avedis_PB,c_whl,avedis_whl,c_par,avedis_par,c_naive,avedis_naive,c_whl_naive,avedis_whl_naive,c_par_naive,avedis_par_naive,scl)
  names(result)<- c("coverage_PB", "length_PB", "coverage_whole_HPPB","length_whole_HPPB","coverage_partial_HPPB","length_partial_HPPB","coverage_naive",
                    "length_naive","coverage_whole_naive","length_whole_naive","coverage_partial_naive","length_partial_naive",
                    "scale_of_noise")
  return(result)
}

############################ application ############################
OLS_application = function(rep,B,level,e,r,group_method,k_par){
  data = OLS_data_real(group_method)
  D = data$D
  Y = data$Y
  Y_logit = data$Y_logit
  
  BS_true = 0
  
  l_PB = matrix(0,rep,n_r)
  l_whl = matrix(0,rep,n_r)
  l_par = matrix(0,rep,n_r)
  l_whl_naive = 0
  l_par_naive = 0
  scl = 0
  tmp = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
  best_subgroup = tmp$best_subgroup
  l_naive = tmp$nonpri_naive
  for (j in 1:rep){
    result = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
    scl[j] = abs(result$scale_of_noise)
    l_PB[j,] = result$nonpri_bstp 
    l_whl[j,] = result$priwhole_bstp 
    l_par[j,] = result$pripartial_bstp
    l_whl_naive[j] = result$priwhole_naive
    l_par_naive[j] = result$pripartial_naive
    
  }
  avedis_l_PB = apply(l_PB,2,mean)
  avedis_l_whl = apply(l_whl,2,mean)
  avedis_l_par = apply(l_par,2,mean)
  avedis_l_whl_naive = mean(l_whl_naive)
  avedis_l_par_naive = mean(l_par_naive)
  scl = mean(scl)
  result <- list(avedis_l_PB,avedis_l_whl,avedis_l_par,l_naive,avedis_l_whl_naive,avedis_l_par_naive,best_subgroup,scl)
  names(result) <- c("nonpri_bstp","priwhole_bstp","pripartial_bstp","nonpri_naive","priwhole_naive","pripartial_naive","best_subgroup","scale_of_noise")
  
  return(result)
}

OLS_grid_search = function(e_grid){
  avg_scl_list = 0
  n_e = length(e_grid)
  for (i in 1:n_e){
    e = e_grid[i]
    avg_scl_list[i] = OLS_application(rep,B,level,e,r,group_method,k_par)$scale_of_noise
  }
  result = data.frame(e_grid, avg_scl_list)
  names(result) = c("list_of_scale", "avg_scale_of_noise")
  return(result)
}