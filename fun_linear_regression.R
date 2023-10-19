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

  Lap_noise = Lap_noise_OLS2(k,gamma,zeta,e,k_par,DtY)
  scl = Lap_noise$scale
  
  theta_hat_pri_whl = solve(DtD+Lap_noise$whl_DtD)%*%(DtY+Lap_noise$whl_DtY)
  theta_hat_pri_par = solve(DtD+Lap_noise$par_DtD)%*%(DtY+Lap_noise$par_DtY)
  
  beta_hat_pri_whl = theta_hat_pri_whl[1:k_par]
  beta_hat_pri_par = theta_hat_pri_par[1:k_par]
  
  BS_hat_pri_whl = max(beta_hat_pri_whl)
  BS_hat_pri_par = max(beta_hat_pri_par)

  Q_hat = 1/n * DtD
  Q_hat_pri_whl = Q_hat + 1/n * Lap_noise$whl_DtD
  Q_hat_pri_par = Q_hat + 1/n * Lap_noise$par_DtD

  width_term = max((zeta - sum(-gamma * abs(theta_hat))) ** 2,(-zeta - sum(gamma * abs(theta_hat))) ** 2)
  Delta_sigma_sq = 1/(n-k) * width_term
  sigma_sq_hat = 1/(n-k) * sum((Y - D%*%theta_hat)^2)
  sigma_sq_hat_pri = sigma_sq_hat + rlaplace(1,0, Delta_sigma_sq/(e/3))
  if (sigma_sq_hat_pri<0){
    sigma_sq_hat_pri = 0.1
  }

  cov_matrix_whl = sigma_sq_hat_pri * Q_hat_pri_whl
  cov_matrix_par = sigma_sq_hat_pri * Q_hat_pri_par
  if (is.positive.definite(cov_matrix_whl) == FALSE){
    cov_matrix_whl = make.positive.definite(cov_matrix_whl, tol=1e-3)
  }
  if (is.positive.definite(cov_matrix_par) == FALSE){
    cov_matrix_par = make.positive.definite(cov_matrix_par, tol=1e-3)
  }
  
  
  lb_naive = beta_hat[s_max] - qt(1-level,df = n_total-1)*sqrt(1.0 * sigma_sq_hat * DtD_inv[s_max,s_max])
  lb_naive_whl = beta_hat_pri_whl[s_max] - qt(1-level,df = n_total-1)*sqrt(1.0 * sigma_sq_hat_pri * DtD_inv[s_max,s_max])
  lb_naive_par = beta_hat_pri_par[s_max] - qt(1-level,df = n_total-1)*sqrt(1.0 * sigma_sq_hat_pri * DtD_inv[s_max,s_max])

  n_r = length(r)
  d = matrix(0,n_r,k_par)
  d_whl = matrix(0,n_r,k_par)
  d_par = matrix(0,n_r,k_par)

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
  #dis_whl = BS_true - lb_whl
  # dis_par = BS_true - lb_par
  # dis_naive = BS_true - lb_naive
  # dis_naive_whl = BS_true - lb_naive_whl
  # dis_naive_par = BS_true - lb_naive_par
  result <- list(lb,lb_whl,lb_par,lb_naive,lb_naive_whl,lb_naive_par,s_max,scl)
  names(result) <- c("nonpri_bstp","priwhole_bstp","pripartial_bstp","nonpri_naive","priwhole_naive","pripartial_naive","best_subgroup","scale_of_noise")
  
  #result_dis <- list(dis_lb,dis_whl,dis_par)
  #names(result_dis) <- c("nonpri_bstp","priwhole_bstp","pripartial_bstp","nonpri_naive","priwhole_naive","pripartial_naive")
  #result_cov <- list(lb,lb_whl,lb_par,lb_naive,lb_naive_whl,lb_naive_par)
  #names(result_cov) <- c("nonpri_bstp","priwhole_bstp","pripartial_bstp","nonpri_naive","priwhole_naive","pripartial_naive")
  #result = c(result_dis,result_cov,s_max,scl)
  #names(result) = c("distance","coverage","best_subgroup","scale_of_noise")
  return(result)
}
