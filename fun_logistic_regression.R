fun_logistic_regression = function(D,Y_logit,level,r,k_par){
  n = dim(D)[1]
  n_total = n
  
  fit_1 <- glm(formula = Y_logit ~ D-1, family = binomial(link = "logit"))
  fit_summary = summary(fit_1)
  summary_coefficients = fit_summary$coefficients
  beta_hat = summary_coefficients[1:k_par]
  BS_hat = max(beta_hat)
  s_max = which.max(beta_hat)
  scl = Lap_noise$scale
  lb_naive = BS_hat - qt(1-level,df = n-1)*summary_coefficients[s_max,2]
  
  Lap_noise = Lap_noise_Logistic(n, e)

  Y_pri_whl = ((Y + Lap_noise$noise_whl) > 0) 
  Y_pri_par = ((Y + Lap_noise$noise_par) > 0)

  fit <- glm(formula = Y_pri_whl ~ D-1, family = binomial(link = "logit"))
  fit_summary = summary(fit)
  summary_coefficients = fit_summary$coefficients
  beta_hat_pri_whl = summary_coefficients[1:k_par]
  BS_hat_pri_whl = max(beta_hat_pri_whl)
  lb_naive_whl = BS_hat_pri_whl - qt(1-level,df = n-1)*summary_coefficients[s_max,2]
  
  fit <- glm(formula = Y_pri_par ~ D-1, family = binomial(link = "logit"))
  fit_summary = summary(fit)
  summary_coefficients = fit_summary$coefficients
  beta_hat_pri_par = summary_coefficients[1:k_par]
  BS_hat_pri_par = max(beta_hat_pri_par)
  lb_naive_par = BS_hat_pri_par - qt(1-level,df = n-1)*summary_coefficients[s_max,2]
  
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
    Y_boot = predict(fit_1, as.data.frame(D),type = "response")
    Y_boot = (Y_boot > 0.5)
    
    Lap_noise = Lap_noise_Logistic(n, e)
    Y_pri_whl_boot = ((Y_boot + Lap_noise$noise_whl) > 0) 
    Y_pri_par_boot = ((Y_boot + Lap_noise$noise_par) > 0)
    
    fit <- glm(formula = Y_boot ~ D-1, family = binomial(link = "logit"))
    fit_summary = summary(fit)
    summary_coefficients = fit_summary$coefficients
    beta_boot = summary_coefficients[1:k_par]
    
    fit <- glm(formula = Y_pri_whl_boot ~ D-1, family = binomial(link = "logit"))
    fit_summary = summary(fit)
    summary_coefficients = fit_summary$coefficients
    beta_boot_whl = summary_coefficients[1:k_par]
    
    fit <- glm(formula = Y_pri_par_boot ~ D-1, family = binomial(link = "logit"))
    fit_summary = summary(fit)
    summary_coefficients = fit_summary$coefficients
    beta_boot_par = summary_coefficients[1:k_par]

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

OLS_application_logistic = function(rep,B,level,e,r,group_method,k_par){
  data = OLS_data_real(group_method)
  D = data$D
  Y_logit = data$Y_logit
  BS_true = 0
  
  l_PB = matrix(0,rep,n_r)
  l_whl = matrix(0,rep,n_r)
  l_par = matrix(0,rep,n_r)
  l_whl_naive = 0
  l_par_naive = 0
  scl = 0
  tmp = fun_logistic_regression(D,Y_logit,level,r,k_par)
  best_subgroup = tmp$best_subgroup
  l_naive = tmp$nonpri_naive
  for (j in 1:rep){
    result = fun_logistic_regression(D,Y_logit,level,r,k_par)
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
