Lap_noise_OLS2 = function(k,gamma,zeta,e,k_par,DtY){
  sensi_DtY = k * 2*gamma * zeta
  w = rlaplace(k,0,sensi_DtY/(e/3))
  # w = rep(0,k)
  # for (i in 1:k){
  #   w[i] = rlaplace(1,0,sensi_DtY/(e/3))
  # }
  scl = max(abs((w/DtY)))
  sensi_DtD = k*(k + 1) * gamma^2
  V = matrix(0,k,k)
  for (i in 1:k){
    for (j in 1:i){
      V[i,j] = rlaplace(1,0,sensi_DtD/(e/3))
      V[j,i] = V[i,j]
    }
  }
  
  w_p = w
  V_p = V
  if (k_par != k){
    for (i in (k_par +1) : k){
      w_p[i] = 0
      for (j in (k_par+1) :k){
        V_p[i,j] = 0
      }
    }
  }
  
  Lap_noise = list(w,V,w_p,V_p,scl)
  names(Lap_noise) = c("whl_DtY","whl_DtD","par_DtY","par_DtD","scale") # whl for whole
  return(Lap_noise)
} # OLS加噪声

make_pos_def = function(x){
  small_positive = 0.1
  svd <- svd(x)
  svd_u = svd$u
  svd_s = svd$d
  svd_v = svd$v
  
  ev <- eigen(x)
  l = ev$val
  w = ev$vec
  
  if (l[order(l,decreasing = FALSE)[1]]>0){
    xPSD = x
  }
  else {
    l_prime = l
    for (i in 1:length(l)){
      if (l[i] < 0 ){
        l_prime[i] = small_positive
      }
    }
    xPSD = w%*%diag(l_prime)%*%t(x)
  }
  return(xPSD)
}

Lap_noise_Logistic = function(n, e){
  noise_whl = rlaplace(n, 1/e)
  noise_par = noise_whl
  for (i in (k_par+1):length(noise_par)){
    noise_par[i] = 0
  }
  scl = max(noise_whl)
  result = list(noise_whl,noise_par,scl)
  names(result) = c("noise_whl","noise_par","scale")
  return(result)
}

library(plyr)
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}

Lap_noise_mv_Gaussian = function(k,e,k_par,theta_hat,sensi_mean,gamma){
  w = rlaplace(k,0,sensi_mean/(e/2))
  scl = max(abs(w/theta_hat))
  sensi_DtD = k*(k + 1) * gamma^2
  V = matrix(0,k,k)
  for (i in 1:k){
    for (j in 1:i){
      V[i,j] = rlaplace(1,0,sensi_DtD/(e/3))
      V[j,i] = V[i,j]
    }
  }
  
  w_p = w
  V_p = V

  if (k_par != k){
    for (i in (k_par +1) : k){
      w_p[i] = 0
      for (j in (k_par+1) :k){
        V_p[i,j] = 0
      }
    }
  }
  
  Lap_noise = list(w,V,w_p,V_p,scl)
  names(Lap_noise) = c("whl_mean","whl_cov","par_mean","par_cov","scale") 
  return(Lap_noise)
}

Lap_noise_Gaussian = function(k,e,k_par,theta_hat,sensi_mean,sensi_sd){
  w = rlaplace(k,0,sensi_mean/(e/2))
  scl = max(abs(w/theta_hat))
  V = rlaplace(k,0,sensi_sd/(e/2))
  
  w_p = w
  V_p = V
  
  if (k_par != k){
    for (i in (k_par +1) : k){
      w_p[i] = 0
      V_p[i] = 0
    }
  }
  
  Lap_noise = list(w,V,w_p,V_p,scl)
  names(Lap_noise) = c("whl_mean","whl_cov","par_mean","par_cov","scale") 
  return(Lap_noise)
} # 未完待续

estimate_Gaussian = function(e,theta_hat,pri){
  Lap_noise = Lap_noise_Gaussian(k,e,k_par,theta_hat,sensi_mean,sensi_sd)
  scl = Lap_noise$scale
  lb_naive = beta_hat[s_max] - qt(1-level,df = n_total-1)*sd_hat[s_max]/sqrt(n)
  if (pri == 'whl'){
    theta_hat_pri = theta_hat + Lap_noise$whl_mean
    sd_hat_pri = sd_hat + Lap_noise$whl_cov
    lb_naive = beta_hat_pri_whl[s_max] - qt(1-level,df = n_total-1)*sd_hat_pri_whl[s_max]/sqrt(n)
  }
  if (pri == 'par'){
    theta_hat_pri = theta_hat + Lap_noise$par_mean
    sd_hat_pri = sd_hat + Lap_noise$par_cov
  }
  beta_hat_pri = theta_hat_pri
  BS_hat_pri = max(beta_hat_pri)
}