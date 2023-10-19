library(MASS)
library(ATE)
Gaussian_data = function(n,theta_true,cov_matrix){
  D = mvrnorm(n,theta_true, cov_matrix)
  return(D)
}

OLS_data_1 = function(n,theta_true){
  k = length(theta_true)
  D = matrix(0,n,k)
  for (i in 1:n){
    index = sample(1:k,1)
    D[i,index] = 1
  }
  
  Y = D %*% theta_true + rnorm(n)
  
  result = list(D,Y)
  names(result) = c("D","Y")
  return(result)
}

OLS_data_2 = function(n,k_par,theta_true){
  # k: 总维数   k_par: 关心的向量的维数
  k = length(theta_true)
  Sigma = matrix(0,k,k)
  for (j in 1:k){
    for (m in 1:j){
      Sigma[j,m] = 0.5^(abs(j-m))
      Sigma[m,j] = Sigma[j,m]
    }
  }
  
  W = mvrnorm(n, rep(0,k), Sigma)
  X = W[,(k_par+1):k]
  Z = W[,1:k_par]
  for (m in 1:n){
    for (l in 1:k_par){
      Z[m,l] = (Z[m,l]>0.5)
    }
  }
  
  D= cbind(Z,X)
  Y = D %*% theta_true +rnorm(n)
  
  result = list(D,Y)
  names(result) = c("D","Y")
  return(result)
}

OLS_data_real = function(group_method){
  ################## general settings  ################## 
  data(nsw)
  nsw_new = nsw
  nsw_new$race = nsw_new$black * 2 + nsw_new$hisp
  nsw_new <- subset(nsw_new, select = -c(black,hisp))
  nsw_new<-subset(nsw_new,select=c(1,2,3,8,4,5,6,7))
  nsw_new$race = nsw_new$race - 1
  nsw_beifen = nsw
  nsw = nsw_new
  
  X<- nsw[,c(-1,-8)]
  Ti<- nsw$treat
  
  n = dim(nsw)[1]
  Y <- nsw$re78-nsw$re75
  Y = (Y- mean(Y))/sd(Y)
  Y_logit <- (nsw$re78 > nsw$re75)
  ################## try different group_methods   ################## 
  if (group_method == 1){
    n_subgroup <- 2 # 5.3改
    k_n <- 2^n_subgroup
    
    Z <- matrix(0,n,k_n)
    for (i in 1:n){ # subgrou = 3: dim1 = 0, dim2 = 1, which is: race!=0, re75 = 1
      indicator_2base = rep(0,n_subgroup) 
      if (X$race[i] == 0) indicator_2base[1] = indicator_2base[1]+1 # dim1
      if (X$ed[i] > median(X$ed)) indicator_2base[2] = indicator_2base[2]+1  # dim2
      indicator_10base = indicator_2base[2]+indicator_2base[1]*2
      Z[i,indicator_10base+1] = 1
    }
  }
  
  if (group_method == 2){
    n_subgroup <- 2 # 5.3改
    k_n <- 2^n_subgroup
    
    Z <- matrix(0,n,k_n)
    for (i in 1:n){ # subgrou = 3: dim1 = 0, dim2 = 1, which is: race!=0, re75 = 1
      indicator_2base = rep(0,n_subgroup) 
      if (X$race[i] == 0) indicator_2base[1] = indicator_2base[1]+1 # dim1
      if (X$re75[i] > median(X$re75)) indicator_2base[2] = indicator_2base[2]+1  # dim2
      indicator_10base = indicator_2base[2]+indicator_2base[1]*2
      Z[i,indicator_10base+1] = 1
    }
  }
  
  if (group_method == 3){
    n_subgroup <- 2
    k_n <- 2^n_subgroup
    
    Z <- matrix(0,n,k_n)
    for (i in 1:n){ # subgrou = 3: dim1 = 0, dim2 = 1, which is: ed low, re75 high
      indicator_2base = rep(0,n_subgroup) 
      if (X$ed[i] > median(X$ed)) indicator_2base[1] = indicator_2base[1]+1 # dim1
      if (X$re75[i] > median(X$re75)) indicator_2base[2] = indicator_2base[2]+1  # dim2
      indicator_10base = indicator_2base[2]+indicator_2base[1]*2
      Z[i,indicator_10base+1] = 1
    }
  }
  
  if (group_method == 4){
    n_subgroup <- 3
    k_n <- 2^n_subgroup
    
    Z <- matrix(0,n,k_n)
    for (i in 1:n){ # subgrou = 3: dim1 = 0, dim2 = 1, dim3 = 0,which is: ed low, 
      indicator_2base = rep(0,n_subgroup) 
      if (X$ed[i] > median(X$ed)) indicator_2base[1] = indicator_2base[1]+1 # dim1
      if (X$race[i] == 0) indicator_2base[2] = indicator_2base[2]+1 # dim2
      if (X$re75[i] > median(X$re75)) indicator_2base[3] = indicator_2base[3]+1  # dim3
      indicator_10base = indicator_2base[3]+indicator_2base[2]*2+indicator_2base[1]*4
      Z[i,indicator_10base+1] = 1
    }
  }
  ################## general settings  ################## 
  DZ <- matrix(0,n,k_n)
  for (i in 1:n){
    DZ[i,] = Ti[i]*Z[i,]
  }
  
  X <- as.matrix(X)
  X = cbind(Z,X) 
  for (j in 1:dim(X)[2]){
    X[,j] = (X[,j]-mean(X[,j]))/sd(X[,j])
  }
  
  D = cbind(DZ,X) 
  
  result = list(D,Y,Y_logit)
  names(result) = c("D","Y","Y_logit")
  return(result)
}