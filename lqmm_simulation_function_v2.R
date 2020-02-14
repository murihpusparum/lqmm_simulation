#--------------------------FUNCTION TO CALCULATE RELATIVE BIAS FOR FIXED EFFECTS AND VARIANCE OF RANDOM EFFECTS
#1 --------------------------ERROR_E ~ N(0,1), ERROR_U ~ N(0,1)

lqmm_int_enorm <- function(B, N, ni, beta_0, sdi, sd, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rnorm(N, 0, sdi)
    u_i = rep(u_ij, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")

  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))

  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)

  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)

  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
    
  return(list(MSE = MSE, ind_ci = individual_ci))
}


#2 --------------------------ERROR_E ~ NORM(0,1), ERROR_U ~ LNORM(0,1)

lqmm_int_ulnorm <- function(B, N, ni, beta_0, sdi, sd, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  u_ij = rlnorm(N, 0, sdi)
  u_i = rep(u_ij, each = ni)
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))
  
  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)
  
  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
  
  return(list(MSE = MSE, ind_ci = individual_ci))
}


#3 --------------------------ERROR_E ~ LNORM(0,1), ERROR_U ~ LNORM(0,1)

lqmm_int_eulnorm <- function(B, N, ni, beta_0, sdi, sd, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  set.seed(101)
  u_ij = rlnorm(N, 0, sdi)
  u_i = rep(u_ij, each = ni)
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rlnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))
  
  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)
  
  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
  
  return(list(MSE = MSE, ind_ci = individual_ci))
}


#4 --------------------------ERROR_E ~ NORM(0,1), ERROR_U ~ UNIFORM(1,2)

lqmm_int_uunif <- function(B, N, ni, beta_0, sd, min, max, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  u_ij = runif(N, min = min, max = max)
  u_i = rep(u_ij, each = ni)
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))
  
  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)
  
  #calculate relative bias
  var_ui = 1/12*((max - min)^2)
  sdi = sqrt(var_ui)
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
  
  return(list(MSE = MSE, ind_ci = individual_ci))
}


#5 --------------------------ERROR_E ~ NORM(0,1), ERROR_U ~ CHISQ(4)

lqmm_int_uchisq <- function(B, N, ni, beta_0, sd, dfi, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  u_ij = rchisq(N, dfi)
  u_i = rep(u_ij, each = ni)
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))
  
  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)
  
  #calculate relative bias
  var_ui = 2*dfi
  sdi = sqrt(var_ui)
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
  
  return(list(MSE = MSE, ind_ci = individual_ci))
}


#6 --------------------------ERROR_E ~ CHISQ(4), ERROR_U ~ CHISQ(4)

lqmm_int_euchisq <- function(B, N, ni, beta_0, df, dfi, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  u_ij = rchisq(N, dfi)
  u_i = rep(u_ij, each = ni)
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rchisq(N*ni, df)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))
  
  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)
  
  #calculate relative bias
  var_ui = 2*dfi
  sdi = sqrt(var_ui)
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
  
  return(list(MSE = MSE, ind_ci = individual_ci))
}

#7 --------------------------ERROR_E ~ NORM(0,1), ERROR_U ~ t(3)

lqmm_int_utdist <- function(B, N, ni, beta_0, sd, dfi, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  u_ij = rt(N, dfi)
  u_i = rep(u_ij, each = ni)
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))
  
  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)
  
  #calculate relative bias
  var_ui = dfi / (dfi - 2)
  sdi = sqrt(var_ui)
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
  
  return(list(MSE = MSE, ind_ci = individual_ci))
}


#8 --------------------------ERROR_E ~ NORM(0,1), ERROR_U ~ t(5)

lqmm_int_utdist <- function(B, N, ni, beta_0, sd, dfi, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  u_ij = rt(N, dfi)
  u_i = rep(u_ij, each = ni)
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))
  
  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)
  
  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
  
  return(list(MSE = MSE, ind_ci = individual_ci))
}


#9 --------------------------ERROR_E ~ CHISQ(4), ERROR_U ~ t(3)

lqmm_int_echisq_ut <- function(B, N, ni, beta_0, sd, df, dfi, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  u_ij = rt(N, dfi)
  u_i = rep(u_ij, each = ni)
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rchisq(N*ni, df)
    Measure = beta_0 + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect caclulation
    ranef_lqmm = data.frame(ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(u_ij, result)
  result <- result[, order(names(result))]
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate ci estimate for each subject
  beta_new <- data.frame(matrix(0, ncol= B, nrow = N))
  hat_ci <- data.frame(matrix(0, ncol= B, nrow = N))
  
  for (j in 1:B) {
    beta_new[,j] = rep(estimates$beta0[j], N)
    hat_ci[,j] = beta_new[,j] + result[, j+1]
  }
  colnames(beta_new) <- paste("beta0_hat", 1:B, sep="")
  colnames(hat_ci) <- paste("ci", 1:B, sep="")
  
  true_ci <- beta_0 + result[, 1]
  individual_ci <- data.frame(hat_ci, true_ci)
  
  #calculate relative bias
  var_ui = dfi / (dfi - 2)
  sdi = sqrt(var_ui)
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (beta_0 - estimates$beta0[j])^2
    bias_b0[j] = (sdi - estimates$b0[j])^2
  }
  MSE_beta0 = mean(bias_beta_0)
  MSE_b0 = mean(bias_b0)
  
  #calculate mean of bias for all ui estimate
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+1] #ui - ui_hat
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sumsq[,j] = (result[,j] - result[,j+1])^2 #(ui - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  MSE <- list(beta0 = MSE_beta0, b0 = MSE_b0, ui = MSE_ui, mean_ui = mean_MSE_ui)
  
  return(list(MSE = MSE, ind_ci = individual_ci))
}


#--------------------------CALCULATE RELATIVE WIDTH

relative_width <- function(data1, data2, B, N){
  width <- data.frame(matrix(0, ncol = B, nrow = N))
  for (j in 1:B) {
    width[, j] <- data2[, j] - data1[, j] #upper - lower
  }
  width_true <- data2[, B+1] - data1[, B+1]
  
  average_width<-NULL
  for (k in 1:N) {
    average_width[k] <-sum(width[k,])/B 
  }
  
  eval1<-data.frame(average_width/width_true)
  colnames(eval1) <- "relative_width"
  return(eval1)
}



rt(5, 2)
