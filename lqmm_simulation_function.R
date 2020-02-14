###############################################################################################################
#                                                                                                             #
#                          MODEL WITH ONE COVARIATE AND ONE RANDOM INTERCEPT                                  #
#                                                                                                             #
###############################################################################################################

#--------------------------FUNCTION TO CALCULATE RELATIVE BIAS FOR FIXED EFFECTS AND VARIANCE OF RANDOM EFFECTS
#--------------------------ERROR ~ N(0,1)
lqmm_enorm <- function(B, N, ni, beta_0, beta_1, sdi, sd, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=4))
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_i = rep(rnorm(N, 0, sdi), each = ni)
    Replicates <- rep(1:ni, N)
    x_ij = rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + (beta_1*x_ij) + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, x_ij, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ x_ij, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:3){
        estimates[i, k] <- fit$theta[k]
      }
      estimates[i,4] <- VarCorr(fit)[[1]]
  }
  colnames(estimates) <- c("beta0", "beta1", "b0", "var_b0")

  #calculate relative bias
  bias_beta_0 <- NULL
  bias_beta_1 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)/abs(estimates$beta0[j])
    bias_beta_1[j] = (estimates$beta1[j] - beta_1)/abs(estimates$beta1[j])
    bias_b0[j] = (estimates$b0[j] - sd)/abs(estimates$b0[j])
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_beta1 = mean(bias_beta_1)
  mean_bias_b0 = mean(bias_b0)
  result <- list(estimates = estimates, bias_beta0 = mean_bias_beta0, bias_beta0 = mean_bias_beta1,
                 bias_b0 = mean_bias_b0)
  return(result)
}


#--------------------------FUNCTION TO CALCULATE RELATIVE BIAS FOR RANDOM EFFECTS
#--------------------------ERROR ~ N(0,1)

lqmm_ranef_enorm <- function(B, N, ni, beta_0, beta_1, sdi, sd, tau) {
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rnorm(N, 0, sdi)
    u_i = rep(u_ij, each = ni)
    Replicates <- rep(1:ni, N)
    x_ij = rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + (beta_1*x_ij) + u_i +error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, x_ij, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ x_ij, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
  mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(results)
}


###############################################################################################################
#                                                                                                             #
#                          MODEL WITHOUT COVARIATE AND ONE RANDOM INTERCEPT                                   #
#                                                                                                             #
###############################################################################################################
 
library(lqmm)
#--------------------------FUNCTION TO CALCULATE RELATIVE BIAS FOR FIXED EFFECTS AND VARIANCE OF RANDOM EFFECTS
#--------------------------ERROR_E ~ N(0,1), ERROR_U ~ N(0,1)

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
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")

  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)

  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}
lqmm_int_enorm(2, 5, 3, 0, 1,1, 0.025)

#--------------------------ERROR_E ~ LNORM(0,1), ERROR_U ~ N(0,1)

lqmm_int_elnorm <- function(B, N, ni, beta_0, sdi, sd, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rnorm(N, 0, sdi)
    u_i = rep(u_ij, each = ni)
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
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)
  
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}


#--------------------------ERROR_E ~ NORM(0,1), ERROR_U ~ LNORM(0,1)

lqmm_int_ulnorm <- function(B, N, ni, beta_0, sdi, sd, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rlnorm(N, 0, sdi)
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
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)
  
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}


#--------------------------ERROR_E ~ LNORM(0,1), ERROR_U ~ LNORM(0,1)

lqmm_int_eulnorm <- function(B, N, ni, beta_0, sdi, sd, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rlnorm(N, 0, sdi)
    u_i = rep(u_ij, each = ni)
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
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)
  
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}


#--------------------------ERROR_E ~ UNIFORM(1,2), ERROR_U ~ N(0,1)

lqmm_int_eunif <- function(B, N, ni, beta_0, sdi, min, max, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rnorm(N, 0, sdi)
    u_i = rep(u_ij, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = runif(N*ni, min, max)
    Measure = beta_0 + u_i + error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)
  
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}


#--------------------------ERROR_E ~ NORM(0,1), ERROR_U ~ UNIFORM(1,2)

lqmm_int_uunif <- function(B, N, ni, beta_0, sd, min, max, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = runif(N, min = min, max = max)
    u_i = rep(u_ij, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i + error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  ?runif
  #calculate relative bias
  var_ui = 1/12*((max - min)^2)
  sdi = sqrt(var_ui)
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)
  
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}


#--------------------------ERROR_E ~ CHISQ(4), ERROR_U ~ NORM(0,1)

lqmm_int_echisq <- function(B, N, ni, beta_0, sdi, df, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rnorm(N, 0, sdi)
    u_i = rep(u_ij, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rchisq(N*ni, df)
    Measure = beta_0 + u_i + error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate relative bias
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)
  
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}


#--------------------------ERROR_E ~ NORM(0,1), ERROR_U ~ CHISQ(4)

lqmm_int_uchisq <- function(B, N, ni, beta_0, sd, df, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rchisq(N, df)
    u_i = rep(u_ij, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rnorm(N*ni, 0, sd)
    Measure = beta_0 + u_i + error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate relative bias
  var_ui = 2*df
  sdi = sqrt(var_ui)
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)
  
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}


#--------------------------ERROR_E ~ CHISQ(4), ERROR_U ~ CHISQ(4)

lqmm_int_euchisq <- function(B, N, ni, beta_0, df, dfi, tau){
  estimates <- data.frame(matrix(0, nrow = B, ncol=3))
  result <- list()
  for (i in 1:B) {
    Subject <- rep(1:N, each = ni)
    u_ij = rchisq(N, dfi)
    u_i = rep(u_ij, each = ni)
    Replicates <- rep(1:ni, N)
    error_1 = rchisq(N*ni, df)
    Measure = beta_0 + u_i + error_1
    data_model_1 = data.frame(Subject, u_i, Replicates, error_1, Measure)
    fit <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau,
                nK = 7, type = "normal", data = data_model_1, control = lqmmControl(LP_max_iter = 1000,
                                                                                    LP_tol_ll = 1e-7))
    for(k in 1:2){
      estimates[i, k] <- fit$theta[k]
    }
    estimates[i,3] <- VarCorr(fit)[[1]]
    
    #random effect calclulation
    ranef_lqmm = data.frame(u_ij, ranef(fit))
    colnames(ranef_lqmm) <- c("u_ij", "u_ij_hat")  
    #store the result in a new list defined above
    result <- append(result, ranef_lqmm)
  }
  colnames(estimates) <- c("beta0", "b0", "var_b0")
  
  #calculate relative bias
  var_ui = 2*dfi
  sdi = sqrt(var_ui)
  bias_beta_0 <- NULL
  bias_b0 <- NULL
  for (j in 1:B) {
    bias_beta_0[j] = (estimates$beta0[j] - beta_0)
    bias_b0[j] = (estimates$b0[j] - sdi)
  }
  mean_bias_beta0 = mean(bias_beta_0)
  mean_bias_b0 = mean(bias_b0)
  bias <- list(mean_bias_beta0 = mean_bias_beta0, mean_bias_b0 = mean_bias_b0)
  
  result <- data.frame(result)
  result <- result[, order(names(result))]
  u_i_sum <- data.frame(matrix(0, ncol = B, nrow =  N))
  mean_bias_i <- NULL
  for (k in 1:N) {
    for (j in 1:B) {
      u_i_sum[,j] = result[,j] - result[,j+B]
    }
    mean_bias_i[k] =sum(u_i_sum[k,])/B
  }
  avg_bias <- mean(mean_bias_i)
  results <- list(bias_i = mean_bias_i, mean_bias = avg_bias)
  return(list(results_ranef = results, bias_estimates = bias))
}