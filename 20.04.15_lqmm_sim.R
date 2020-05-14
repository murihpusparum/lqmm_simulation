lqmm_cal <- function(B, N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2, FUN){
  estimates1 <- data.frame(matrix(0, nrow = B, ncol=3))
  estimates2 <- data.frame(matrix(0, nrow = B, ncol=3))
  
  warn1a <- NA
  warn_out1 <- NA
  warn2a <- NA
  warn_out2 <- NA
  
  ranef1 <- list()
  ranef2 <- list()
  uij <- list()
  
  catch_lqmm <- function(theta, cov_name){
    out <- tryCatch(covHandling(theta = theta, n = 1, cov_name = cov_name, quad_type = "normal"), 
                    error=function(e){
                      message("catch error")
                      print(e)
                    },
                    warning=function(w){
                      message("catch warning")
                      print(w)
                    })
    out <- inherits(out, what = "warning")
    return(out)
  }
  
  for (i in 1:B) {
    data <- FUN(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2)
    
    fit1 <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau1,
                 nK = 11, type = "normal", data = data$d, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    
    fit2 <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau2,
                 nK = 11, type = "normal", data = data$d, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    
    warn_out1[i] <- catch_lqmm(theta = fit1$theta_z, cov_name = fit1$cov_name)
    warn_out2[i] <- catch_lqmm(theta = fit2$theta_z, cov_name = fit2$cov_name)
    
    warn1a[i] <- catch_lqmm(theta = fit1$theta_z, cov_name = fit1$cov_name)
    warn2a[i] <- catch_lqmm(theta = fit2$theta_z, cov_name = fit2$cov_name)
    
    warn1 <- catch_lqmm(theta = fit1$theta_z, cov_name = fit1$cov_name)
    warn2 <- catch_lqmm(theta = fit2$theta_z, cov_name = fit2$cov_name)
    
    j <- 1
    while (warn1[j] == TRUE | warn2[j] == TRUE) {
#      rm(.Random.seed)
      data <- FUN(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2)
      
      fit1 <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau1,
                   nK = 11, type = "normal", data = data$d, 
                   control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                         beta = 0.5, gamma = 1))
      fit2 <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau2,
                   nK = 11, type = "normal", data = data$d, 
                   control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                         beta = 0.5, gamma = 1))
      
      j <- j+1
      warn1[j] <- catch_lqmm(theta = fit1$theta_z, cov_name = fit1$cov_name)
      warn2[j] <- catch_lqmm(theta = fit2$theta_z, cov_name = fit2$cov_name)
      
    }
    warn1a <- append(warn1a, warn1)
    warn2a <- append(warn2a, warn2)
    
    
    for(k in 1:2){
      estimates1[i, k] <- fit1$theta[k]
      estimates2[i, k] <- fit2$theta[k]
    }
    estimates1[i,3] <- VarCorr(fit1)[[1]]
    estimates2[i,3] <- VarCorr(fit2)[[1]]
    
    #store real random effect
    u_ij <- data$u
    u_ij = data.frame(u_ij)
    uij <- append(uij, u_ij)
    uij <- data.frame(uij)
    
    #random effect caclulation
    ranef_lqmm1 = data.frame(ranef(fit1))
    colnames(ranef_lqmm1) <- c("u_ij_hat")  
    #store the result in a new list defined above
    ranef1 <- append(ranef1, ranef_lqmm1)
    
    #random effect caclulation
    ranef_lqmm2 = data.frame(ranef(fit2))
    colnames(ranef_lqmm2) <- c("u_ij_hat")  
    #store the result in a new list defined above
    ranef2 <- append(ranef2, ranef_lqmm2)
  }
  
  colnames(estimates1) <- c("beta01", "b01", "var_b01")
  colnames(estimates2) <- c("beta02", "b02", "var_b02")
  
  ranef1 <- data.frame(ranef1)
  ranef2 <- data.frame(ranef2)
  
  #modify data with outliers
  estimates <- cbind.data.frame(estimates1, estimates2)
  estimates$z1 <- (estimates$beta01 - mean(estimates$beta01))/sd(estimates$beta01)
  estimates$z2 <- (estimates$beta02 - mean(estimates$beta02))/sd(estimates$beta02)
  est <- estimates[estimates$z1 >= -2 & estimates$z1 <= 2 & estimates$z2 >= -2 & estimates$z2 <= 2, ]
  len_est <- length(est$beta01)
  
  rem <- NULL
  for (i in 1:B) {
    if (estimates$z1[i] >= -2 & estimates$z1[i] <= 2 & estimates$z2[i] >= -2 & estimates$z2[i] <= 2){
      rem[i] = NA
    }
    else {rem[i] = i}
  }
  rem <- as.numeric(na.omit(rem))
  if(length(rem != 0)){
    ranef1 = ranef1[, -rem]
    ranef2 = ranef2[, -rem]
    uij    = uij[, -rem]
  } else {
    ranef1 = ranef1
    ranef2 = ranef2
    uij    = uij
  }
  
  len_ran1 <- ncol(ranef1)
  len_ran2 <- ncol(ranef2)
  
  #calculate ci estimate for each subject
  B1 <- length(est$beta01)
  beta_new1 <- data.frame(matrix(0, ncol= B1, nrow = N))
  lower <- data.frame(matrix(0, ncol= B1, nrow = N))
  beta_new2 <- data.frame(matrix(0, ncol= B1, nrow = N))
  upper <- data.frame(matrix(0, ncol= B1, nrow = N))
  width_ci <- data.frame(matrix(0, ncol= B1, nrow = N))
  
  
  for (j in 1:B1) {
    beta_new1[,j] = rep(est$beta01[j], N)
    lower[,j] = beta_new1[,j] + ranef1[, j]
    beta_new2[,j] = rep(est$beta02[j], N)
    upper[,j] = beta_new2[,j] + ranef2[, j]
    width_ci[, j] <- upper[, j] - lower[, j]
  }
  
  colnames(beta_new1) <- paste("beta0_hat", 1:B1, sep="")
  colnames(lower) <- paste("low", 1:B1, sep="")
  
  colnames(beta_new2) <- paste("beta0_hat", 1:B1, sep="")
  colnames(upper) <- paste("up", 1:B1, sep="")
  
  colnames(width_ci) <- paste("width_ci", 1:B1, sep="")
  
  all_width_ci <- data.frame(width_ci)
  
  WR_j <- NULL
  for (i in 1:N) {
    WR_j[i] = (sum(all_width_ci[i, ]) / B1) / data$d$true_width[1]
  }
  WR <- mean(WR_j)
  
  
  #calculate relative bias and MSE
  bias_beta_01 <- NULL
  bias_b01 <- NULL
  bias_beta_02 <- NULL
  bias_b02 <- NULL
  
  rel_bias_beta01 <- NULL
  rel_bias_beta02 <- NULL
  for (j in 1:B1) {
    bias_beta_01[j] = (data$d$beta_01[1] - est$beta01[j])^2
    bias_b01[j] = (data$s - est$b01[j])^2 #sdi - hat(sdi)
    bias_beta_02[j] = (data$d$beta_02[1] - est$beta02[j])^2
    bias_b02[j] = (data$s - est$b02[j])^2
    
    rel_bias_beta01[j] = (est$beta01[j] - data$d$beta_01[1]) / abs(data$d$beta_01[1])
    rel_bias_beta02[j] = (est$beta02[j] - data$d$beta_02[1]) / abs(data$d$beta_02[1])
  }
  lower_MSE_beta0 = mean(bias_beta_01)
  lower_MSE_b0 = mean(bias_b01)
  upper_MSE_beta0 = mean(bias_beta_02)
  upper_MSE_b0 = mean(bias_b02)
  
  lower_relbias_beta0 = mean(rel_bias_beta01)
  upper_relbias_beta0 = mean(rel_bias_beta02)
  
  
  #calculate MSE/mean of MSE for all ui estimate
  u_i_sumsq <- data.frame(matrix(0, ncol = B1, nrow =  N))
  MSE_ui <- NULL
  for (k in 1:N) {
    for (j in 1:B1) {
      u_i_sumsq[,j] = (uij[, j] - ranef1[, j])^2 #(uij - ui_hat)^2
    }
    MSE_ui[k] =sum(u_i_sumsq[k,])/B1
  }
  mean_MSE_ui <- mean(MSE_ui)
  
  lower_MSE <- list(beta0 = lower_MSE_beta0, b0 = lower_MSE_b0)
  upper_MSE <- list(beta0 = upper_MSE_beta0, b0 = upper_MSE_b0)
  return(list(MSE_tau_0.025 = lower_MSE, MSE_tau_0.975 = upper_MSE, rel_width = WR, 
              relbias_beta0_tau0.025 = lower_relbias_beta0, relbias_beta0_tau0.975 = upper_relbias_beta0, 
              w1=warn_out1, w1a=warn1a, w2=warn_out2, w2a=warn2a, est = est,
              B1 = B1))
  
}


#1 for both e_ij and u_i follow t
sim_t <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  Subject <- rep(1:N, each = ni)
  u_ijs = rt(N, df=dfi)
  u_ij = u_ijs - 0
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = rt(N*ni, df = df) - 0
  beta_01 <- alpha_0 + qt(tau1, df=df)
  beta_02 <- alpha_0 + qt(tau2, df=df)
  Measure = alpha_0 + e_ij + u_i                   
  
  true_width <- qt(tau2, df=df) - qt(tau1, df=df)
  var_ui = dfi/(dfi - 2)
  sdi = sqrt(var_ui)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#1a for both e_ij and u_i follow t - scaled
sim_t_sc <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  var_ui = dfi/(dfi - 2)
  sdi = sqrt(var_ui)
  
  var_e = df/(df - 2)
  sd = sqrt(var_e)
  
  Subject <- rep(1:N, each = ni)
  u_ijs = rt(N, df=dfi)
  u_ij = (u_ijs - 0)/sdi
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = (rt(N*ni, df = df) - 0)/sd
  beta_01 <- alpha_0 + qt(tau1, df=df)
  beta_02 <- alpha_0 + qt(tau2, df=df)
  Measure = alpha_0 + e_ij + u_i                   
  
  true_width <- qt(tau2, df=df) - qt(tau1, df=df)
 
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#2 for both e_ij and u_i follow N
sim_N <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  Subject <- rep(1:N, each = ni)
  Replicates <- rep(1:ni, N)
  
  u_ijs = rnorm(N, mean=meanij, sd=sdi)                 
  u_ij = u_ijs - meanij                                 
  u_i = rep(u_ij, each = ni)                            
  e_ij = rnorm(N*ni, mean = meanx, sd = sd) - meanx
  Measure = alpha_0 + e_ij + u_i                   
  
  beta_01 <- alpha_0 + qnorm(tau1, mean=meanx, sd=sd)    
  beta_02 <- alpha_0 + qnorm(tau2, mean=meanx, sd=sd)   
  
  true_width <- qnorm(tau2, mean=meanx, sd=sd) - qnorm(tau1, mean=meanx, sd=sd)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#3 for e_ij follow N and u_i follow t
sim_N_t <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  var_ui = dfi/(dfi - 2)
  sdi = sqrt(var_ui)
  
  Subject <- rep(1:N, each = ni)
  u_ijs = rt(N, df=dfi)
  u_ij = (u_ijs - 0)/sdi
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = rnorm(N*ni, mean = meanx, sd = sd) - meanx
  beta_01 <- alpha_0 + qnorm(tau1, mean=meanx, sd=sd)
  beta_02 <- alpha_0 + qnorm(tau2, mean=meanx, sd=sd)
  Measure = alpha_0 + e_ij + u_i                   
  
  true_width <- qnorm(tau2, mean=meanx, sd=sd) - qnorm(tau1, mean=meanx, sd=sd)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#4 for e_ij follow t and u_i follow N
sim_t_N <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  var_e = df/(df - 2)
  sd = sqrt(var_e)
  
  Subject <- rep(1:N, each = ni)
  u_ijs = rnorm(N, mean=meanij, sd=sdi)
  u_ij = u_ijs - meanij
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = (rt(N*ni, df = df) - 0)/sd
  beta_01 <- alpha_0 + qt(tau1, df=df)
  beta_02 <- alpha_0 + qt(tau2, df=df)
  Measure = alpha_0 + e_ij + u_i   
  
  true_width <- qt(tau2, df=df) - qt(tau1, df=df)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#5 for both e_ij and u_i follow Chisq
sim_chi <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  Subject <- rep(1:N, each = ni)
  Replicates <- rep(1:ni, N)
  
  u_ijs = rchisq(N, df=dfi)
  u_ij = u_ijs - dfi
  u_i = rep(u_ij, each = ni)
  e_ij = rchisq(N*ni, df = df) - df
  Measure = alpha_0 + e_ij + u_i   
  
  beta_01 <- alpha_0 + qchisq(tau1, df=df)
  beta_02 <- alpha_0 + qchisq(tau2, df=df)
  true_width <- qchisq(tau2, df=df) - qchisq(tau1, df=df)
  
  var_ui = 2*dfi
  sdi = sqrt(var_ui)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#5a for both e_ij and u_i follow Chisq
sim_chi_sc <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  Subject <- rep(1:N, each = ni)
  Replicates <- rep(1:ni, N)
  
  u_ijs = rchisq(N, df=dfi)
  u_ij = (u_ijs - dfi)/(sqrt(2*dfi))
  u_i = rep(u_ij, each = ni)
  e_ij = (rchisq(N*ni, df = df) - df)/(sqrt(2*df))
  Measure = alpha_0 + e_ij + u_i   
  
  beta_01 <- alpha_0 + qchisq(tau1, df=df)
  beta_02 <- alpha_0 + qchisq(tau2, df=df)
  true_width <- qchisq(tau2, df=df) - qchisq(tau1, df=df)
  
  var_ui = 2*dfi
  sdi = sqrt(var_ui)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#6 for e_ij follow N and u_i follow Chisq
sim_N_chi <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  Subject <- rep(1:N, each = ni)
  u_ijs = rchisq(N, df=dfi)
  u_ij = (u_ijs - dfi)/(sqrt(2*dfi))
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = rnorm(N*ni, mean = meanx, sd = sd) - meanx
  beta_01 <- alpha_0 + qnorm(tau1, mean=meanx, sd=sd)
  beta_02 <- alpha_0 + qnorm(tau2, mean=meanx, sd=sd)
  Measure = alpha_0 + e_ij + u_i   
  
  true_width <- qnorm(tau2, mean=meanx, sd=sd) - qnorm(tau1, mean=meanx, sd=sd)
  
  var_ui = 2*dfi
  sdi = sqrt(var_ui)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#7 for e_ij follow Chisq and u_i follow N
sim_chi_N <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  Subject <- rep(1:N, each = ni)
  u_ijs = rnorm(N, mean=meanij, sd=sdi)
  u_ij = u_ijs - meanij
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = (rchisq(N*ni, df = 4) - df)/(sqrt(2*df))
  beta_01 <- alpha_0 + qchisq(tau1, df=df)
  beta_02 <- alpha_0 + qchisq(tau2, df=df)
  Measure = alpha_0 + e_ij + u_i   
  
  true_width <- qchisq(tau2, df=df) - qchisq(tau1, df=df)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#8 for both e_ij and u_i follow Lognormal(0, 1)
sim_lnorm <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  Subject <- rep(1:N, each = ni)
  u_ijs = rlnorm(N, mean=meanij, sd=sdi)
  u_ij = u_ijs - exp(meanij + ((sdi^2)/2))
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = rlnorm(N*ni, mean=meanx, sd=sd) - exp(meanx + ((sd^2)/2))
  beta_01 <- alpha_0 + qlnorm(tau1, mean=meanx, sd=sd)
  beta_02 <- alpha_0 + qlnorm(tau2, mean=meanx, sd=sd)
  Measure = alpha_0 + e_ij + u_i   
  
  true_width <- qlnorm(tau2, mean=meanx, sd=sd) - qlnorm(tau1, mean=meanx, sd=sd)
  
  var_ui = exp((2*meanij)+(sdi^2))*(exp(sdi^2)-1)
  sdi = sqrt(var_ui)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi))
}

#8a for both e_ij and u_i follow Lognormal(0, 1) - scaled
sim_lnorm_sc <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  
  var_ui = exp((2*meanij)+(sdi^2))*(exp(sdi^2)-1)
  sdi_sc = sqrt(var_ui)
  
  var_e = exp((2*meanx)+(sd^2))*(exp(sd^2)-1)
  sd_sc = sqrt(var_e)
  
  Subject <- rep(1:N, each = ni)
  u_ijs = rlnorm(N, mean=meanij, sd=sdi)
  u_ij = (u_ijs - exp(meanij + ((sdi^2)/2)))/sdi_sc
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = (rlnorm(N*ni, mean=meanx, sd=sd) - exp(meanx + ((sd^2)/2)))/sd_sc
  beta_01 <- alpha_0 + qlnorm(tau1, mean=meanx, sd=sd)
  beta_02 <- alpha_0 + qlnorm(tau2, mean=meanx, sd=sd)
  Measure = alpha_0 + e_ij + u_i   
  
  true_width <- qlnorm(tau2, mean=meanx, sd=sd) - qlnorm(tau1, mean=meanx, sd=sd)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi_sc))
}

#9 for e_ij follow N and u_i follow Lognormal
sim_N_lnorm <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  var_ui = exp((2*meanij)+(sdi^2))*(exp(sdi^2)-1)
  sdi_sc = sqrt(var_ui)
  
  Subject <- rep(1:N, each = ni)
  u_ijs = rlnorm(N, mean=meanij, sd=sdi)
  u_ij = (u_ijs - exp(meanij + ((sdi^2)/2)))/sdi_sc
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = rnorm(N*ni, mean = meanx, sd = sd) - meanx
  beta_01 <- alpha_0 + qnorm(tau1, mean=meanx, sd=sd)
  beta_02 <- alpha_0 + qnorm(tau2, mean=meanx, sd=sd)
  Measure = alpha_0 + e_ij + u_i   
  
  true_width <- qnorm(tau2, mean=meanx, sd=sd) - qnorm(tau1, mean=meanx, sd=sd)
  
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi_sc))
}

#10 for e_ij follow Lognormal and u_i follow t
sim_lnorm_t <- function(N, ni, alpha_0, meanx, meanij, sd, sdi, df, dfi, tau1, tau2){
  var_ui = dfi/(dfi - 2)
  sdi_sc = sqrt(var_ui)
  
  var_e = exp((2*meanx)+(sd^2))*(exp(sd^2)-1)
  sd_sc = sqrt(var_e)
  
  Subject <- rep(1:N, each = ni)
  u_ijs = rt(N, df=dfi)
  u_ij = (u_ijs - 0)/sdi_sc
  u_i = rep(u_ij, each = ni)
  Replicates <- rep(1:ni, N)
  e_ij = (rlnorm(N*ni, mean=meanx, sd=sd) - exp(meanx + ((sd^2)/2)))/sd_sc
  beta_01 <- alpha_0 + qlnorm(tau1, mean=meanx, sd=sd)
  beta_02 <- alpha_0 + qlnorm(tau2, mean=meanx, sd=sd)
  Measure = alpha_0 + e_ij + u_i   
  
  true_width <- qlnorm(tau2, mean=meanx, sd=sd) - qlnorm(tau1, mean=meanx, sd=sd)
  
 
  data_model_1 = data.frame(Subject, u_i, e_ij, Replicates, Measure, beta_01, beta_02, true_width)
  return(list(d=data_model_1, u=u_ij, s=sdi_sc))
}

library(lqmm)
inp <- read.csv("sim_input.csv", sep=";", header=T)
inp$B[11] <- 100

ss <- lqmm_cal(B=inp$B[11], N=inp$N[11], ni=inp$ni[11], alpha_0=inp$alpha_0[11], meanx=inp$meanx[11], meanij=inp$meanij[11], 
               sd=inp$sd[11], sdi=inp$sdi[11], df=inp$df[11], dfi=inp$dfi[11], tau1=inp$tau1[11], tau2=inp$tau2[11], FUN = sim_chi_sc)
ss
ss$est1$beta0

ss <- lqmm_cal(B=inp$B[1], N=inp$N[1], ni=inp$ni[1], alpha_0=inp$alpha_0[1], meanx=inp$meanx[1], meanij=inp$meanij[1], 
               sd=inp$sd[1], sdi=inp$sdi[1], df=inp$df[1], dfi=inp$dfi[1], tau1=inp$tau1[1], tau2=inp$tau2[1], FUN = sim_lnorm_sc)
ss

(1.8 - mean(ss$estimates$beta01))/sd(ss$estimates$beta01)
qnorm(0.025, 0, 1)
r1 <- sim_chi(N=inp$N[11], ni=inp$ni[11], alpha_0=inp$alpha_0[11], meanx=inp$meanx[11], meanij=inp$meanij[11], 
              sd=inp$sd[11], sdi=inp$sdi[11],
              df=inp$df[11], dfi=inp$dfi[11], tau1=inp$tau1[11], tau2=inp$tau2[11])
r1$s

return(list(est = est, estimates = estimates, rem = rem, ranef1 = ranef1, ranef2 = ranef2,
            len = len_est, len_ran1 = len_ran1, len_ran2 = len_ran1))

}
