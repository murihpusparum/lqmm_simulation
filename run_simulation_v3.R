getwd()
path <- ("Y:\\Unit_MRG\\BREM&G\\Projecten\\1810018 I Am Pilot cohort\\05 Workdoc_Results\\WP4 Systems data analysis\\Statistical Analysis of Clinical Data\\PhD_Murih")
setwd(paste0(path, "\\03_Data_Analysis\\LQMM_Simulation\\Data"))
setwd("D:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\15.04.20")

a <- list()
b <- 1
b <- c(1,2)
b <- as.data.frame(b)
c <- c(1:3)
b <- as.data.frame(c)
a<-append(a,b)
a$b
a$c
data.frame(a)

inp <- read.csv("sim_input.csv", sep=";", header=T)
MSE_beta0_0.025 <- NULL
MSE_b0_0.025    <- NULL
MSE_beta0_0.975 <- NULL
MSE_b0_0.975    <- NULL
rel_width       <- NULL
B1 <- NULL
relbias_beta0_tau0.025 <- NULL
relbias_beta0_tau0.975 <- NULL
est11 <- list()
est21 <- list()
est12 <- list()
est22 <- list()
est13 <- list()
est23 <- list()
est14 <- list()
est24 <- list()
est15 <- list()
est25 <- list()
est16 <- list()
est26 <- list()
est17 <- list()
est27 <- list()
est18 <- list()
est28 <- list()
est19 <- list()
est29 <- list()
est110 <- list()
est210 <- list()
est111 <- list()
est211 <- list()
est112 <- list()
est212 <- list()
est113 <- list()
est213 <- list()
est114 <- list()
est214 <- list()
est115 <- list()
est215 <- list()

start_time <- Sys.time()
for (i in 1:2) {
  r1 <- lqmm_cal(B=inp$B[i], N=inp$N[i], ni=inp$ni[i], alpha_0=inp$alpha_0[i], meanx=inp$meanx[i], meanij=inp$meanij[i], 
                 sd=inp$sd[i], sdi=inp$sdi[i], df=inp$df[i], dfi=inp$dfi[i], tau1=inp$tau1[i], tau2=inp$tau2[i], FUN = sim_N)
  MSE_beta0_0.025[i]  = r1$MSE_tau_0.025$beta0
  MSE_b0_0.025[i]     = r1$MSE_tau_0.025$b0
  MSE_beta0_0.975[i]  = r1$MSE_tau_0.975$beta0
  MSE_b0_0.975[i]     = r1$MSE_tau_0.975$b0
  rel_width[i]        = r1$rel_width
  B1[i]  = r1$B1
  relbias_beta0_tau0.025[i] = r1$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i] = r1$relbias_beta0_tau0.975
  
  r1$est$beta01 <- data.frame(r1$est$beta01)
  r1$est$beta02 <- data.frame(r1$est$beta02)
  est11 <- append(est11, r1$est$beta01)
  est21 <- append(est21, r1$est$beta02)
#  est11 <- data.frame(est11)
#  est21 <- data.frame(est21)
  
  
  r2 <- lqmm_cal(B=inp$B[i+2], N=inp$N[i+2], ni=inp$ni[i+2], alpha_0=inp$alpha_0[i+2], meanx=inp$meanx[i+2], meanij=inp$meanij[i+2],
              sd=inp$sd[i+2], sdi=inp$sdi[i+2], df=inp$df[i+2], dfi=inp$dfi[i+2], tau1=inp$tau1[i+2], tau2=inp$tau2[i+2], FUN = sim_N_t)
  MSE_beta0_0.025[i+2]  = r2$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+2]     = r2$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+2]  = r2$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+2]     = r2$MSE_tau_0.975$b0
  rel_width[i+2]        = r2$rel_width
  B1[i+2]  = r2$B1
  relbias_beta0_tau0.025[i+2] = r2$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+2] = r2$relbias_beta0_tau0.975
  
  r2$est$beta01 <- data.frame(r2$est$beta01)
  r2$est$beta02 <- data.frame(r2$est$beta02)
  est12 <- append(est12, r2$est$beta01)
  est22 <- append(est22, r2$est$beta02)
# est12 <- data.frame(est12)
# est22 <- data.frame(est22)

  
  r3 <- lqmm_cal(B=inp$B[i+4], N=inp$N[i+4], ni=inp$ni[i+4], alpha_0=inp$alpha_0[i+4], meanx=inp$meanx[i+4], meanij=inp$meanij[i+4],
                 sd=inp$sd[i+4], sdi=inp$sdi[i+4], df=inp$df[i+4], dfi=inp$dfi[i+4], tau1=inp$tau1[i+4], tau2=inp$tau2[i+4], FUN = sim_t_sc)
  MSE_beta0_0.025[i+4]  = r3$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+4]     = r3$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+4]  = r3$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+4]     = r3$MSE_tau_0.975$b0
  rel_width[i+4]        = r3$rel_width
  B1[i+4]  = r3$B1
  relbias_beta0_tau0.025[i+4] = r3$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+4] = r3$relbias_beta0_tau0.975
  
  r3$est$beta01 <- data.frame(r3$est$beta01)
  r3$est$beta02 <- data.frame(r3$est$beta02)
  est13 <- append(est13, r3$est$beta01)
  est23 <- append(est23, r3$est$beta02)
#  est13 <- data.frame(est13)
#  est23 <- data.frame(est23)
  
  

  r4 <- lqmm_cal(B=inp$B[i+6], N=inp$N[i+6], ni=inp$ni[i+6], alpha_0=inp$alpha_0[i+6], meanx=inp$meanx[i+6], meanij=inp$meanij[i+6],
                 sd=inp$sd[i+6], sdi=inp$sdi[i+6], df=inp$df[i+6], dfi=inp$dfi[i+6], tau1=inp$tau1[i+6], tau2=inp$tau2[i+6], FUN = sim_t_sc)
  MSE_beta0_0.025[i+6]  = r4$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+6]     = r4$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+6]  = r4$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+6]     = r4$MSE_tau_0.975$b0
  rel_width[i+6]        = r4$rel_width
  B1[i+6]  = r4$B1
  relbias_beta0_tau0.025[i+6] = r4$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+6] = r4$relbias_beta0_tau0.975

  r4$est$beta01 <- data.frame(r4$est$beta01)
  r4$est$beta02 <- data.frame(r4$est$beta02)
  est14 <- append(est14, r4$est$beta01)
  est24 <- append(est24, r4$est$beta02)
#  est14 <- data.frame(est14)
#  est24 <- data.frame(est24)
  
  
  
  
  r5 <- lqmm_cal(B=inp$B[i+8], N=inp$N[i+8], ni=inp$ni[i+8], alpha_0=inp$alpha_0[i+8], meanx=inp$meanx[i+8], meanij=inp$meanij[i+8],
                 sd=inp$sd[i+8], sdi=inp$sdi[i+8], df=inp$df[i+8], dfi=inp$dfi[i+8], tau1=inp$tau1[i+8], tau2=inp$tau2[i+8], FUN = sim_t_N)
  MSE_beta0_0.025[i+8]  = r5$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+8]     = r5$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+8]  = r5$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+8]     = r5$MSE_tau_0.975$b0
  rel_width[i+8]        = r5$rel_width
  B1[i+8]  = r5$B1
  relbias_beta0_tau0.025[i+8] = r5$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+8] = r5$relbias_beta0_tau0.975

  r5$est$beta01 <- data.frame(r5$est$beta01)
  r5$est$beta02 <- data.frame(r5$est$beta02)
  est15 <- append(est15, r5$est$beta01)
  est25 <- append(est25, r5$est$beta02)
#  est15 <- data.frame(est15)
#  est25 <- data.frame(est25)
  
  

  r6 <- lqmm_cal(B=inp$B[i+10], N=inp$N[i+10], ni=inp$ni[i+10], alpha_0=inp$alpha_0[i+10], meanx=inp$meanx[i+10], meanij=inp$meanij[i+10],
                 sd=inp$sd[i+10], sdi=inp$sdi[i+10], df=inp$df[i+10], dfi=inp$dfi[i+10], tau1=inp$tau1[i+10], tau2=inp$tau2[i+10], FUN = sim_chi_sc)
  MSE_beta0_0.025[i+10]  = r6$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+10]     = r6$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+10]  = r6$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+10]     = r6$MSE_tau_0.975$b0
  rel_width[i+10]        = r6$rel_width
  B1[i+10]  = r6$B1
  relbias_beta0_tau0.025[i+10] = r6$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+10] = r6$relbias_beta0_tau0.975

  r6$est$beta01 <- data.frame(r6$est$beta01)
  r6$est$beta02 <- data.frame(r6$est$beta02)
  est16 <- append(est16, r6$est$beta01)
  est26 <- append(est26, r6$est$beta02)
#  est16 <- data.frame(est16)
#  est26 <- data.frame(est26)
  
 
  r7 <- lqmm_cal(B=inp$B[i+12], N=inp$N[i+12], ni=inp$ni[i+12], alpha_0=inp$alpha_0[i+12], meanx=inp$meanx[i+12], meanij=inp$meanij[i+12],
                 sd=inp$sd[i+12], sdi=inp$sdi[i+12], df=inp$df[i+12], dfi=inp$dfi[i+12], tau1=inp$tau1[i+12], tau2=inp$tau2[i+12], FUN = sim_N_chi)
  MSE_beta0_0.025[i+12]  = r7$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+12]     = r7$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+12]  = r7$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+12]     = r7$MSE_tau_0.975$b0
  rel_width[i+12]        = r7$rel_width
  B1[i+12]  = r7$B1
  relbias_beta0_tau0.025[i+12] = r7$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+12] = r7$relbias_beta0_tau0.975

  r7$est$beta01 <- data.frame(r7$est$beta01)
  r7$est$beta02 <- data.frame(r7$est$beta02)
  est17 <- append(est17, r7$est$beta01)
  est27 <- append(est27, r7$est$beta02)
#  est17 <- data.frame(est17)
#  est27 <- data.frame(est27)
  
  
 
  r8 <- lqmm_cal(B=inp$B[i+14], N=inp$N[i+14], ni=inp$ni[i+14], alpha_0=inp$alpha_0[i+14], meanx=inp$meanx[i+14], meanij=inp$meanij[i+14],
                 sd=inp$sd[i+14], sdi=inp$sdi[i+14], df=inp$df[i+14], dfi=inp$dfi[i+14], tau1=inp$tau1[i+14], tau2=inp$tau2[i+14], FUN = sim_chi_N)
  MSE_beta0_0.025[i+14]  = r8$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+14]     = r8$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+14]  = r8$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+14]     = r8$MSE_tau_0.975$b0
  rel_width[i+14]        = r8$rel_width
  B1[i+14]  = r8$B1
  relbias_beta0_tau0.025[i+14] = r8$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+14] = r8$relbias_beta0_tau0.975

  r8$est$beta01 <- data.frame(r8$est$beta01)
  r8$est$beta02 <- data.frame(r8$est$beta02)
  est18 <- append(est18, r8$est$beta01)
  est28 <- append(est28, r8$est$beta02)
#  est18 <- data.frame(est18)
#  est28 <- data.frame(est28)
  
  
  
  r9 <- lqmm_cal(B=inp$B[i+16], N=inp$N[i+16], ni=inp$ni[i+16], alpha_0=inp$alpha_0[i+16], meanx=inp$meanx[i+16], meanij=inp$meanij[i+16],
                 sd=inp$sd[i+16], sdi=inp$sdi[i+16], df=inp$df[i+16], dfi=inp$dfi[i+16], tau1=inp$tau1[i+16], tau2=inp$tau2[i+16], FUN = sim_N_lnorm)
  MSE_beta0_0.025[i+16]  = r9$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+16]     = r9$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+16]  = r9$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+16]     = r9$MSE_tau_0.975$b0
  rel_width[i+16]        = r9$rel_width
  B1[i+16]  = r9$B1
  relbias_beta0_tau0.025[i+16] = r9$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+16] = r9$relbias_beta0_tau0.975

  r9$est$beta01 <- data.frame(r9$est$beta01)
  r9$est$beta02 <- data.frame(r9$est$beta02)
  est19 <- append(est19, r9$est$beta01)
  est29 <- append(est19, r9$est$beta02)
#  est19 <- data.frame(est19)
#  est29 <- data.frame(est29)
  
  
  
  r10 <- lqmm_cal(B=inp$B[i+18], N=inp$N[i+18], ni=inp$ni[i+18], alpha_0=inp$alpha_0[i+18], meanx=inp$meanx[i+18], meanij=inp$meanij[i+18],
                 sd=inp$sd[i+18], sdi=inp$sdi[i+18], df=inp$df[i+18], dfi=inp$dfi[i+18], tau1=inp$tau1[i+18], tau2=inp$tau2[i+18], FUN = sim_lnorm_t)
  MSE_beta0_0.025[i+18]  = r10$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+18]     = r10$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+18]  = r10$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+18]     = r10$MSE_tau_0.975$b0
  rel_width[i+18]        = r10$rel_width
  B1[i+18]  = r10$B1
  relbias_beta0_tau0.025[i+18] = r10$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+18] = r10$relbias_beta0_tau0.975
  
  r10$est$beta01 <- data.frame(r10$est$beta01)
  r10$est$beta02 <- data.frame(r10$est$beta02)
  est110 <- append(est110, r10$est$beta01)
  est210 <- append(est210, r10$est$beta02)
  #  est110 <- data.frame(est110)
  #  est210 <- data.frame(est210)
  
  
  
  r11 <- lqmm_cal(B=inp$B[i+20], N=inp$N[i+20], ni=inp$ni[i+20], alpha_0=inp$alpha_0[i+20], meanx=inp$meanx[i+20], meanij=inp$meanij[i+20],
                  sd=inp$sd[i+20], sdi=inp$sdi[i+20], df=inp$df[i+20], dfi=inp$dfi[i+20], tau1=inp$tau1[i+20], tau2=inp$tau2[i+20], FUN = sim_lnorm_sc)
  MSE_beta0_0.025[i+20]  = r11$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+20]     = r11$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+20]  = r11$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+20]     = r11$MSE_tau_0.975$b0
  rel_width[i+20]        = r11$rel_width
  B1[i+20]  = r11$B1
  relbias_beta0_tau0.025[i+20] = r11$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+20] = r11$relbias_beta0_tau0.975
  
  r11$est$beta01 <- data.frame(r11$est$beta01)
  r11$est$beta02 <- data.frame(r11$est$beta02)
  est111 <- append(est111, r11$est$beta01)
  est211 <- append(est211, r11$est$beta02)
  #  est111 <- data.frame(est111)
  #  est211 <- data.frame(est211)
  
  
  
  r12 <- lqmm_cal(B=inp$B[i+22], N=inp$N[i+22], ni=inp$ni[i+22], alpha_0=inp$alpha_0[i+22], meanx=inp$meanx[i+22], meanij=inp$meanij[i+22],
                  sd=inp$sd[i+22], sdi=inp$sdi[i+22], df=inp$df[i+22], dfi=inp$dfi[i+22], tau1=inp$tau1[i+22], tau2=inp$tau2[i+22], FUN = sim_t)
  MSE_beta0_0.025[i+22]  = r12$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+22]     = r12$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+22]  = r12$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+22]     = r12$MSE_tau_0.975$b0
  rel_width[i+22]        = r12$rel_width
  B1[i+22]  = r12$B1
  relbias_beta0_tau0.025[i+22] = r12$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+22] = r12$relbias_beta0_tau0.975
  
  r12$est$beta01 <- data.frame(r12$est$beta01)
  r12$est$beta02 <- data.frame(r12$est$beta02)
  est112 <- append(est112, r12$est$beta01)
  est212 <- append(est212, r12$est$beta02)
  #  est112 <- data.frame(est112)
  #  est212 <- data.frame(est212)
  
  
  
  r13 <- lqmm_cal(B=inp$B[i+24], N=inp$N[i+24], ni=inp$ni[i+24], alpha_0=inp$alpha_0[i+24], meanx=inp$meanx[i+24], meanij=inp$meanij[i+24],
                  sd=inp$sd[i+24], sdi=inp$sdi[i+24], df=inp$df[i+24], dfi=inp$dfi[i+24], tau1=inp$tau1[i+24], tau2=inp$tau2[i+24], FUN = sim_t)
  MSE_beta0_0.025[i+24]  = r13$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+24]     = r13$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+24]  = r13$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+24]     = r13$MSE_tau_0.975$b0
  rel_width[i+24]        = r13$rel_width
  B1[i+24]  = r13$B1
  relbias_beta0_tau0.025[i+24] = r13$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+24] = r13$relbias_beta0_tau0.975
  
  r13$est$beta01 <- data.frame(r13$est$beta01)
  r13$est$beta02 <- data.frame(r13$est$beta02)
  est113 <- append(est113, r13$est$beta01)
  est213 <- append(est213, r13$est$beta02)
  #  est113 <- data.frame(est113)
  #  est213 <- data.frame(est213)
  
  
  
  r14 <- lqmm_cal(B=inp$B[i+26], N=inp$N[i+26], ni=inp$ni[i+26], alpha_0=inp$alpha_0[i+26], meanx=inp$meanx[i+26], meanij=inp$meanij[i+26],
                  sd=inp$sd[i+26], sdi=inp$sdi[i+26], df=inp$df[i+26], dfi=inp$dfi[i+26], tau1=inp$tau1[i+26], tau2=inp$tau2[i+26], FUN = sim_chi)
  MSE_beta0_0.025[i+26]  = r14$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+26]     = r14$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+26]  = r14$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+26]     = r14$MSE_tau_0.975$b0
  rel_width[i+26]        = r14$rel_width
  B1[i+26]  = r14$B1
  relbias_beta0_tau0.025[i+26] = r14$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+26] = r14$relbias_beta0_tau0.975
  
  r14$est$beta01 <- data.frame(r14$est$beta01)
  r14$est$beta02 <- data.frame(r14$est$beta02)
  est114 <- append(est114, r14$est$beta01)
  est214 <- append(est214, r14$est$beta02)
  #  est114 <- data.frame(est114)
  #  est214 <- data.frame(est214)
  
  
  r15 <- lqmm_cal(B=inp$B[i+28], N=inp$N[i+28], ni=inp$ni[i+28], alpha_0=inp$alpha_0[i+28], meanx=inp$meanx[i+28], meanij=inp$meanij[i+28],
                  sd=inp$sd[i+28], sdi=inp$sdi[i+28], df=inp$df[i+28], dfi=inp$dfi[i+28], tau1=inp$tau1[i+28], tau2=inp$tau2[i+28], FUN = sim_lnorm)
  MSE_beta0_0.025[i+28]  = r15$MSE_tau_0.025$beta0
  MSE_b0_0.025[i+28]     = r15$MSE_tau_0.025$b0
  MSE_beta0_0.975[i+28]  = r15$MSE_tau_0.975$beta0
  MSE_b0_0.975[i+28]     = r15$MSE_tau_0.975$b0
  rel_width[i+28]        = r15$rel_width
  B1[i+28]  = r15$B1
  relbias_beta0_tau0.025[i+28] = r15$relbias_beta0_tau0.025
  relbias_beta0_tau0.975[i+28] = r15$relbias_beta0_tau0.975
  
  r15$est$beta01 <- data.frame(r15$est$beta01)
  r15$est$beta02 <- data.frame(r15$est$beta02)
  est115 <- append(est115, r15$est$beta01)
  est215 <- append(est215, r15$est$beta02)
  #  est115 <- data.frame(est115)
  #  est215 <- data.frame(est215)
  
  
  #  est1[,i+28] = r15$est$beta01
  #  est2[,i+28] = r15$est$beta02
  
}

#result <- data.frame(MSE_beta0_0.025, MSE_b0_0.025, MSE_beta0_0.975, MSE_b0_0.975, rel_width)
inp2 <- inp
inp2 <- inp[-30,]

inp2$MSE_beta0_0.025 <- MSE_beta0_0.025[1:22]
inp2$MSE_b0_0.025    <- MSE_b0_0.025[1:22]
inp2$MSE_beta0_0.975 <- MSE_beta0_0.975[1:22]
inp2$MSE_b0_0.975    <- MSE_b0_0.975[1:22]
inp2$rel_width       <- rel_width[1:22]
inp2$B1              <- B1[1:22]
inp2$relbias_beta0_tau0.025 <- relbias_beta0_tau0.025[1:22]
inp2$relbias_beta0_tau0.975 <- relbias_beta0_tau0.975[1:22]

inp2$MSE_beta0_0.025 <- MSE_beta0_0.025
inp2$MSE_b0_0.025    <- MSE_b0_0.025
inp2$MSE_beta0_0.975 <- MSE_beta0_0.975
inp2$MSE_b0_0.975    <- MSE_b0_0.975
inp2$rel_width       <- rel_width
inp2$B1              <- B1
inp2$relbias_beta0_tau0.025 <- relbias_beta0_tau0.025
inp2$relbias_beta0_tau0.975 <- relbias_beta0_tau0.975

setwd("D:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\15.04.20\\beta_est\\sim1")
getwd()
length(est11[[1]])
length(est11[[2]])

est11 <- list()
write.csv(est11[[1]], "1tau25_30.csv")
write.csv(est11[[2]], "1tau25_100.csv")
est21 <- list()
write.csv(est21[[1]], "1tau975_30.csv")
write.csv(est21[[2]], "1tau975_100.csv")

est12 <- list()
write.csv(est12[[1]], "2tau25_30.csv")
write.csv(est12[[2]], "2tau25_100.csv")
est22 <- list()
write.csv(est22[[1]], "2tau975_30.csv")
write.csv(est22[[2]], "2tau975_100.csv")

est13 <- list()
write.csv(est13[[1]], "3tau25_30.csv")
write.csv(est13[[2]], "3tau25_100.csv")
est23 <- list()
write.csv(est23[[1]], "3tau975_30.csv")
write.csv(est23[[2]], "3tau975_100.csv")

est14 <- list()
write.csv(est14[[1]], "4tau25_30.csv")
write.csv(est14[[2]], "4tau25_100.csv")
est24 <- list()
write.csv(est24[[1]], "4tau975_30.csv")
write.csv(est24[[2]], "4tau975_100.csv")

est15 <- list()
write.csv(est15[[1]], "5tau25_30.csv")
write.csv(est15[[2]], "5tau25_100.csv")
est25 <- list()
write.csv(est25[[1]], "5tau975_30.csv")
write.csv(est25[[2]], "5tau975_100.csv")

est16 <- list()
write.csv(est16[[1]], "6tau25_30.csv")
write.csv(est16[[2]], "6tau25_100.csv")
est26 <- list()
write.csv(est26[[1]], "6tau975_30.csv")
write.csv(est26[[2]], "6tau975_100.csv")

est17 <- list()
write.csv(est17[[1]], "7tau25_30.csv")
write.csv(est17[[2]], "7tau25_100.csv")
est27 <- list()
write.csv(est27[[1]], "7tau975_30.csv")
write.csv(est27[[2]], "7tau975_100.csv")

est18 <- list()
write.csv(est18[[1]], "8tau25_30.csv")
write.csv(est18[[2]], "8tau25_100.csv")
est28 <- list()
write.csv(est28[[1]], "8tau975_30.csv")
write.csv(est28[[2]], "8tau975_100.csv")

est19 <- list()
write.csv(est19[[1]], "9tau25_30.csv")
write.csv(est19[[2]], "9tau25_100.csv")
est29 <- list()
write.csv(est29[[1]], "9tau975_30.csv")
write.csv(est29[[2]], "9tau975_100.csv")

est110 <- list()
write.csv(est110[[1]], "10tau25_30.csv")
write.csv(est110[[2]], "10tau25_100.csv")
est210 <- list()
write.csv(est210[[1]], "10tau975_30.csv")
write.csv(est210[[2]], "10tau975_100.csv")

est111 <- list()
write.csv(est111[[1]], "11tau25_30.csv")
write.csv(est111[[2]], "11tau25_100.csv")
est211 <- list()
write.csv(est211[[1]], "11tau975_30.csv")
write.csv(est211[[2]], "11tau975_100.csv")

est112 <- list()
write.csv(est112[[1]], "12tau25_30.csv")
#write.csv(est112[[2]], "12tau25_100.csv") no estimates - NA
est212 <- list()
write.csv(est212[[1]], "12tau975_30.csv")
#write.csv(est212[[2]], "12tau975_100.csv") no estimates - NA

est113 <- list()
write.csv(est113[[1]], "13tau25_30.csv")
#write.csv(est113[[2]], "13tau25_100.csv") no estimates - NA
est213 <- list()
write.csv(est213[[1]], "13tau975_30.csv")
#write.csv(est213[[2]], "13tau975_100.csv") no estimates - NA

est114 <- list()
write.csv(est114[[1]], "14tau25_30.csv")
#write.csv(est114[[2]], "14tau25_100.csv") no estimates - NA
est214 <- list()
write.csv(est214[[1]], "14tau975_30.csv")
#write.csv(est214[[2]], "14tau975_100.csv") no estimates - NA

est115 <- list()
write.csv(est115[[1]], "15tau25_30.csv")
#write.csv(est115[[2]], "15tau25_100.csv") no estimates - NA
est215 <- list()
write.csv(est215[[1]], "15tau975_30.csv")
#write.csv(est215[[2]], "15tau975_100.csv") no estimates - NA





#colnames(est1) <- inp$Dist
#colnames(est2) <- inp$Dist

write.csv(inp2, "16.04.20 - lqmm_sim_result1.csv")
end_time - start_time


write.csv(est1, "15.04.20 - lower_est3.csv")
write.csv(est2, "15.04.20 - upper_est3.csv")
end_time <- Sys.time()

getwd()
est214 <- as.data.frame(est214)
est11 <- data.frame(est11)

write.csv(est215)