setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
input1 <- read.csv("input1.csv", header=T, sep=",")

MSE <- data.frame(Setting=character(), MSE_beta0=integer(), MSE_b0=integer(), 
                  MSE_mean_ui=integer(), stringsAsFactors=FALSE)
rel_width <- NULL
#run simulation
for (i in seq(from=1, to=length(input1$B), by=2)) {
  result1 <- lqmm_int_eulnorm(B=input1[i, ]$B, N=input1[i, ]$N, ni=input1[i, ]$ni, beta_0 =input1[i,]$beta_0, 
                            sdi =input1[i, ]$sdi, sd=input1[i, ]$sd, tau=input1[i, ]$tau)
  result2 <- lqmm_int_eulnorm(B=input1[i+1, ]$B, N=input1[i+1, ]$N, ni=input1[i+1, ]$ni, beta_0 =input1[i+1,]$beta_0, 
                            sdi =input1[i+1, ]$sdi, sd=input1[i+1, ]$sd, tau=input1[i+1, ]$tau)
  
  MSE[i,] <- c(paste(paste0("N", input1[i, ]$N),
                     paste0("ni", input1[i, ]$ni),
                     paste0("tau", input1[i, ]$tau), sep="_"), 
               result1$MSE$beta0, result1$MSE$b0, result1$MSE$mean_ui)
  MSE[i+1,] <- c(paste(paste0("N", input1[i+1, ]$N),
                       paste0("ni", input1[i+1, ]$ni),
                       paste0("tau", input1[i+1, ]$tau), sep="_"), 
                 result2$MSE$beta0, result2$MSE$b0, result2$MSE$mean_ui)
  
  rel_width[i] <- relative_width(data1 = result1$ind_ci, data2 = result2$ind_ci, 
                                 B = input1[i+1, ]$B, N = input1[i+1, ]$N)
}

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
input1 <- read.csv("input2.csv", header=T, sep=",")
MSE <- data.frame(Setting=character(), MSE_beta0=integer(), MSE_b0=integer(), 
                  MSE_mean_ui=integer(), stringsAsFactors=FALSE)
rel_width <- NULL

#run simulation
for (i in seq(from=1, to=length(input1$B), by=2)) {
  result1 <- lqmm_int_echisq_ut(B=input1[i, ]$B, N=input1[i, ]$N, ni=input1[i, ]$ni, beta_0 =input1[i,]$beta_0, 
                             df = input1[i, ]$dfi, dfi=input1[i, ]$df1, tau=input1[i, ]$tau)
  result2 <- lqmm_int_echisq_ut(B=input1[i+1, ]$B, N=input1[i+1, ]$N, ni=input1[i+1, ]$ni, beta_0 =input1[i+1,]$beta_0, 
                             df = input1[i+1, ]$dfi, dfi=input1[i+1, ]$df1, tau=input1[i+1, ]$tau)
  
  MSE[i,] <- c(paste(paste0("N", input1[i, ]$N),
                     paste0("ni", input1[i, ]$ni),
                     paste0("tau", input1[i, ]$tau), sep="_"), 
               result1$MSE$beta0, result1$MSE$b0, result1$MSE$mean_ui)
  MSE[i+1,] <- c(paste(paste0("N", input1[i+1, ]$N),
                       paste0("ni", input1[i+1, ]$ni),
                       paste0("tau", input1[i+1, ]$tau), sep="_"), 
                 result2$MSE$beta0, result2$MSE$b0, result2$MSE$mean_ui)
  
  rel_width[i] <- relative_width(data1 = result1$ind_ci, data2 = result2$ind_ci, 
                                 B = input1[i+1, ]$B, N = input1[i+1, ]$N)
}

#remove null
rel_width <- rel_width[!sapply(rel_width, is.null)]
names(rel_width) <- paste("ci", seq_along(rel_width), sep="")
#wide format 30
rwidth_30 <- data.frame(rel_width$ci1, rel_width$ci2, rel_width$ci5, rel_width$ci6)
#long format 30
colnames(rwidth_30) <- c("N30_ni9_ci95", "N30_ni18_ci95", "N30_ni9_ci90", "N30_ni18_ci90")
rwidth_30_long <- data.frame(c(rel_width$ci1, rel_width$ci2, rel_width$ci5, rel_width$ci6),
                             c(rep("N30_ni9_ci95", 30),
                               rep("N30_ni18_ci95", 30),
                               rep("N30_ni9_ci90", 30),
                               rep("N30_ni18_ci90", 30)))
colnames(rwidth_30_long) <- c("ci_width", "setting")

#wide format 30
rwidth_100 <- data.frame(rel_width$ci3, rel_width$ci4, rel_width$ci7, rel_width$ci8)
#long format 30
colnames(rwidth_100) <- c("N100_ni9_ci95", "N100_ni18_ci95", "N100_ni9_ci90", "N100_ni18_ci90")
rwidth_100_long <- data.frame(c(rel_width$ci3, rel_width$ci4, rel_width$ci7, rel_width$ci8),
                             c(rep("N100_ni9_ci95", 100),
                               rep("N100_ni18_ci95", 100),
                               rep("N100_ni9_ci90", 100),
                               rep("N100_ni18_ci90", 100)))
colnames(rwidth_100_long) <- c("ci_width", "setting")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Result CSV 2")
write.csv(MSE, "echisq_ut_MSE.csv")
write.csv(rwidth_30, "echisq_ut_rwidth_30.csv")
write.csv(rwidth_30_long, "echisq_ut_rwidth_30_long.csv")
write.csv(rwidth_100, "echisq_ut_rwidth_100.csv")
write.csv(rwidth_100_long, "echisq_ut_rwidth_100_long.csv")
