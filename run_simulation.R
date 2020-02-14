##set working directory
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")

##clear global environment
rm(list = ls())
ls()

##call all libraries
library(lqmm)
library(lme4)
library(ggplot2)
library(directlabels)

##set source of functions
source("./lqmm_simulation_function.R")

#call the function
###################################################
#-------------E ~ N(0,1), U ~ N(0,1)--------------#
###################################################
#######################################
#                                     #
#   error_i ~ N(0,1), u_i ~ N(0,1)    #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
norm1 <- lqmm_int_enorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
norm2 <- lqmm_int_enorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
norm3 <- lqmm_int_enorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
norm4 <- lqmm_int_enorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
norm5 <- lqmm_int_enorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)

#######################################
#                                     #
#   error_i ~ N(0,1), u_i ~ N(0,1)    #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
norm6 <- lqmm_int_enorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
norm7 <- lqmm_int_enorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
norm8 <- lqmm_int_enorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
norm9 <- lqmm_int_enorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
norm10 <- lqmm_int_enorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)

#######################################
#                                     #
#   error_i ~ N(0,1), u_i ~ N(0,1)    #
#             tau = 0.5              #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
norm11 <- lqmm_int_enorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
norm12 <- lqmm_int_enorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
norm13 <- lqmm_int_enorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
norm14 <- lqmm_int_enorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
norm15 <- lqmm_int_enorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)

#######################################
#                                     #
#   error_i ~ N(0,1), u_i ~ N(0,1)    #
#             tau = 0.975              #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
norm16 <- lqmm_int_enorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
norm17 <- lqmm_int_enorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
norm18 <- lqmm_int_enorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
norm19 <- lqmm_int_enorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
norm20 <- lqmm_int_enorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)

### FOR N=30 ###
bias_ui_norm_N30 <- data.frame(norm1$results_ranef$bias_i, norm2$results_ranef$bias_i,
                               norm6$results_ranef$bias_i, norm7$results_ranef$bias_i,
                               norm11$results_ranef$bias_i, norm12$results_ranef$bias_i,
                               norm16$results_ranef$bias_i, norm17$results_ranef$bias_i)
colnames(bias_ui_norm_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                "ni9_t0.25", "ni18_t0.25",
                                "ni9_t0.5", "ni18_t0.5",
                                "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_norm_N30, "bias_ui_norm_N30.csv")
head(bias_ui_norm_N30)

bias_ui_norm_N30_long <- reshape(bias_ui_norm_N30,
                 direction = "long",
                 varying = list(names(bias_ui_norm_N30)[1:8]),
                 v.names = "bias",
                 timevar = "Setting",
                 times = c("ni9_t0.05", "ni18_t0.05",
                           "ni9_t0.25", "ni18_t0.25",
                           "ni9_t0.5", "ni18_t0.5",
                           "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_norm_N100 <- data.frame(norm3$results_ranef$bias_i, norm4$results_ranef$bias_i,
                               norm8$results_ranef$bias_i, norm9$results_ranef$bias_i,
                               norm13$results_ranef$bias_i, norm14$results_ranef$bias_i,
                               norm18$results_ranef$bias_i, norm19$results_ranef$bias_i)
colnames(bias_ui_norm_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                "ni9_t0.25", "ni18_t0.25",
                                "ni9_t0.5", "ni18_t0.5",
                                "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_norm_N100, "bias_ui_norm_N100.csv")
head(bias_ui_norm_N100)

bias_ui_norm_N100_long <- reshape(bias_ui_norm_N100,
                                 direction = "long",
                                 varying = list(names(bias_ui_norm_N100)[1:8]),
                                 v.names = "bias",
                                 timevar = "Setting",
                                 times = c("ni9_t0.05", "ni18_t0.05",
                                           "ni9_t0.25", "ni18_t0.25",
                                           "ni9_t0.5", "ni18_t0.5",
                                           "ni9_t0.975", "ni18_t0.975"))

mean_bias_norm <- rbind(
cbind(norm1$results_ranef$mean_bias,
      norm1$bias_estimates$mean_bias_beta0,
      norm1$bias_estimates$mean_bias_b0),
cbind(norm2$results_ranef$mean_bias,
      norm2$bias_estimates$mean_bias_beta0,
      norm2$bias_estimates$mean_bias_b0),
cbind(norm3$results_ranef$mean_bias,
      norm3$bias_estimates$mean_bias_beta0,
      norm3$bias_estimates$mean_bias_b0),
cbind(norm4$results_ranef$mean_bias,
      norm4$bias_estimates$mean_bias_beta0,
      norm4$bias_estimates$mean_bias_b0),
cbind(norm5$results_ranef$mean_bias,
      norm5$bias_estimates$mean_bias_beta0,
      norm5$bias_estimates$mean_bias_b0),
cbind(norm6$results_ranef$mean_bias,
      norm6$bias_estimates$mean_bias_beta0,
      norm6$bias_estimates$mean_bias_b0),
cbind(norm7$results_ranef$mean_bias,
      norm7$bias_estimates$mean_bias_beta0,
      norm7$bias_estimates$mean_bias_b0),
cbind(norm8$results_ranef$mean_bias,
      norm8$bias_estimates$mean_bias_beta0,
      norm8$bias_estimates$mean_bias_b0),
cbind(norm9$results_ranef$mean_bias,
      norm9$bias_estimates$mean_bias_beta0,
      norm9$bias_estimates$mean_bias_b0),
cbind(norm10$results_ranef$mean_bias,
      norm10$bias_estimates$mean_bias_beta0,
      norm10$bias_estimates$mean_bias_b0),
cbind(norm11$results_ranef$mean_bias,
      norm11$bias_estimates$mean_bias_beta0,
      norm11$bias_estimates$mean_bias_b0),
cbind(norm12$results_ranef$mean_bias,
      norm12$bias_estimates$mean_bias_beta0,
      norm12$bias_estimates$mean_bias_b0),
cbind(norm13$results_ranef$mean_bias,
      norm13$bias_estimates$mean_bias_beta0,
      norm13$bias_estimates$mean_bias_b0),
cbind(norm14$results_ranef$mean_bias,
      norm14$bias_estimates$mean_bias_beta0,
      norm14$bias_estimates$mean_bias_b0),
cbind(norm15$results_ranef$mean_bias,
      norm15$bias_estimates$mean_bias_beta0,
      norm15$bias_estimates$mean_bias_b0),
cbind(norm16$results_ranef$mean_bias,
      norm16$bias_estimates$mean_bias_beta0,
      norm16$bias_estimates$mean_bias_b0),
cbind(norm17$results_ranef$mean_bias,
      norm17$bias_estimates$mean_bias_beta0,
      norm17$bias_estimates$mean_bias_b0),
cbind(norm18$results_ranef$mean_bias,
      norm18$bias_estimates$mean_bias_beta0,
      norm18$bias_estimates$mean_bias_b0),
cbind(norm19$results_ranef$mean_bias,
      norm19$bias_estimates$mean_bias_beta0,
      norm19$bias_estimates$mean_bias_b0),
cbind(norm20$results_ranef$mean_bias,
      norm20$bias_estimates$mean_bias_beta0,
      norm20$bias_estimates$mean_bias_b0))


rownames(mean_bias_norm) <- c("norm1", "norm2", "norm3", "norm4", "norm5", "norm6", "norm7",
                 "norm8", "norm9", "norm10", "norm11", "norm12", "norm13", "norm14", 
                 "norm15", "norm16", "norm17", "norm18", "norm19", "norm20")
colnames(mean_bias_norm) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_norm, "mean_bias_norm.csv")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")

png("Bias_ui_N30_norm.png", width = 600, height = 500)
ggplot(bias_ui_norm_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ N(0,1), e_ij ~ N(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_norm.png", width = 600, height = 500)
ggplot(bias_ui_norm_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ N(0,1), e_ij ~ N(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()


###################################################
#-----------E ~ LNORM(0,1), U ~ N(0,1)------------#
###################################################
#######################################
#                                     #
# error_i ~ LNORM(0,1), u_i ~ N(0,1)  #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
elnorm1 <- lqmm_int_elnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
elnorm2 <- lqmm_int_elnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
elnorm3 <- lqmm_int_elnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
elnorm4 <- lqmm_int_elnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
elnorm5 <- lqmm_int_elnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)

#######################################
#                                     #
# error_i ~ LNORM(0,1), u_i ~ N(0,1)  #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
elnorm6 <- lqmm_int_elnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
elnorm7 <- lqmm_int_elnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
elnorm8 <- lqmm_int_elnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
elnorm9 <- lqmm_int_elnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
elnorm10 <- lqmm_int_elnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)

#######################################
#                                     #
# error_i ~ LNORM(0,1), u_i ~ N(0,1)  #
#             tau = 0.5               #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
elnorm11 <- lqmm_int_elnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
elnorm12 <- lqmm_int_elnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
elnorm13 <- lqmm_int_elnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
elnorm14 <- lqmm_int_elnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
elnorm15 <- lqmm_int_elnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)

#######################################
#                                     #
#  error_i ~ LNORM(0,1), u_i ~ N(0,1) #
#             tau = 0.975             #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
elnorm16 <- lqmm_int_elnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
elnorm17 <- lqmm_int_elnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
elnorm18 <- lqmm_int_elnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
elnorm19 <- lqmm_int_elnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
elnorm20 <- lqmm_int_elnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)

### FOR N=30 ###
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")

bias_ui_elnorm_N30 <- data.frame(elnorm1$results_ranef$bias_i, elnorm2$results_ranef$bias_i,
                                   elnorm6$results_ranef$bias_i, elnorm7$results_ranef$bias_i,
                                   elnorm11$results_ranef$bias_i, elnorm12$results_ranef$bias_i,
                                   elnorm16$results_ranef$bias_i, elnorm17$results_ranef$bias_i)
colnames(bias_ui_elnorm_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                  "ni9_t0.25", "ni18_t0.25",
                                  "ni9_t0.5", "ni18_t0.5",
                                  "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_elnorm_N30, "bias_ui_elnorm_N30.csv")
head(bias_ui_elnorm_N30)

bias_ui_elnorm_N30_long <- reshape(bias_ui_elnorm_N30,
                                   direction = "long",
                                   varying = list(names(bias_ui_elnorm_N30)[1:8]),
                                   v.names = "bias",
                                   timevar = "Setting",
                                   times = c("ni9_t0.05", "ni18_t0.05",
                                             "ni9_t0.25", "ni18_t0.25",
                                             "ni9_t0.5", "ni18_t0.5",
                                             "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_elnorm_N100 <- data.frame(elnorm3$results_ranef$bias_i, elnorm4$results_ranef$bias_i,
                                  elnorm8$results_ranef$bias_i, elnorm9$results_ranef$bias_i,
                                  elnorm13$results_ranef$bias_i, elnorm14$results_ranef$bias_i,
                                  elnorm18$results_ranef$bias_i, elnorm19$results_ranef$bias_i)
colnames(bias_ui_elnorm_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                   "ni9_t0.25", "ni18_t0.25",
                                   "ni9_t0.5", "ni18_t0.5",
                                   "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_elnorm_N100, "bias_ui_elnorm_N100.csv")
head(bias_ui_elnorm_N100)

bias_ui_elnorm_N100_long <- reshape(bias_ui_elnorm_N100,
                                    direction = "long",
                                    varying = list(names(bias_ui_elnorm_N100)[1:8]),
                                    v.names = "bias",
                                    timevar = "Setting",
                                    times = c("ni9_t0.05", "ni18_t0.05",
                                              "ni9_t0.25", "ni18_t0.25",
                                              "ni9_t0.5", "ni18_t0.5",
                                              "ni9_t0.975", "ni18_t0.975"))

mean_bias_elnorm <- rbind(
  cbind(elnorm1$results_ranef$mean_bias,
        elnorm1$bias_estimates$mean_bias_beta0,
        elnorm1$bias_estimates$mean_bias_b0),
  cbind(elnorm2$results_ranef$mean_bias,
        elnorm2$bias_estimates$mean_bias_beta0,
        elnorm2$bias_estimates$mean_bias_b0),
  cbind(elnorm3$results_ranef$mean_bias,
        elnorm3$bias_estimates$mean_bias_beta0,
        elnorm3$bias_estimates$mean_bias_b0),
  cbind(elnorm4$results_ranef$mean_bias,
        elnorm4$bias_estimates$mean_bias_beta0,
        elnorm4$bias_estimates$mean_bias_b0),
  cbind(elnorm5$results_ranef$mean_bias,
        elnorm5$bias_estimates$mean_bias_beta0,
        elnorm5$bias_estimates$mean_bias_b0),
  cbind(elnorm6$results_ranef$mean_bias,
        elnorm6$bias_estimates$mean_bias_beta0,
        elnorm6$bias_estimates$mean_bias_b0),
  cbind(elnorm7$results_ranef$mean_bias,
        elnorm7$bias_estimates$mean_bias_beta0,
        elnorm7$bias_estimates$mean_bias_b0),
  cbind(elnorm8$results_ranef$mean_bias,
        elnorm8$bias_estimates$mean_bias_beta0,
        elnorm8$bias_estimates$mean_bias_b0),
  cbind(elnorm9$results_ranef$mean_bias,
        elnorm9$bias_estimates$mean_bias_beta0,
        elnorm9$bias_estimates$mean_bias_b0),
  cbind(elnorm10$results_ranef$mean_bias,
        elnorm10$bias_estimates$mean_bias_beta0,
        elnorm10$bias_estimates$mean_bias_b0),
  cbind(elnorm11$results_ranef$mean_bias,
        elnorm11$bias_estimates$mean_bias_beta0,
        elnorm11$bias_estimates$mean_bias_b0),
  cbind(elnorm12$results_ranef$mean_bias,
        elnorm12$bias_estimates$mean_bias_beta0,
        elnorm12$bias_estimates$mean_bias_b0),
  cbind(elnorm13$results_ranef$mean_bias,
        elnorm13$bias_estimates$mean_bias_beta0,
        elnorm13$bias_estimates$mean_bias_b0),
  cbind(elnorm14$results_ranef$mean_bias,
        elnorm14$bias_estimates$mean_bias_beta0,
        elnorm14$bias_estimates$mean_bias_b0),
  cbind(elnorm15$results_ranef$mean_bias,
        elnorm15$bias_estimates$mean_bias_beta0,
        elnorm15$bias_estimates$mean_bias_b0),
  cbind(elnorm16$results_ranef$mean_bias,
        elnorm16$bias_estimates$mean_bias_beta0,
        elnorm16$bias_estimates$mean_bias_b0),
  cbind(elnorm17$results_ranef$mean_bias,
        elnorm17$bias_estimates$mean_bias_beta0,
        elnorm17$bias_estimates$mean_bias_b0),
  cbind(elnorm18$results_ranef$mean_bias,
        elnorm18$bias_estimates$mean_bias_beta0,
        elnorm18$bias_estimates$mean_bias_b0),
  cbind(elnorm19$results_ranef$mean_bias,
        elnorm19$bias_estimates$mean_bias_beta0,
        elnorm19$bias_estimates$mean_bias_b0),
  cbind(elnorm20$results_ranef$mean_bias,
        elnorm20$bias_estimates$mean_bias_beta0,
        elnorm20$bias_estimates$mean_bias_b0))


rownames(mean_bias_elnorm) <- c("elnorm1", "elnorm2", "elnorm3", "elnorm4", "elnorm5", "elnorm6", "elnorm7",
                                "elnorm8", "elnorm9", "elnorm10", "elnorm11", "elnorm12", "elnorm13", "elnorm14", 
                                "elnorm15", "elnorm16", "elnorm17", "elnorm18", "elnorm19", "elnorm20")
colnames(mean_bias_elnorm) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_elnorm, "mean_bias_elnorm.csv")


setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")

png("Bias_ui_N30_elnorm.png", width = 600, height = 500)
ggplot(bias_ui_elnorm_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ N(0,1), e_ij ~ LNORM(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_elnorm.png", width = 600, height = 500)
ggplot(bias_ui_elnorm_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ N(0,1), e_ij ~ LNORM(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()



#--------------------------------------E ~ N(0,1), U ~ LNORM(0,1)





###################################################
#-----------E ~ N(0,1), U ~ LNORM(0,1)------------#
###################################################
#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ LNORM(0,1)  #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
ulnorm1 <- lqmm_int_ulnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
ulnorm2 <- lqmm_int_ulnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
ulnorm3 <- lqmm_int_ulnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
ulnorm4 <- lqmm_int_ulnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
ulnorm5 <- lqmm_int_ulnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)

#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ LNORM(0,1)  #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
ulnorm6 <- lqmm_int_ulnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
ulnorm7 <- lqmm_int_ulnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
ulnorm8 <- lqmm_int_ulnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
ulnorm9 <- lqmm_int_ulnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
ulnorm10 <- lqmm_int_ulnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)

#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ LNORM(0,1)  #
#             tau = 0.5               #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
ulnorm11 <- lqmm_int_ulnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
ulnorm12 <- lqmm_int_ulnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
ulnorm13 <- lqmm_int_ulnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
ulnorm14 <- lqmm_int_ulnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
ulnorm15 <- lqmm_int_ulnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)

#######################################
#                                     #
#  error_i ~ N(0,1), u_i ~ LNORM(0,1) #
#             tau = 0.975             #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
ulnorm16 <- lqmm_int_ulnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
ulnorm17 <- lqmm_int_ulnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
ulnorm18 <- lqmm_int_ulnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
ulnorm19 <- lqmm_int_ulnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
ulnorm20 <- lqmm_int_ulnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)


### FOR N=30 ###
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
bias_ui_ulnorm_N30 <- data.frame(ulnorm1$results_ranef$bias_i, ulnorm2$results_ranef$bias_i,
                                 ulnorm6$results_ranef$bias_i, ulnorm7$results_ranef$bias_i,
                                 ulnorm11$results_ranef$bias_i, ulnorm12$results_ranef$bias_i,
                                 ulnorm16$results_ranef$bias_i, ulnorm17$results_ranef$bias_i)
colnames(bias_ui_ulnorm_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                  "ni9_t0.25", "ni18_t0.25",
                                  "ni9_t0.5", "ni18_t0.5",
                                  "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_ulnorm_N30, "bias_ui_ulnorm_N30.csv")
head(bias_ui_ulnorm_N30)

bias_ui_ulnorm_N30_long <- reshape(bias_ui_ulnorm_N30,
                                   direction = "long",
                                   varying = list(names(bias_ui_ulnorm_N30)[1:8]),
                                   v.names = "bias",
                                   timevar = "Setting",
                                   times = c("ni9_t0.05", "ni18_t0.05",
                                             "ni9_t0.25", "ni18_t0.25",
                                             "ni9_t0.5", "ni18_t0.5",
                                             "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_ulnorm_N100 <- data.frame(ulnorm3$results_ranef$bias_i, ulnorm4$results_ranef$bias_i,
                                  ulnorm8$results_ranef$bias_i, ulnorm9$results_ranef$bias_i,
                                  ulnorm13$results_ranef$bias_i, ulnorm14$results_ranef$bias_i,
                                  ulnorm18$results_ranef$bias_i, ulnorm19$results_ranef$bias_i)
colnames(bias_ui_ulnorm_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                   "ni9_t0.25", "ni18_t0.25",
                                   "ni9_t0.5", "ni18_t0.5",
                                   "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_ulnorm_N100, "bias_ui_ulnorm_N100.csv")
head(bias_ui_ulnorm_N100)

bias_ui_ulnorm_N100_long <- reshape(bias_ui_ulnorm_N100,
                                    direction = "long",
                                    varying = list(names(bias_ui_ulnorm_N100)[1:8]),
                                    v.names = "bias",
                                    timevar = "Setting",
                                    times = c("ni9_t0.05", "ni18_t0.05",
                                              "ni9_t0.25", "ni18_t0.25",
                                              "ni9_t0.5", "ni18_t0.5",
                                              "ni9_t0.975", "ni18_t0.975"))


mean_bias_ulnorm <- rbind(
  cbind(ulnorm1$results_ranef$mean_bias,
        ulnorm1$bias_estimates$mean_bias_beta0,
        ulnorm1$bias_estimates$mean_bias_b0),
  cbind(ulnorm2$results_ranef$mean_bias,
        ulnorm2$bias_estimates$mean_bias_beta0,
        ulnorm2$bias_estimates$mean_bias_b0),
  cbind(ulnorm3$results_ranef$mean_bias,
        ulnorm3$bias_estimates$mean_bias_beta0,
        ulnorm3$bias_estimates$mean_bias_b0),
  cbind(ulnorm4$results_ranef$mean_bias,
        ulnorm4$bias_estimates$mean_bias_beta0,
        ulnorm4$bias_estimates$mean_bias_b0),
#  cbind(ulnorm5$results_ranef$mean_bias,
#        ulnorm5$bias_estimates$mean_bias_beta0,
#        ulnorm5$bias_estimates$mean_bias_b0),
  cbind(ulnorm6$results_ranef$mean_bias,
        ulnorm6$bias_estimates$mean_bias_beta0,
        ulnorm6$bias_estimates$mean_bias_b0),
  cbind(ulnorm7$results_ranef$mean_bias,
        ulnorm7$bias_estimates$mean_bias_beta0,
        ulnorm7$bias_estimates$mean_bias_b0),
  cbind(ulnorm8$results_ranef$mean_bias,
        ulnorm8$bias_estimates$mean_bias_beta0,
        ulnorm8$bias_estimates$mean_bias_b0),
  cbind(ulnorm9$results_ranef$mean_bias,
        ulnorm9$bias_estimates$mean_bias_beta0,
        ulnorm9$bias_estimates$mean_bias_b0),
#  cbind(ulnorm10$results_ranef$mean_bias,
#        ulnorm10$bias_estimates$mean_bias_beta0,
#        ulnorm10$bias_estimates$mean_bias_b0),
  cbind(ulnorm11$results_ranef$mean_bias,
        ulnorm11$bias_estimates$mean_bias_beta0,
        ulnorm11$bias_estimates$mean_bias_b0),
  cbind(ulnorm12$results_ranef$mean_bias,
        ulnorm12$bias_estimates$mean_bias_beta0,
        ulnorm12$bias_estimates$mean_bias_b0),
  cbind(ulnorm13$results_ranef$mean_bias,
        ulnorm13$bias_estimates$mean_bias_beta0,
        ulnorm13$bias_estimates$mean_bias_b0),
  cbind(ulnorm14$results_ranef$mean_bias,
        ulnorm14$bias_estimates$mean_bias_beta0,
        ulnorm14$bias_estimates$mean_bias_b0),
#  cbind(ulnorm15$results_ranef$mean_bias,
#        ulnorm15$bias_estimates$mean_bias_beta0,
#        ulnorm15$bias_estimates$mean_bias_b0),
  cbind(ulnorm16$results_ranef$mean_bias,
        ulnorm16$bias_estimates$mean_bias_beta0,
        ulnorm16$bias_estimates$mean_bias_b0),
  cbind(ulnorm17$results_ranef$mean_bias,
        ulnorm17$bias_estimates$mean_bias_beta0,
        ulnorm17$bias_estimates$mean_bias_b0),
  cbind(ulnorm18$results_ranef$mean_bias,
        ulnorm18$bias_estimates$mean_bias_beta0,
        ulnorm18$bias_estimates$mean_bias_b0),
  cbind(ulnorm19$results_ranef$mean_bias,
        ulnorm19$bias_estimates$mean_bias_beta0,
        ulnorm19$bias_estimates$mean_bias_b0))
#  cbind(ulnorm20$results_ranef$mean_bias,
#        ulnorm20$bias_estimates$mean_bias_beta0,
#        ulnorm20$bias_estimates$mean_bias_b0))


rownames(mean_bias_ulnorm) <- c("ulnorm1", "ulnorm2", "ulnorm3", "ulnorm4", 
                                "ulnorm6", "ulnorm7", "ulnorm8", "ulnorm9", 
                                "ulnorm11", "ulnorm12", "ulnorm13", "ulnorm14", 
                                "ulnorm16", "ulnorm17", "ulnorm18", "ulnorm19")

colnames(mean_bias_ulnorm) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_ulnorm, "mean_bias_ulnorm.csv")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")
png("Bias_ui_N30_ulnorm.png", width = 600, height = 500)
ggplot(bias_ui_ulnorm_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ LNORM(0,1), e_ij ~ N(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_ulnorm.png", width = 600, height = 500)
ggplot(bias_ui_ulnorm_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ LNORM(0,1), e_ij ~ N(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()




###################################################
#-----------E ~ LNORM(0,1), U ~ LNORM(0,1)------------#
###################################################
#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ LNORM(0,1)  #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
eulnorm1 <- lqmm_int_eulnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
eulnorm2 <- lqmm_int_eulnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
eulnorm3 <- lqmm_int_eulnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
eulnorm4 <- lqmm_int_eulnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
eulnorm5 <- lqmm_int_eulnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.05)

#######################################
#                                     #
# error_i ~ LNORM(0,1), u_i ~ LNORM(0,1)  #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
eulnorm6 <- lqmm_int_eulnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
eulnorm7 <- lqmm_int_eulnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
eulnorm8 <- lqmm_int_eulnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
eulnorm9 <- lqmm_int_eulnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
eulnorm10 <- lqmm_int_eulnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.25)

#######################################
#                                     #
# error_i ~ LNORM(0,1), u_i ~ LNORM(0,1)  #
#             tau = 0.5               #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
eulnorm11 <- lqmm_int_eulnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
eulnorm12 <- lqmm_int_eulnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
eulnorm13 <- lqmm_int_eulnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
eulnorm14 <- lqmm_int_eulnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
eulnorm15 <- lqmm_int_eulnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.5)

#######################################
#                                     #
#  error_i ~ LNORM(0,1), u_i ~ LNORM(0,1) #
#             tau = 0.975             #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
eulnorm16 <- lqmm_int_eulnorm(B=1000, N=30, ni=9, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
eulnorm17 <- lqmm_int_eulnorm(B=1000, N=30, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
eulnorm18 <- lqmm_int_eulnorm(B=1000, N=100, ni=9, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
eulnorm19 <- lqmm_int_eulnorm(B=1000, N=100, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
eulnorm20 <- lqmm_int_eulnorm(B=1000, N=500, ni=18, beta_0=2, sdi=1, sd=1, tau=0.975)

### FOR N=30 ###
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
bias_ui_eulnorm_N30 <- data.frame(eulnorm1$results_ranef$bias_i, eulnorm2$results_ranef$bias_i,
                                  eulnorm6$results_ranef$bias_i, eulnorm7$results_ranef$bias_i,
                                  eulnorm11$results_ranef$bias_i, eulnorm12$results_ranef$bias_i,
                                  eulnorm16$results_ranef$bias_i, eulnorm17$results_ranef$bias_i)
colnames(bias_ui_eulnorm_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                   "ni9_t0.25", "ni18_t0.25",
                                   "ni9_t0.5", "ni18_t0.5",
                                   "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_eulnorm_N30, "bias_ui_eulnorm_N30.csv")
head(bias_ui_eulnorm_N30)

bias_ui_eulnorm_N30_long <- reshape(bias_ui_eulnorm_N30,
                                    direction = "long",
                                    varying = list(names(bias_ui_eulnorm_N30)[1:8]),
                                    v.names = "bias",
                                    timevar = "Setting",
                                    times = c("ni9_t0.05", "ni18_t0.05",
                                              "ni9_t0.25", "ni18_t0.25",
                                              "ni9_t0.5", "ni18_t0.5",
                                              "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_eulnorm_N100 <- data.frame(eulnorm3$results_ranef$bias_i, eulnorm4$results_ranef$bias_i,
                                   eulnorm8$results_ranef$bias_i, eulnorm9$results_ranef$bias_i,
                                   eulnorm13$results_ranef$bias_i, eulnorm14$results_ranef$bias_i,
                                   eulnorm18$results_ranef$bias_i, eulnorm19$results_ranef$bias_i)
colnames(bias_ui_eulnorm_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                    "ni9_t0.25", "ni18_t0.25",
                                    "ni9_t0.5", "ni18_t0.5",
                                    "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_eulnorm_N100, "bias_ui_eulnorm_N100.csv")
head(bias_ui_eulnorm_N100)

bias_ui_eulnorm_N100_long <- reshape(bias_ui_eulnorm_N100,
                                     direction = "long",
                                     varying = list(names(bias_ui_eulnorm_N100)[1:8]),
                                     v.names = "bias",
                                     timevar = "Setting",
                                     times = c("ni9_t0.05", "ni18_t0.05",
                                               "ni9_t0.25", "ni18_t0.25",
                                               "ni9_t0.5", "ni18_t0.5",
                                               "ni9_t0.975", "ni18_t0.975"))

mean_bias_eulnorm <- rbind(
  cbind(eulnorm1$results_ranef$mean_bias,
        eulnorm1$bias_estimates$mean_bias_beta0,
        eulnorm1$bias_estimates$mean_bias_b0),
  cbind(eulnorm2$results_ranef$mean_bias,
        eulnorm2$bias_estimates$mean_bias_beta0,
        eulnorm2$bias_estimates$mean_bias_b0),
  cbind(eulnorm3$results_ranef$mean_bias,
        eulnorm3$bias_estimates$mean_bias_beta0,
        eulnorm3$bias_estimates$mean_bias_b0),
  cbind(eulnorm4$results_ranef$mean_bias,
        eulnorm4$bias_estimates$mean_bias_beta0,
        eulnorm4$bias_estimates$mean_bias_b0),
  #  cbind(eulnorm5$results_ranef$mean_bias,
  #        eulnorm5$bias_estimates$mean_bias_beta0,
  #        eulnorm5$bias_estimates$mean_bias_b0),
  cbind(eulnorm6$results_ranef$mean_bias,
        eulnorm6$bias_estimates$mean_bias_beta0,
        eulnorm6$bias_estimates$mean_bias_b0),
  cbind(eulnorm7$results_ranef$mean_bias,
        eulnorm7$bias_estimates$mean_bias_beta0,
        eulnorm7$bias_estimates$mean_bias_b0),
  cbind(eulnorm8$results_ranef$mean_bias,
        eulnorm8$bias_estimates$mean_bias_beta0,
        eulnorm8$bias_estimates$mean_bias_b0),
  cbind(eulnorm9$results_ranef$mean_bias,
        eulnorm9$bias_estimates$mean_bias_beta0,
        eulnorm9$bias_estimates$mean_bias_b0),
  #  cbind(eulnorm10$results_ranef$mean_bias,
  #        eulnorm10$bias_estimates$mean_bias_beta0,
  #        eulnorm10$bias_estimates$mean_bias_b0),
  cbind(eulnorm11$results_ranef$mean_bias,
        eulnorm11$bias_estimates$mean_bias_beta0,
        eulnorm11$bias_estimates$mean_bias_b0),
  cbind(eulnorm12$results_ranef$mean_bias,
        eulnorm12$bias_estimates$mean_bias_beta0,
        eulnorm12$bias_estimates$mean_bias_b0),
  cbind(eulnorm13$results_ranef$mean_bias,
        eulnorm13$bias_estimates$mean_bias_beta0,
        eulnorm13$bias_estimates$mean_bias_b0),
  cbind(eulnorm14$results_ranef$mean_bias,
        eulnorm14$bias_estimates$mean_bias_beta0,
        eulnorm14$bias_estimates$mean_bias_b0),
  #  cbind(eulnorm15$results_ranef$mean_bias,
  #        eulnorm15$bias_estimates$mean_bias_beta0,
  #        eulnorm15$bias_estimates$mean_bias_b0),
  cbind(eulnorm16$results_ranef$mean_bias,
        eulnorm16$bias_estimates$mean_bias_beta0,
        eulnorm16$bias_estimates$mean_bias_b0),
  cbind(eulnorm17$results_ranef$mean_bias,
        eulnorm17$bias_estimates$mean_bias_beta0,
        eulnorm17$bias_estimates$mean_bias_b0),
  cbind(eulnorm18$results_ranef$mean_bias,
        eulnorm18$bias_estimates$mean_bias_beta0,
        eulnorm18$bias_estimates$mean_bias_b0),
  cbind(eulnorm19$results_ranef$mean_bias,
        eulnorm19$bias_estimates$mean_bias_beta0,
        eulnorm19$bias_estimates$mean_bias_b0))
#  cbind(eulnorm20$results_ranef$mean_bias,
#        eulnorm20$bias_estimates$mean_bias_beta0,
#        eulnorm20$bias_estimates$mean_bias_b0))


rownames(mean_bias_eulnorm) <- c("eulnorm1", "eulnorm2", "eulnorm3", "eulnorm4", 
                                 "eulnorm6", "eulnorm7", "eulnorm8", "eulnorm9", 
                                 "eulnorm11", "eulnorm12", "eulnorm13", "eulnorm14", 
                                 "eulnorm16", "eulnorm17", "eulnorm18", "eulnorm19")

colnames(mean_bias_eulnorm) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_eulnorm, "mean_bias_eulnorm.csv")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")
png("Bias_ui_N30_eulnorm.png", width = 600, height = 500)
ggplot(bias_ui_eulnorm_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ LNORM(0,1), e_ij ~ LNORM(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_eulnorm.png", width = 600, height = 500)
ggplot(bias_ui_eulnorm_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ LNORM(0,1), e_ij ~ LNORM(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()






###################################################
#-----------E ~ UNIF(1,2), U ~ N(0,1)------------#
###################################################
#######################################
#                                     #
# error_i ~ UNIF(1,2), u_i ~ N(0,1)  #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
eunif1 <- lqmm_int_eunif(B=1000, N=30, ni=9, beta_0=2, sdi=1, min=1, max=2, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
eunif2 <- lqmm_int_eunif(B=1000, N=30, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
eunif3 <- lqmm_int_eunif(B=1000, N=100, ni=9, beta_0=2, sdi=1, min=1, max=2, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
eunif4 <- lqmm_int_eunif(B=1000, N=100, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
eunif5 <- lqmm_int_eunif(B=1000, N=500, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.05)

#######################################
#                                     #
# error_i ~ UNIF(1,2), u_i ~ NORM(0,1)  #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
eunif6 <- lqmm_int_eunif(B=1000, N=30, ni=9, beta_0=2, sdi=1, min=1, max=2, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
eunif7 <- lqmm_int_eunif(B=1000, N=30, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
eunif8 <- lqmm_int_eunif(B=1000, N=100, ni=9, beta_0=2, sdi=1, min=1, max=2, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
eunif9 <- lqmm_int_eunif(B=1000, N=100, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
eunif10 <- lqmm_int_eunif(B=1000, N=500, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.25)

#######################################
#                                     #
# error_i ~ UNIF(1,2), u_i ~ NORM(0,1)  #
#             tau = 0.5               #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
eunif11 <- lqmm_int_eunif(B=1000, N=30, ni=9, beta_0=2, sdi=1, min=1, max=2, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
eunif12 <- lqmm_int_eunif(B=1000, N=30, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
eunif13 <- lqmm_int_eunif(B=1000, N=100, ni=9, beta_0=2, sdi=1, min=1, max=2, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
eunif14 <- lqmm_int_eunif(B=1000, N=100, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
eunif15 <- lqmm_int_eunif(B=1000, N=500, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.5)

#######################################
#                                     #
#  error_i ~ UNIF(1,2), u_i ~ NORM(0,1) #
#             tau = 0.975             #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
eunif16 <- lqmm_int_eunif(B=1000, N=30, ni=9, beta_0=2, sdi=1, min=1, max=2, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
eunif17 <- lqmm_int_eunif(B=1000, N=30, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
eunif18 <- lqmm_int_eunif(B=1000, N=100, ni=9, beta_0=2, sdi=1, min=1, max=2, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
eunif19 <- lqmm_int_eunif(B=1000, N=100, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
eunif20 <- lqmm_int_eunif(B=1000, N=500, ni=18, beta_0=2, sdi=1, min=1, max=2, tau=0.975)

### FOR N=30 ###
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
bias_ui_eunif_N30 <- data.frame(eunif1$results_ranef$bias_i, eunif2$results_ranef$bias_i,
                                eunif6$results_ranef$bias_i, eunif7$results_ranef$bias_i,
                                eunif11$results_ranef$bias_i, eunif12$results_ranef$bias_i,
                                eunif16$results_ranef$bias_i, eunif17$results_ranef$bias_i)
colnames(bias_ui_eunif_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                 "ni9_t0.25", "ni18_t0.25",
                                 "ni9_t0.5", "ni18_t0.5",
                                 "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_eunif_N30, "bias_ui_eunif_N30.csv")
head(bias_ui_eunif_N30)

bias_ui_eunif_N30_long <- reshape(bias_ui_eunif_N30,
                                  direction = "long",
                                  varying = list(names(bias_ui_eunif_N30)[1:8]),
                                  v.names = "bias",
                                  timevar = "Setting",
                                  times = c("ni9_t0.05", "ni18_t0.05",
                                            "ni9_t0.25", "ni18_t0.25",
                                            "ni9_t0.5", "ni18_t0.5",
                                            "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_eunif_N100 <- data.frame(eunif3$results_ranef$bias_i, eunif4$results_ranef$bias_i,
                                 eunif8$results_ranef$bias_i, eunif9$results_ranef$bias_i,
                                 eunif13$results_ranef$bias_i, eunif14$results_ranef$bias_i,
                                 eunif18$results_ranef$bias_i, eunif19$results_ranef$bias_i)
colnames(bias_ui_eunif_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                  "ni9_t0.25", "ni18_t0.25",
                                  "ni9_t0.5", "ni18_t0.5",
                                  "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_eunif_N100, "bias_ui_eunif_N100.csv")
head(bias_ui_eunif_N100)

bias_ui_eunif_N100_long <- reshape(bias_ui_eunif_N100,
                                   direction = "long",
                                   varying = list(names(bias_ui_eunif_N100)[1:8]),
                                   v.names = "bias",
                                   timevar = "Setting",
                                   times = c("ni9_t0.05", "ni18_t0.05",
                                             "ni9_t0.25", "ni18_t0.25",
                                             "ni9_t0.5", "ni18_t0.5",
                                             "ni9_t0.975", "ni18_t0.975"))

mean_bias_eunif <- rbind(
  cbind(eunif1$results_ranef$mean_bias,
        eunif1$bias_estimates$mean_bias_beta0,
        eunif1$bias_estimates$mean_bias_b0),
  cbind(eunif2$results_ranef$mean_bias,
        eunif2$bias_estimates$mean_bias_beta0,
        eunif2$bias_estimates$mean_bias_b0),
  cbind(eunif3$results_ranef$mean_bias,
        eunif3$bias_estimates$mean_bias_beta0,
        eunif3$bias_estimates$mean_bias_b0),
  cbind(eunif4$results_ranef$mean_bias,
        eunif4$bias_estimates$mean_bias_beta0,
        eunif4$bias_estimates$mean_bias_b0),
  #  cbind(eunif5$results_ranef$mean_bias,
  #        eunif5$bias_estimates$mean_bias_beta0,
  #        eunif5$bias_estimates$mean_bias_b0),
  cbind(eunif6$results_ranef$mean_bias,
        eunif6$bias_estimates$mean_bias_beta0,
        eunif6$bias_estimates$mean_bias_b0),
  cbind(eunif7$results_ranef$mean_bias,
        eunif7$bias_estimates$mean_bias_beta0,
        eunif7$bias_estimates$mean_bias_b0),
  cbind(eunif8$results_ranef$mean_bias,
        eunif8$bias_estimates$mean_bias_beta0,
        eunif8$bias_estimates$mean_bias_b0),
  cbind(eunif9$results_ranef$mean_bias,
        eunif9$bias_estimates$mean_bias_beta0,
        eunif9$bias_estimates$mean_bias_b0),
  #  cbind(eunif10$results_ranef$mean_bias,
  #        eunif10$bias_estimates$mean_bias_beta0,
  #        eunif10$bias_estimates$mean_bias_b0),
  cbind(eunif11$results_ranef$mean_bias,
        eunif11$bias_estimates$mean_bias_beta0,
        eunif11$bias_estimates$mean_bias_b0),
  cbind(eunif12$results_ranef$mean_bias,
        eunif12$bias_estimates$mean_bias_beta0,
        eunif12$bias_estimates$mean_bias_b0),
  cbind(eunif13$results_ranef$mean_bias,
        eunif13$bias_estimates$mean_bias_beta0,
        eunif13$bias_estimates$mean_bias_b0),
  cbind(eunif14$results_ranef$mean_bias,
        eunif14$bias_estimates$mean_bias_beta0,
        eunif14$bias_estimates$mean_bias_b0),
  #  cbind(eunif15$results_ranef$mean_bias,
  #        eunif15$bias_estimates$mean_bias_beta0,
  #        eunif15$bias_estimates$mean_bias_b0),
  cbind(eunif16$results_ranef$mean_bias,
        eunif16$bias_estimates$mean_bias_beta0,
        eunif16$bias_estimates$mean_bias_b0),
  cbind(eunif17$results_ranef$mean_bias,
        eunif17$bias_estimates$mean_bias_beta0,
        eunif17$bias_estimates$mean_bias_b0),
  cbind(eunif18$results_ranef$mean_bias,
        eunif18$bias_estimates$mean_bias_beta0,
        eunif18$bias_estimates$mean_bias_b0),
  cbind(eunif19$results_ranef$mean_bias,
        eunif19$bias_estimates$mean_bias_beta0,
        eunif19$bias_estimates$mean_bias_b0))
#  cbind(eunif20$results_ranef$mean_bias,
#        eunif20$bias_estimates$mean_bias_beta0,
#        eunif20$bias_estimates$mean_bias_b0))


rownames(mean_bias_eunif) <- c("eunif1", "eunif2", "eunif3", "eunif4", 
                               "eunif6", "eunif7", "eunif8", "eunif9", 
                               "eunif11", "eunif12", "eunif13", "eunif14", 
                               "eunif16", "eunif17", "eunif18", "eunif19")

colnames(mean_bias_eunif) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_eunif, "mean_bias_eunif.csv")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")
png("Bias_ui_N30_eunif.png", width = 600, height = 500)
ggplot(bias_ui_eunif_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ NORM(0,1), e_ij ~ UNIFORM(1,2)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_eunif.png", width = 600, height = 500)
ggplot(bias_ui_eunif_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ NORM(0,1), e_ij ~ UNIFORM(1,2)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()





###################################################
#-----------E ~N(0,1) , U ~ UNIF(1,2)------------#
###################################################
#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ UNIF(1,2)  #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
uunif1 <- lqmm_int_uunif(B=1000, N=30, ni=9, beta_0=2, sd=1, min=1, max=2, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
uunif2 <- lqmm_int_uunif(B=1000, N=30, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
uunif3 <- lqmm_int_uunif(B=1000, N=100, ni=9, beta_0=2, sd=1, min=1, max=2, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
uunif4 <- lqmm_int_uunif(B=1000, N=100, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
uunif5 <- lqmm_int_uunif(B=1000, N=500, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.05)

#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ UNIF(1,2)  #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
uunif6 <- lqmm_int_uunif(B=1000, N=30, ni=9, beta_0=2, sd=1, min=1, max=2, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
uunif7 <- lqmm_int_uunif(B=1000, N=30, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
uunif8 <- lqmm_int_uunif(B=1000, N=100, ni=9, beta_0=2, sd=1, min=1, max=2, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
uunif9 <- lqmm_int_uunif(B=1000, N=100, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
uunif10 <- lqmm_int_uunif(B=1000, N=500, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.25)

#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ UNIF(1,2)  #
#             tau = 0.5               #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
uunif11 <- lqmm_int_uunif(B=1000, N=30, ni=9, beta_0=2, sd=1, min=1, max=2, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
uunif12 <- lqmm_int_uunif(B=1000, N=30, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
uunif13 <- lqmm_int_uunif(B=1000, N=100, ni=9, beta_0=2, sd=1, min=1, max=2, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
uunif14 <- lqmm_int_uunif(B=1000, N=100, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
uunif15 <- lqmm_int_uunif(B=1000, N=500, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.5)

#######################################
#                                     #
#  error_i ~ N(0,1), u_i ~ UNIF(1,2) #
#             tau = 0.975             #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
uunif16 <- lqmm_int_uunif(B=1000, N=30, ni=9, beta_0=2, sd=1, min=1, max=2, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
uunif17 <- lqmm_int_uunif(B=1000, N=30, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
uunif18 <- lqmm_int_uunif(B=1000, N=100, ni=9, beta_0=2, sd=1, min=1, max=2, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
uunif19 <- lqmm_int_uunif(B=1000, N=100, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
uunif20 <- lqmm_int_uunif(B=1000, N=500, ni=18, beta_0=2, sd=1, min=1, max=2, tau=0.975)


### FOR N=30 ###
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
bias_ui_uunif_N30 <- data.frame(uunif1$results_ranef$bias_i, uunif2$results_ranef$bias_i,
                                uunif6$results_ranef$bias_i, uunif7$results_ranef$bias_i,
                                uunif11$results_ranef$bias_i, uunif12$results_ranef$bias_i,
                                uunif16$results_ranef$bias_i, uunif17$results_ranef$bias_i)
colnames(bias_ui_uunif_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                 "ni9_t0.25", "ni18_t0.25",
                                 "ni9_t0.5", "ni18_t0.5",
                                 "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_uunif_N30, "bias_ui_uunif_N30.csv")
head(bias_ui_uunif_N30)

bias_ui_uunif_N30_long <- reshape(bias_ui_uunif_N30,
                                  direction = "long",
                                  varying = list(names(bias_ui_uunif_N30)[1:8]),
                                  v.names = "bias",
                                  timevar = "Setting",
                                  times = c("ni9_t0.05", "ni18_t0.05",
                                            "ni9_t0.25", "ni18_t0.25",
                                            "ni9_t0.5", "ni18_t0.5",
                                            "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_uunif_N100 <- data.frame(uunif3$results_ranef$bias_i, uunif4$results_ranef$bias_i,
                                 uunif8$results_ranef$bias_i, uunif9$results_ranef$bias_i,
                                 uunif13$results_ranef$bias_i, uunif14$results_ranef$bias_i,
                                 uunif18$results_ranef$bias_i, uunif19$results_ranef$bias_i)
colnames(bias_ui_uunif_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                  "ni9_t0.25", "ni18_t0.25",
                                  "ni9_t0.5", "ni18_t0.5",
                                  "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_uunif_N100, "bias_ui_uunif_N100.csv")
head(bias_ui_uunif_N100)

bias_ui_uunif_N100_long <- reshape(bias_ui_uunif_N100,
                                   direction = "long",
                                   varying = list(names(bias_ui_uunif_N100)[1:8]),
                                   v.names = "bias",
                                   timevar = "Setting",
                                   times = c("ni9_t0.05", "ni18_t0.05",
                                             "ni9_t0.25", "ni18_t0.25",
                                             "ni9_t0.5", "ni18_t0.5",
                                             "ni9_t0.975", "ni18_t0.975"))

mean_bias_uunif <- rbind(
  cbind(uunif1$results_ranef$mean_bias,
        uunif1$bias_estimates$mean_bias_beta0,
        uunif1$bias_estimates$mean_bias_b0),
  cbind(uunif2$results_ranef$mean_bias,
        uunif2$bias_estimates$mean_bias_beta0,
        uunif2$bias_estimates$mean_bias_b0),
  cbind(uunif3$results_ranef$mean_bias,
        uunif3$bias_estimates$mean_bias_beta0,
        uunif3$bias_estimates$mean_bias_b0),
  cbind(uunif4$results_ranef$mean_bias,
        uunif4$bias_estimates$mean_bias_beta0,
        uunif4$bias_estimates$mean_bias_b0),
  #  cbind(uunif5$results_ranef$mean_bias,
  #        uunif5$bias_estimates$mean_bias_beta0,
  #        uunif5$bias_estimates$mean_bias_b0),
  cbind(uunif6$results_ranef$mean_bias,
        uunif6$bias_estimates$mean_bias_beta0,
        uunif6$bias_estimates$mean_bias_b0),
  cbind(uunif7$results_ranef$mean_bias,
        uunif7$bias_estimates$mean_bias_beta0,
        uunif7$bias_estimates$mean_bias_b0),
  cbind(uunif8$results_ranef$mean_bias,
        uunif8$bias_estimates$mean_bias_beta0,
        uunif8$bias_estimates$mean_bias_b0),
  cbind(uunif9$results_ranef$mean_bias,
        uunif9$bias_estimates$mean_bias_beta0,
        uunif9$bias_estimates$mean_bias_b0),
  #  cbind(uunif10$results_ranef$mean_bias,
  #        uunif10$bias_estimates$mean_bias_beta0,
  #        uunif10$bias_estimates$mean_bias_b0),
  cbind(uunif11$results_ranef$mean_bias,
        uunif11$bias_estimates$mean_bias_beta0,
        uunif11$bias_estimates$mean_bias_b0),
  cbind(uunif12$results_ranef$mean_bias,
        uunif12$bias_estimates$mean_bias_beta0,
        uunif12$bias_estimates$mean_bias_b0),
  cbind(uunif13$results_ranef$mean_bias,
        uunif13$bias_estimates$mean_bias_beta0,
        uunif13$bias_estimates$mean_bias_b0),
  cbind(uunif14$results_ranef$mean_bias,
        uunif14$bias_estimates$mean_bias_beta0,
        uunif14$bias_estimates$mean_bias_b0),
  #  cbind(uunif15$results_ranef$mean_bias,
  #        uunif15$bias_estimates$mean_bias_beta0,
  #        uunif15$bias_estimates$mean_bias_b0),
  cbind(uunif16$results_ranef$mean_bias,
        uunif16$bias_estimates$mean_bias_beta0,
        uunif16$bias_estimates$mean_bias_b0),
  cbind(uunif17$results_ranef$mean_bias,
        uunif17$bias_estimates$mean_bias_beta0,
        uunif17$bias_estimates$mean_bias_b0),
  cbind(uunif18$results_ranef$mean_bias,
        uunif18$bias_estimates$mean_bias_beta0,
        uunif18$bias_estimates$mean_bias_b0),
  cbind(uunif19$results_ranef$mean_bias,
        uunif19$bias_estimates$mean_bias_beta0,
        uunif19$bias_estimates$mean_bias_b0))
#  cbind(uunif20$results_ranef$mean_bias,
#        uunif20$bias_estimates$mean_bias_beta0,
#        uunif20$bias_estimates$mean_bias_b0))


rownames(mean_bias_uunif) <- c("uunif1", "uunif2", "uunif3", "uunif4", 
                               "uunif6", "uunif7", "uunif8", "uunif9", 
                               "uunif11", "uunif12", "uunif13", "uunif14", 
                               "uunif16", "uunif17", "uunif18", "uunif19")

colnames(mean_bias_uunif) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_uunif, "mean_bias_uunif.csv")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")
png("Bias_ui_N30_uunif.png", width = 600, height = 500)
ggplot(bias_ui_uunif_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ UNIFORM(1,2), e_ij ~ N(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_uunif.png", width = 600, height = 500)
ggplot(bias_ui_uunif_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ UNIFORM(1,2), e_ij ~ NORM(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()






###################################################
#-----------E ~ CHISQ(4) , U ~ NORM(0,1)------------#
###################################################
#######################################
#                                     #
# error_i ~ CHISQ(4), u_i ~ N(0,1)  #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
echisq1 <- lqmm_int_echisq(B=1000, N=30, ni=9, beta_0=2, sdi=1, df=4, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
echisq2 <- lqmm_int_echisq(B=1000, N=30, ni=18, beta_0=2, sdi=1, df=4, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
echisq3 <- lqmm_int_echisq(B=1000, N=100, ni=9, beta_0=2, sdi=1, df=4, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
echisq4 <- lqmm_int_echisq(B=1000, N=100, ni=18, beta_0=2, sdi=1, df=4, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
echisq5 <- lqmm_int_echisq(B=1000, N=500, ni=18, beta_0=2, sdi=1, df=4, tau=0.05)

#######################################
#                                     #
# error_i ~ CHISQ(4), u_i ~ N(0,1)  #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
echisq6 <- lqmm_int_echisq(B=1000, N=30, ni=9, beta_0=2, sdi=1, df=4, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
echisq7 <- lqmm_int_echisq(B=1000, N=30, ni=18, beta_0=2, sdi=1, df=4, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
echisq8 <- lqmm_int_echisq(B=1000, N=100, ni=9, beta_0=2, sdi=1, df=4, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
echisq9 <- lqmm_int_echisq(B=1000, N=100, ni=18, beta_0=2, sdi=1, df=4, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
echisq10 <- lqmm_int_echisq(B=1000, N=500, ni=18, beta_0=2, sdi=1, df=4, tau=0.25)

#######################################
#                                     #
# error_i ~ CHISQ(4), u_i ~ N(0,1)  #
#             tau = 0.5               #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
echisq11 <- lqmm_int_echisq(B=1000, N=30, ni=9, beta_0=2, sdi=1, df=4, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
echisq12 <- lqmm_int_echisq(B=1000, N=30, ni=18, beta_0=2, sdi=1, df=4, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
echisq13 <- lqmm_int_echisq(B=1000, N=100, ni=9, beta_0=2, sdi=1, df=4, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
echisq14 <- lqmm_int_echisq(B=1000, N=100, ni=18, beta_0=2, sdi=1, df=4, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
echisq15 <- lqmm_int_echisq(B=1000, N=500, ni=18, beta_0=2, sdi=1, df=4, tau=0.5)

#######################################
#                                     #
#  error_i ~ CHISQ(4), u_i ~ N(0,1) #
#             tau = 0.975             #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
echisq16 <- lqmm_int_echisq(B=1000, N=30, ni=9, beta_0=2, sdi=1, df=4, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
echisq17 <- lqmm_int_echisq(B=1000, N=30, ni=18, beta_0=2, sdi=1, df=4, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
echisq18 <- lqmm_int_echisq(B=1000, N=100, ni=9, beta_0=2, sdi=1, df=4, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
echisq19 <- lqmm_int_echisq(B=1000, N=100, ni=18, beta_0=2, sdi=1, df=4, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
echisq20 <- lqmm_int_echisq(B=1000, N=500, ni=18, beta_0=2, sdi=1, df=4, tau=0.975)


### FOR N=30 ###
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
bias_ui_echisq_N30 <- data.frame(echisq1$results_ranef$bias_i, echisq2$results_ranef$bias_i,
                                 echisq6$results_ranef$bias_i, echisq7$results_ranef$bias_i,
                                 echisq11$results_ranef$bias_i, echisq12$results_ranef$bias_i,
                                 echisq16$results_ranef$bias_i, echisq17$results_ranef$bias_i)
colnames(bias_ui_echisq_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                  "ni9_t0.25", "ni18_t0.25",
                                  "ni9_t0.5", "ni18_t0.5",
                                  "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_echisq_N30, "bias_ui_echisq_N30.csv")
head(bias_ui_echisq_N30)

bias_ui_echisq_N30_long <- reshape(bias_ui_echisq_N30,
                                   direction = "long",
                                   varying = list(names(bias_ui_echisq_N30)[1:8]),
                                   v.names = "bias",
                                   timevar = "Setting",
                                   times = c("ni9_t0.05", "ni18_t0.05",
                                             "ni9_t0.25", "ni18_t0.25",
                                             "ni9_t0.5", "ni18_t0.5",
                                             "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_echisq_N100 <- data.frame(echisq3$results_ranef$bias_i, echisq4$results_ranef$bias_i,
                                  echisq8$results_ranef$bias_i, echisq9$results_ranef$bias_i,
                                  echisq13$results_ranef$bias_i, echisq14$results_ranef$bias_i,
                                  echisq18$results_ranef$bias_i, echisq19$results_ranef$bias_i)
colnames(bias_ui_echisq_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                   "ni9_t0.25", "ni18_t0.25",
                                   "ni9_t0.5", "ni18_t0.5",
                                   "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_echisq_N100, "bias_ui_echisq_N100.csv")
head(bias_ui_echisq_N100)

bias_ui_echisq_N100_long <- reshape(bias_ui_echisq_N100,
                                    direction = "long",
                                    varying = list(names(bias_ui_echisq_N100)[1:8]),
                                    v.names = "bias",
                                    timevar = "Setting",
                                    times = c("ni9_t0.05", "ni18_t0.05",
                                              "ni9_t0.25", "ni18_t0.25",
                                              "ni9_t0.5", "ni18_t0.5",
                                              "ni9_t0.975", "ni18_t0.975"))

mean_bias_echisq <- rbind(
  cbind(echisq1$results_ranef$mean_bias,
        echisq1$bias_estimates$mean_bias_beta0,
        echisq1$bias_estimates$mean_bias_b0),
  cbind(echisq2$results_ranef$mean_bias,
        echisq2$bias_estimates$mean_bias_beta0,
        echisq2$bias_estimates$mean_bias_b0),
  cbind(echisq3$results_ranef$mean_bias,
        echisq3$bias_estimates$mean_bias_beta0,
        echisq3$bias_estimates$mean_bias_b0),
  cbind(echisq4$results_ranef$mean_bias,
        echisq4$bias_estimates$mean_bias_beta0,
        echisq4$bias_estimates$mean_bias_b0),
  #  cbind(echisq5$results_ranef$mean_bias,
  #        echisq5$bias_estimates$mean_bias_beta0,
  #        echisq5$bias_estimates$mean_bias_b0),
  cbind(echisq6$results_ranef$mean_bias,
        echisq6$bias_estimates$mean_bias_beta0,
        echisq6$bias_estimates$mean_bias_b0),
  cbind(echisq7$results_ranef$mean_bias,
        echisq7$bias_estimates$mean_bias_beta0,
        echisq7$bias_estimates$mean_bias_b0),
  cbind(echisq8$results_ranef$mean_bias,
        echisq8$bias_estimates$mean_bias_beta0,
        echisq8$bias_estimates$mean_bias_b0),
  cbind(echisq9$results_ranef$mean_bias,
        echisq9$bias_estimates$mean_bias_beta0,
        echisq9$bias_estimates$mean_bias_b0),
  #  cbind(echisq10$results_ranef$mean_bias,
  #        echisq10$bias_estimates$mean_bias_beta0,
  #        echisq10$bias_estimates$mean_bias_b0),
  cbind(echisq11$results_ranef$mean_bias,
        echisq11$bias_estimates$mean_bias_beta0,
        echisq11$bias_estimates$mean_bias_b0),
  cbind(echisq12$results_ranef$mean_bias,
        echisq12$bias_estimates$mean_bias_beta0,
        echisq12$bias_estimates$mean_bias_b0),
  cbind(echisq13$results_ranef$mean_bias,
        echisq13$bias_estimates$mean_bias_beta0,
        echisq13$bias_estimates$mean_bias_b0),
  cbind(echisq14$results_ranef$mean_bias,
        echisq14$bias_estimates$mean_bias_beta0,
        echisq14$bias_estimates$mean_bias_b0),
  #  cbind(echisq15$results_ranef$mean_bias,
  #        echisq15$bias_estimates$mean_bias_beta0,
  #        echisq15$bias_estimates$mean_bias_b0),
  cbind(echisq16$results_ranef$mean_bias,
        echisq16$bias_estimates$mean_bias_beta0,
        echisq16$bias_estimates$mean_bias_b0),
  cbind(echisq17$results_ranef$mean_bias,
        echisq17$bias_estimates$mean_bias_beta0,
        echisq17$bias_estimates$mean_bias_b0),
  cbind(echisq18$results_ranef$mean_bias,
        echisq18$bias_estimates$mean_bias_beta0,
        echisq18$bias_estimates$mean_bias_b0),
  cbind(echisq19$results_ranef$mean_bias,
        echisq19$bias_estimates$mean_bias_beta0,
        echisq19$bias_estimates$mean_bias_b0))
#  cbind(echisq20$results_ranef$mean_bias,
#        echisq20$bias_estimates$mean_bias_beta0,
#        echisq20$bias_estimates$mean_bias_b0))


rownames(mean_bias_echisq) <- c("echisq1", "echisq2", "echisq3", "echisq4", 
                                "echisq6", "echisq7", "echisq8", "echisq9", 
                                "echisq11", "echisq12", "echisq13", "echisq14", 
                                "echisq16", "echisq17", "echisq18", "echisq19")

colnames(mean_bias_echisq) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_echisq, "mean_bias_echisq.csv")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")
png("Bias_ui_N30_echisq.png", width = 600, height = 500)
ggplot(bias_ui_echisq_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ N(0,1), e_ij ~ CHISQ(4)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_echisq.png", width = 600, height = 500)
ggplot(bias_ui_echisq_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ N(0,1), e_ij ~ CHISQ(4)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()




###################################################
#-----------E ~ NORM(0,1), U ~ CHISQ(4)------------#
###################################################
#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ CHISQ(4)  #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
uchisq1 <- lqmm_int_uchisq(B=1000, N=30, ni=9, beta_0=2, sd=1, df=4, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
uchisq2 <- lqmm_int_uchisq(B=1000, N=30, ni=18, beta_0=2, sd=1, df=4, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
uchisq3 <- lqmm_int_uchisq(B=1000, N=100, ni=9, beta_0=2, sd=1, df=4, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
uchisq4 <- lqmm_int_uchisq(B=1000, N=100, ni=18, beta_0=2, sd=1, df=4, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
uchisq5 <- lqmm_int_uchisq(B=1000, N=500, ni=18, beta_0=2, sd=1, df=4, tau=0.05)

#######################################
#                                     #
# error_i ~ N(0,1), u_i ~ CHISQ(4)  #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
uchisq6 <- lqmm_int_uchisq(B=1000, N=30, ni=9, beta_0=2, sd=1, df=4, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
uchisq7 <- lqmm_int_uchisq(B=1000, N=30, ni=18, beta_0=2, sd=1, df=4, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
uchisq8 <- lqmm_int_uchisq(B=1000, N=100, ni=9, beta_0=2, sd=1, df=4, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
uchisq9 <- lqmm_int_uchisq(B=1000, N=100, ni=18, beta_0=2, sd=1, df=4, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
uchisq10 <- lqmm_int_uchisq(B=1000, N=500, ni=18, beta_0=2, sd=1, df=4, tau=0.25)

#######################################
#                                     #
# error_i ~ NORM(0,1), u_i ~ CHISQ(4)  #
#             tau = 0.5               #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
uchisq11 <- lqmm_int_uchisq(B=1000, N=30, ni=9, beta_0=2, sd=1, df=4, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
uchisq12 <- lqmm_int_uchisq(B=1000, N=30, ni=18, beta_0=2, sd=1, df=4, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
uchisq13 <- lqmm_int_uchisq(B=1000, N=100, ni=9, beta_0=2, sd=1, df=4, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
uchisq14 <- lqmm_int_uchisq(B=1000, N=100, ni=18, beta_0=2, sd=1, df=4, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
uchisq15 <- lqmm_int_uchisq(B=1000, N=500, ni=18, beta_0=2, sd=1, df=4, tau=0.5)

#######################################
#                                     #
#  error_i ~ NORM(0,1), u_i ~ CHISQ(4) #
#             tau = 0.975             #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
uchisq16 <- lqmm_int_uchisq(B=1000, N=30, ni=9, beta_0=2, sd=1, df=4, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
uchisq17 <- lqmm_int_uchisq(B=1000, N=30, ni=18, beta_0=2, sd=1, df=4, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
uchisq18 <- lqmm_int_uchisq(B=1000, N=100, ni=9, beta_0=2, sd=1, df=4, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
uchisq19 <- lqmm_int_uchisq(B=1000, N=100, ni=18, beta_0=2, sd=1, df=4, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
uchisq20 <- lqmm_int_uchisq(B=1000, N=500, ni=18, beta_0=2, sd=1, df=4, tau=0.975)

### FOR N=30 ###
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
bias_ui_uchisq_N30 <- data.frame(uchisq1$results_ranef$bias_i, uchisq2$results_ranef$bias_i,
                                 uchisq6$results_ranef$bias_i, uchisq7$results_ranef$bias_i,
                                 uchisq11$results_ranef$bias_i, uchisq12$results_ranef$bias_i,
                                 uchisq16$results_ranef$bias_i, uchisq17$results_ranef$bias_i)
colnames(bias_ui_uchisq_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                  "ni9_t0.25", "ni18_t0.25",
                                  "ni9_t0.5", "ni18_t0.5",
                                  "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_uchisq_N30, "bias_ui_uchisq_N30.csv")
head(bias_ui_uchisq_N30)

bias_ui_uchisq_N30_long <- reshape(bias_ui_uchisq_N30,
                                   direction = "long",
                                   varying = list(names(bias_ui_uchisq_N30)[1:8]),
                                   v.names = "bias",
                                   timevar = "Setting",
                                   times = c("ni9_t0.05", "ni18_t0.05",
                                             "ni9_t0.25", "ni18_t0.25",
                                             "ni9_t0.5", "ni18_t0.5",
                                             "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_uchisq_N100 <- data.frame(uchisq3$results_ranef$bias_i, uchisq4$results_ranef$bias_i,
                                  uchisq8$results_ranef$bias_i, uchisq9$results_ranef$bias_i,
                                  uchisq13$results_ranef$bias_i, uchisq14$results_ranef$bias_i,
                                  uchisq18$results_ranef$bias_i, uchisq19$results_ranef$bias_i)
colnames(bias_ui_uchisq_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                   "ni9_t0.25", "ni18_t0.25",
                                   "ni9_t0.5", "ni18_t0.5",
                                   "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_uchisq_N100, "bias_ui_uchisq_N100.csv")
head(bias_ui_uchisq_N100)

bias_ui_uchisq_N100_long <- reshape(bias_ui_uchisq_N100,
                                    direction = "long",
                                    varying = list(names(bias_ui_uchisq_N100)[1:8]),
                                    v.names = "bias",
                                    timevar = "Setting",
                                    times = c("ni9_t0.05", "ni18_t0.05",
                                              "ni9_t0.25", "ni18_t0.25",
                                              "ni9_t0.5", "ni18_t0.5",
                                              "ni9_t0.975", "ni18_t0.975"))

mean_bias_uchisq <- rbind(
  cbind(uchisq1$results_ranef$mean_bias,
        uchisq1$bias_estimates$mean_bias_beta0,
        uchisq1$bias_estimates$mean_bias_b0),
  cbind(uchisq2$results_ranef$mean_bias,
        uchisq2$bias_estimates$mean_bias_beta0,
        uchisq2$bias_estimates$mean_bias_b0),
  cbind(uchisq3$results_ranef$mean_bias,
        uchisq3$bias_estimates$mean_bias_beta0,
        uchisq3$bias_estimates$mean_bias_b0),
  cbind(uchisq4$results_ranef$mean_bias,
        uchisq4$bias_estimates$mean_bias_beta0,
        uchisq4$bias_estimates$mean_bias_b0),
  #  cbind(uchisq5$results_ranef$mean_bias,
  #        uchisq5$bias_estimates$mean_bias_beta0,
  #        uchisq5$bias_estimates$mean_bias_b0),
  cbind(uchisq6$results_ranef$mean_bias,
        uchisq6$bias_estimates$mean_bias_beta0,
        uchisq6$bias_estimates$mean_bias_b0),
  cbind(uchisq7$results_ranef$mean_bias,
        uchisq7$bias_estimates$mean_bias_beta0,
        uchisq7$bias_estimates$mean_bias_b0),
  cbind(uchisq8$results_ranef$mean_bias,
        uchisq8$bias_estimates$mean_bias_beta0,
        uchisq8$bias_estimates$mean_bias_b0),
  cbind(uchisq9$results_ranef$mean_bias,
        uchisq9$bias_estimates$mean_bias_beta0,
        uchisq9$bias_estimates$mean_bias_b0),
  #  cbind(uchisq10$results_ranef$mean_bias,
  #        uchisq10$bias_estimates$mean_bias_beta0,
  #        uchisq10$bias_estimates$mean_bias_b0),
  cbind(uchisq11$results_ranef$mean_bias,
        uchisq11$bias_estimates$mean_bias_beta0,
        uchisq11$bias_estimates$mean_bias_b0),
  cbind(uchisq12$results_ranef$mean_bias,
        uchisq12$bias_estimates$mean_bias_beta0,
        uchisq12$bias_estimates$mean_bias_b0),
  cbind(uchisq13$results_ranef$mean_bias,
        uchisq13$bias_estimates$mean_bias_beta0,
        uchisq13$bias_estimates$mean_bias_b0),
  cbind(uchisq14$results_ranef$mean_bias,
        uchisq14$bias_estimates$mean_bias_beta0,
        uchisq14$bias_estimates$mean_bias_b0),
  #  cbind(uchisq15$results_ranef$mean_bias,
  #        uchisq15$bias_estimates$mean_bias_beta0,
  #        uchisq15$bias_estimates$mean_bias_b0),
  cbind(uchisq16$results_ranef$mean_bias,
        uchisq16$bias_estimates$mean_bias_beta0,
        uchisq16$bias_estimates$mean_bias_b0),
  cbind(uchisq17$results_ranef$mean_bias,
        uchisq17$bias_estimates$mean_bias_beta0,
        uchisq17$bias_estimates$mean_bias_b0),
  cbind(uchisq18$results_ranef$mean_bias,
        uchisq18$bias_estimates$mean_bias_beta0,
        uchisq18$bias_estimates$mean_bias_b0),
  cbind(uchisq19$results_ranef$mean_bias,
        uchisq19$bias_estimates$mean_bias_beta0,
        uchisq19$bias_estimates$mean_bias_b0))
#  cbind(uchisq20$results_ranef$mean_bias,
#        uchisq20$bias_estimates$mean_bias_beta0,
#        uchisq20$bias_estimates$mean_bias_b0))


rownames(mean_bias_uchisq) <- c("uchisq1", "uchisq2", "uchisq3", "uchisq4", 
                                "uchisq6", "uchisq7", "uchisq8", "uchisq9", 
                                "uchisq11", "uchisq12", "uchisq13", "uchisq14", 
                                "uchisq16", "uchisq17", "uchisq18", "uchisq19")

colnames(mean_bias_uchisq) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_uchisq, "mean_bias_uchisq.csv")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")
png("Bias_ui_N30_uchisq.png", width = 600, height = 500)
ggplot(bias_ui_uchisq_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ CHISQ(4), e_ij ~ N(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_uchisq.png", width = 600, height = 500)
ggplot(bias_ui_uchisq_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ CHISQ(4), e_ij ~ N(0,1)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()




###################################################
#-----------E ~ CHISQ(4), U ~ CHISQ(4)------------#
###################################################
#######################################
#                                     #
# error_i ~ CHISQ(4), u_i ~ CHISQ(4)  #
#             tau = 0.05              #
#######################################

#Setting 1: N=30, ni=9, B=1000, beta_0=2
euchisq1 <- lqmm_int_euchisq(B=1000, N=30, ni=9, beta_0=2, df=4, dfi=4, tau=0.05)
#Setting 2: N=30, ni=18, B=1000, beta_0=2
euchisq2 <- lqmm_int_euchisq(B=1000, N=30, ni=18, beta_0=2, df=4, dfi=4, tau=0.05)
#Setting 3: N=100, ni=9, B=1000, beta_0=2
euchisq3 <- lqmm_int_euchisq(B=1000, N=100, ni=9, beta_0=2, df=4, dfi=4, tau=0.05)
#Setting 4: N=100, ni=18, B=1000, beta_0=2
euchisq4 <- lqmm_int_euchisq(B=1000, N=100, ni=18, beta_0=2, df=4, dfi=4, tau=0.05)
#Setting 5: N=500, ni=18, B=1000, beta_0=2
euchisq5 <- lqmm_int_euchisq(B=1000, N=500, ni=18, beta_0=2, df=4, dfi=4, tau=0.05)

#######################################
#                                     #
# error_i ~ CHISQ(4), u_i ~ CHISQ(4)  #
#             tau = 0.25              #
#######################################

#Setting 6: N=30, ni=9, B=1000, beta_0=2
euchisq6 <- lqmm_int_euchisq(B=1000, N=30, ni=9, beta_0=2, df=4, dfi=4, tau=0.25)
#Setting 7: N=30, ni=18, B=1000, beta_0=2
euchisq7 <- lqmm_int_euchisq(B=1000, N=30, ni=18, beta_0=2, df=4, dfi=4, tau=0.25)
#Setting 8: N=100, ni=9, B=1000, beta_0=2
euchisq8 <- lqmm_int_euchisq(B=1000, N=100, ni=9, beta_0=2, df=4, dfi=4, tau=0.25)
#Setting 9: N=100, ni=18, B=1000, beta_0=2
euchisq9 <- lqmm_int_euchisq(B=1000, N=100, ni=18, beta_0=2, df=4, dfi=4, tau=0.25)
#Setting 10: N=500, ni=18, B=1000, beta_0=2
euchisq10 <- lqmm_int_euchisq(B=1000, N=500, ni=18, beta_0=2, df=4, dfi=4, tau=0.25)

#######################################
#                                     #
# error_i ~ CHISQ(4), u_i ~ CHISQ(4)  #
#             tau = 0.5               #
#######################################

#Setting 11: N=30, ni=9, B=1000, beta_0=2
euchisq11 <- lqmm_int_euchisq(B=1000, N=30, ni=9, beta_0=2, df=4, dfi=4, tau=0.5)
#Setting 12: N=30, ni=18, B=1000, beta_0=2
euchisq12 <- lqmm_int_euchisq(B=1000, N=30, ni=18, beta_0=2, df=4, dfi=4, tau=0.5)
#Setting 13: N=100, ni=9, B=1000, beta_0=2
euchisq13 <- lqmm_int_euchisq(B=1000, N=100, ni=9, beta_0=2, df=4, dfi=4, tau=0.5)
#Setting 14: N=100, ni=18, B=1000, beta_0=2
euchisq14 <- lqmm_int_euchisq(B=1000, N=100, ni=18, beta_0=2, df=4, dfi=4, tau=0.5)
#Setting 15: N=500, ni=18, B=1000, beta_0=2
euchisq15 <- lqmm_int_euchisq(B=1000, N=500, ni=18, beta_0=2, df=4, dfi=4, tau=0.5)

#######################################
#                                     #
#  error_i ~ CHISQ(4), u_i ~ CHISQ(4) #
#             tau = 0.975             #
#######################################

#Setting 16: N=30, ni=9, B=1000, beta_0=2
euchisq16 <- lqmm_int_euchisq(B=1000, N=30, ni=9, beta_0=2, df=4, dfi=4, tau=0.975)
#Setting 17: N=30, ni=18, B=1000, beta_0=2
euchisq17 <- lqmm_int_euchisq(B=1000, N=30, ni=18, beta_0=2, df=4, dfi=4, tau=0.975)
#Setting 18: N=100, ni=9, B=1000, beta_0=2
euchisq18 <- lqmm_int_euchisq(B=1000, N=100, ni=9, beta_0=2, df=4, dfi=4, tau=0.975)
#Setting 19: N=100, ni=18, B=1000, beta_0=2
euchisq19 <- lqmm_int_euchisq(B=1000, N=100, ni=18, beta_0=2, df=4, dfi=4, tau=0.975)
#Setting 20: N=500, ni=18, B=1000, beta_0=2
euchisq20 <- lqmm_int_euchisq(B=1000, N=500, ni=18, beta_0=2, df=4, dfi=4, tau=0.975)

### FOR N=30 ###
setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation")
bias_ui_euchisq_N30 <- data.frame(euchisq1$results_ranef$bias_i, euchisq2$results_ranef$bias_i,
                                  euchisq6$results_ranef$bias_i, euchisq7$results_ranef$bias_i,
                                  euchisq11$results_ranef$bias_i, euchisq12$results_ranef$bias_i,
                                  euchisq16$results_ranef$bias_i, euchisq17$results_ranef$bias_i)
colnames(bias_ui_euchisq_N30) <- c("ni9_t0.05", "ni18_t0.05",
                                   "ni9_t0.25", "ni18_t0.25",
                                   "ni9_t0.5", "ni18_t0.5",
                                   "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_euchisq_N30, "bias_ui_euchisq_N30.csv")
head(bias_ui_euchisq_N30)

bias_ui_euchisq_N30_long <- reshape(bias_ui_euchisq_N30,
                                    direction = "long",
                                    varying = list(names(bias_ui_euchisq_N30)[1:8]),
                                    v.names = "bias",
                                    timevar = "Setting",
                                    times = c("ni9_t0.05", "ni18_t0.05",
                                              "ni9_t0.25", "ni18_t0.25",
                                              "ni9_t0.5", "ni18_t0.5",
                                              "ni9_t0.975", "ni18_t0.975"))


### FOR N=100 ###
bias_ui_euchisq_N100 <- data.frame(euchisq3$results_ranef$bias_i, euchisq4$results_ranef$bias_i,
                                   euchisq8$results_ranef$bias_i, euchisq9$results_ranef$bias_i,
                                   euchisq13$results_ranef$bias_i, euchisq14$results_ranef$bias_i,
                                   euchisq18$results_ranef$bias_i, euchisq19$results_ranef$bias_i)
colnames(bias_ui_euchisq_N100) <- c("ni9_t0.05", "ni18_t0.05",
                                    "ni9_t0.25", "ni18_t0.25",
                                    "ni9_t0.5", "ni18_t0.5",
                                    "ni9_t0.975", "ni18_t0.975")
write.csv(bias_ui_euchisq_N100, "bias_ui_euchisq_N100.csv")
head(bias_ui_euchisq_N100)

bias_ui_euchisq_N100_long <- reshape(bias_ui_euchisq_N100,
                                     direction = "long",
                                     varying = list(names(bias_ui_euchisq_N100)[1:8]),
                                     v.names = "bias",
                                     timevar = "Setting",
                                     times = c("ni9_t0.05", "ni18_t0.05",
                                               "ni9_t0.25", "ni18_t0.25",
                                               "ni9_t0.5", "ni18_t0.5",
                                               "ni9_t0.975", "ni18_t0.975"))

mean_bias_euchisq <- rbind(
  cbind(euchisq1$results_ranef$mean_bias,
        euchisq1$bias_estimates$mean_bias_beta0,
        euchisq1$bias_estimates$mean_bias_b0),
  cbind(euchisq2$results_ranef$mean_bias,
        euchisq2$bias_estimates$mean_bias_beta0,
        euchisq2$bias_estimates$mean_bias_b0),
  cbind(euchisq3$results_ranef$mean_bias,
        euchisq3$bias_estimates$mean_bias_beta0,
        euchisq3$bias_estimates$mean_bias_b0),
  cbind(euchisq4$results_ranef$mean_bias,
        euchisq4$bias_estimates$mean_bias_beta0,
        euchisq4$bias_estimates$mean_bias_b0),
  #  cbind(euchisq5$results_ranef$mean_bias,
  #        euchisq5$bias_estimates$mean_bias_beta0,
  #        euchisq5$bias_estimates$mean_bias_b0),
  cbind(euchisq6$results_ranef$mean_bias,
        euchisq6$bias_estimates$mean_bias_beta0,
        euchisq6$bias_estimates$mean_bias_b0),
  cbind(euchisq7$results_ranef$mean_bias,
        euchisq7$bias_estimates$mean_bias_beta0,
        euchisq7$bias_estimates$mean_bias_b0),
  cbind(euchisq8$results_ranef$mean_bias,
        euchisq8$bias_estimates$mean_bias_beta0,
        euchisq8$bias_estimates$mean_bias_b0),
  cbind(euchisq9$results_ranef$mean_bias,
        euchisq9$bias_estimates$mean_bias_beta0,
        euchisq9$bias_estimates$mean_bias_b0),
  #  cbind(euchisq10$results_ranef$mean_bias,
  #        euchisq10$bias_estimates$mean_bias_beta0,
  #        euchisq10$bias_estimates$mean_bias_b0),
  cbind(euchisq11$results_ranef$mean_bias,
        euchisq11$bias_estimates$mean_bias_beta0,
        euchisq11$bias_estimates$mean_bias_b0),
  cbind(euchisq12$results_ranef$mean_bias,
        euchisq12$bias_estimates$mean_bias_beta0,
        euchisq12$bias_estimates$mean_bias_b0),
  cbind(euchisq13$results_ranef$mean_bias,
        euchisq13$bias_estimates$mean_bias_beta0,
        euchisq13$bias_estimates$mean_bias_b0),
  cbind(euchisq14$results_ranef$mean_bias,
        euchisq14$bias_estimates$mean_bias_beta0,
        euchisq14$bias_estimates$mean_bias_b0),
  #  cbind(euchisq15$results_ranef$mean_bias,
  #        euchisq15$bias_estimates$mean_bias_beta0,
  #        euchisq15$bias_estimates$mean_bias_b0),
  cbind(euchisq16$results_ranef$mean_bias,
        euchisq16$bias_estimates$mean_bias_beta0,
        euchisq16$bias_estimates$mean_bias_b0),
  cbind(euchisq17$results_ranef$mean_bias,
        euchisq17$bias_estimates$mean_bias_beta0,
        euchisq17$bias_estimates$mean_bias_b0),
  cbind(euchisq18$results_ranef$mean_bias,
        euchisq18$bias_estimates$mean_bias_beta0,
        euchisq18$bias_estimates$mean_bias_b0),
  cbind(euchisq19$results_ranef$mean_bias,
        euchisq19$bias_estimates$mean_bias_beta0,
        euchisq19$bias_estimates$mean_bias_b0))
#  cbind(euchisq20$results_ranef$mean_bias,
#        euchisq20$bias_estimates$mean_bias_beta0,
#        euchisq20$bias_estimates$mean_bias_b0))


rownames(mean_bias_euchisq) <- c("euchisq1", "euchisq2", "euchisq3", "euchisq4", 
                                 "euchisq6", "euchisq7", "euchisq8", "euchisq9", 
                                 "euchisq11", "euchisq12", "euchisq13", "euchisq14", 
                                 "euchisq16", "euchisq17", "euchisq18", "euchisq19")

colnames(mean_bias_euchisq) <- c("mean_bias_ui", "mean_bias_beta0", "mean_bias_b0")

write.csv(mean_bias_euchisq, "mean_bias_euchisq.csv")

setwd("F:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\Results")
png("Bias_ui_N30_euchisq.png", width = 600, height = 500)
ggplot(bias_ui_euchisq_N30_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=30", subtitle="u_i ~ CHISQ(4), e_ij ~ CHISQ(4)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

png("Bias_ui_N100_euchisq.png", width = 600, height = 500)
ggplot(bias_ui_euchisq_N100_long) +
  geom_line(aes(y=bias, x=id, group = Setting, color = Setting), size=0.8) +
  labs(title="Bias u_i for N=100", subtitle="u_i ~ CHISQ(4), e_ij ~ CHISQ(4)") +
  scale_color_brewer(palette="Set1") + theme_minimal() 
dev.off()

