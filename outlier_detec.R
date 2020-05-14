N = 100
ni = 12
df=4
dfi = 4
meanx = 0
meanij = 0
sd = 1
sdi = 1
alpha_0 = 10
tau1 = 0.025
tau2 = 0.975

Subject <- rep(1:N, each = ni)
Replicates <- rep(1:ni, N)


########################
#for scenario 1
########################

#2 for both e_ij and u_i follow N
var_ui = exp((2*meanij)+(sdi^2))*(exp(sdi^2)-1)
sdi_sc = sqrt(var_ui)

var_e = exp((2*meanx)+(sd^2))*(exp(sd^2)-1)
sd_sc = sqrt(var_e)

u_ijs = rlnorm(N, mean=meanij, sd=sdi)
u_ij = (u_ijs - exp(meanij + ((sdi^2)/2)))/sdi_sc
u_i = rep(u_ij, each = ni)
e_ij = (rlnorm(N*ni, mean=meanx, sd=sd) - exp(meanx + ((sd^2)/2)))/sd_sc
Measure = alpha_0 + e_ij + u_i   


data = data.frame(Measure, Subject, u_i, e_ij)
length(data$Measure)

#boxplot for error
e_ij <- data.frame(e_ij)
head(e_ij)
box1 <- ggplot(e_ij, aes(e_ij)) +
  geom_boxplot() 
box1

#boxplot for random effect
u_ij <- data.frame(u_ij)
box2 <- ggplot(u_ij, aes(u_ij)) +
  geom_boxplot() + labs(x = "u_i")
box2


#fit the model using original data
fit1a <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau1,
              nK = 11, type = "normal", data = data, 
              control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-7,
                                    beta = 0.5, gamma = 1))
fit1a
fit1b <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau2,
             nK = 11, type = "normal", data = data, 
             control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-7,
                                   beta = 0.5, gamma = 1))
fit1b

#new data - error outlier removed
data2 <- data
data2$z1 <- (data2$e_ij - mean(data2$e_ij))/sd(data2$e_ij)
data2a <- data2[data2$z1 >= -2 & data2$z1 <= 2, ]
length(data2a$Measure)

box3 <- ggplot(data2a, aes(e_ij)) +
  geom_boxplot() 
box3

fit2a <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau1,
              nK = 11, type = "normal", data = data2a, 
              control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-7,
                                    beta = 0.5, gamma = 1))
fit2a

fit2b <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau2,
              nK = 11, type = "normal", data = data2a, 
              control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-7,
                                    beta = 0.5, gamma = 1))
fit2b
mean(u_ij$u_ij)
mean(data2$u_i)
#new data - u_i outlier removed
data2$z2 <- (data2$u_i - mean(data2$u_i))/sd(data2$u_i)
data2b <- data2[data2$z2 >= -2 & data2$z2 <= 2, ]
length(data2b$Measure)

library(dplyr)
u_ij2 <- distinct(data2b, u_i)
box4 <- ggplot(u_ij2, aes(u_i)) +
  geom_boxplot() 
box4

#fit the model using new data
fit2a <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau1,
              nK = 11, type = "normal", data = data2b, 
              control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-7,
                                    beta = 0.5, gamma = 1))
fit2a

fit2b <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau2,
             nK = 11, type = "normal", data = data2b, 
             control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-7,
                                   beta = 0.5, gamma = 1))
fit2b

library(ggpubr)
fig1 <- ggarrange(box1, box2, box3, box4,
                  nrow = 2, ncol = 2)
fig1


########################
#try remove IQR outliers
########################

quantile(data$e_ij)
IQR = 0.1287940 - (-0.5335007)
out2 <- 0.1287940 + (1.5*IQR)
out1 <- (-0.5335007) - (1.5*IQR)
#new data - error outlier removed
data3 <- data
data3$out <- (data2$e_ij - mean(data2$e_ij))/sd(data2$e_ij)
data3 <- data2[data3$e_ij >= out1 & data3$e_ij <= out2, ]
length(data3$Measure)

box5 <- ggplot(data3, aes(e_ij)) +
  geom_boxplot() 
box5

fit3a <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau1,
              nK = 11, type = "normal", data = data3, 
              control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-7,
                                    beta = 0.5, gamma = 1))
fit3a

fit3b <- lqmm(fixed = Measure ~ 1, random = ~ 1, group = Subject, tau = tau2,
              nK = 11, type = "normal", data = data3, 
              control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-7,
                                    beta = 0.5, gamma = 1))
fit3b

fig2 <- ggarrange(box1, box2, box5, box4,
                  nrow = 2, ncol = 2)
fig2

