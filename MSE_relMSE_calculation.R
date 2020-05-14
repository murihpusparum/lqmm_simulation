r1_1 <- alpha_0 + qnorm(tau1, mean=meanx, sd=sd)    
r1_2 <- alpha_0 + qnorm(tau2, mean=meanx, sd=sd)   
r2_1 <- alpha_0 + qnorm(tau1, mean=meanx, sd=sd)
r2_2 <- alpha_0 + qnorm(tau2, mean=meanx, sd=sd)
r3_1 <- alpha_0 + qt(tau1, df=3)
r3_2 <- alpha_0 + qt(tau2, df=3)
r4_1 <- alpha_0 + qt(tau1, df=5)
r4_2 <- alpha_0 + qt(tau2, df=5)
r5_1 <- alpha_0 + qt(tau1, df=5)
r5_2 <- alpha_0 + qt(tau2, df=5)
r6_1 <- alpha_0 + qchisq(tau1, df=4)
r6_2 <- alpha_0 + qchisq(tau2, df=4)
r7_1 <- alpha_0 + qnorm(tau1, mean=meanx, sd=sd)
r7_2 <- alpha_0 + qnorm(tau2, mean=meanx, sd=sd)
r8_1 <- alpha_0 + qchisq(tau1, df=4)
r8_2 <- alpha_0 + qchisq(tau2, df=4)
r9_1 <- alpha_0 + qnorm(tau1, mean=meanx, sd=sd)
r9_2 <- alpha_0 + qnorm(tau2, mean=meanx, sd=sd)
r10_1 <- alpha_0 + qlnorm(tau1, mean=meanx, sd=sd)
r10_2 <- alpha_0 + qlnorm(tau2, mean=meanx, sd=sd)
r11_1 <- alpha_0 + qlnorm(tau1, mean=meanx, sd=sd)
r11_2 <- alpha_0 + qlnorm(tau2, mean=meanx, sd=sd)

b0.025 <- c(r1_1, r2_1, r3_1, r4_1, r5_1, r6_1, r7_1, r8_1, r9_1, r10_1, r11_1)
b0.975 <- c(r1_2, r2_2, r3_2, r4_2, r5_2, r6_2, r7_2, r8_2, r9_2, r10_2, r11_2)
b <- data.frame(b0.025, b0.975)

setwd("D:\\BIOSTAT - UHASSELT\\Doctoral Degree\\1st Year\\QR practice\\Simulation\\15.04.20\\beta_est\\sim1")
est11_30 <- read.csv("1tau25_30.csv", sep = ",", header = T)
est21_30 <- read.csv("1tau975_30.csv", sep = ",", header = T)
est12_30 <- read.csv("2tau25_30.csv", sep = ",", header = T)
est22_30 <- read.csv("2tau975_30.csv", sep = ",", header = T)
est13_30 <- read.csv("3tau25_30.csv", sep = ",", header = T)
est23_30 <- read.csv("3tau975_30.csv", sep = ",", header = T)
est14_30 <- read.csv("4tau25_30.csv", sep = ",", header = T)
est24_30 <- read.csv("4tau975_30.csv", sep = ",", header = T)
est15_30 <- read.csv("5tau25_30.csv", sep = ",", header = T)
est25_30 <- read.csv("5tau975_30.csv", sep = ",", header = T)
est16_30 <- read.csv("6tau25_30.csv", sep = ",", header = T)
est26_30 <- read.csv("6tau975_30.csv", sep = ",", header = T)
est17_30 <- read.csv("7tau25_30.csv", sep = ",", header = T)
est27_30 <- read.csv("7tau975_30.csv", sep = ",", header = T)
est18_30 <- read.csv("8tau25_30.csv", sep = ",", header = T)
est28_30 <- read.csv("8tau975_30.csv", sep = ",", header = T)
est19_30 <- read.csv("9tau25_30.csv", sep = ",", header = T)
est29_30 <- read.csv("9tau975_30.csv", sep = ",", header = T)
est110_30 <- read.csv("10tau25_30.csv", sep = ",", header = T)
est210_30 <- read.csv("10tau975_30.csv", sep = ",", header = T)
est111_30 <- read.csv("11tau25_30.csv", sep = ",", header = T)
est211_30 <- read.csv("11tau975_30.csv", sep = ",", header = T)

##### N = 30
#r1 - 25
rb1_25 <- mean((est11_30[,2] - b[1,1])/abs(b[1,1]))
rb2_25 <- mean(abs(est11_30[,2] - b[1,1])/b[1,1])
MSE_25 <- mean((b[1,1] - est11_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est11_30[,2])/est11_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[1,1])/b[1,1])


#r1 - 975
rb1_975 <- mean((est21_30[,2] - b[1,2])/abs(b[1,2]))
rb2_975 <- mean(abs(est21_30[,2] - b[1,2])/b[1,2])
MSE_975 <- mean((b[1,2] - est21_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est21_30[,2])/est21_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[1,2])/b[1,2])

r1 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)

#r2 - 25
rb1_25 <- mean((est12_30[,2] - b[2,1])/abs(b[2,1]))
rb2_25 <- mean(abs(est12_30[,2] - b[2,1])/b[2,1])
MSE_25 <- mean((b[2,1] - est12_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est12_30[,2])/est12_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[2,1])/b[2,1])


#r2 - 975
rb1_975 <- mean((est22_30[,2] - b[2,2])/abs(b[2,2]))
rb2_975 <- mean(abs(est22_30[,2] - b[2,2])/b[2,2])
MSE_975 <- mean((b[2,2] - est22_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est22_30[,2])/est22_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[2,2])/b[2,2])

r2 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)


#r3 - 25
rb1_25 <- mean((est13_30[,2] - b[3,1])/abs(b[3,1]))
rb2_25 <- mean(abs(est13_30[,2] - b[3,1])/b[3,1])
MSE_25 <- mean((b[3,1] - est13_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est13_30[,2])/est13_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[3,1])/b[3,1])


#r3 - 975
rb1_975 <- mean((est23_30[,2] - b[3,2])/abs(b[3,2]))
rb2_975 <- mean(abs(est23_30[,2] - b[3,2])/b[3,2])
MSE_975 <- mean((b[3,2] - est23_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est23_30[,2])/est23_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[3,2])/b[3,2])

r3 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)


#r4 - 25
rb1_25 <- mean((est14_30[,2] - b[4,1])/abs(b[4,1]))
rb2_25 <- mean(abs(est14_30[,2] - b[4,1])/b[4,1])
MSE_25 <- mean((b[4,1] - est14_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est14_30[,2])/est14_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[4,1])/b[4,1])


#r4 - 975
rb1_975 <- mean((est24_30[,2] - b[4,2])/abs(b[4,2]))
rb2_975 <- mean(abs(est24_30[,2] - b[4,2])/b[4,2])
MSE_975 <- mean((b[4,2] - est24_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est24_30[,2])/est24_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[4,2])/b[4,2])

r4 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)

#r5 - 25
rb1_25 <- mean((est15_30[,2] - b[5,1])/abs(b[5,1]))
rb2_25 <- mean(abs(est15_30[,2] - b[5,1])/b[5,1])
MSE_25 <- mean((b[5,1] - est15_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est15_30[,2])/est15_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[5,1])/b[5,1])


#r5 - 975
rb1_975 <- mean((est25_30[,2] - b[5,2])/abs(b[5,2]))
rb2_975 <- mean(abs(est25_30[,2] - b[5,2])/b[5,2])
MSE_975 <- mean((b[5,2] - est25_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est25_30[,2])/est25_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[5,2])/b[5,2])

r5 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)

#r6 - 25
rb1_25 <- mean((est16_30[,2] - b[6,1])/abs(b[6,1]))
rb2_25 <- mean(abs(est16_30[,2] - b[6,1])/b[6,1])
MSE_25 <- mean((b[6,1] - est16_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est16_30[,2])/est16_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[6,1])/b[6,1])


#r6 - 975
rb1_975 <- mean((est26_30[,2] - b[6,2])/abs(b[6,2]))
rb2_975 <- mean(abs(est26_30[,2] - b[6,2])/b[6,2])
MSE_975 <- mean((b[6,2] - est26_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est26_30[,2])/est26_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[6,2])/b[6,2])

r6 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)


#r7 - 25
rb1_25 <- mean((est17_30[,2] - b[7,1])/abs(b[7,1]))
rb2_25 <- mean(abs(est17_30[,2] - b[7,1])/b[7,1])
MSE_25 <- mean((b[7,1] - est17_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est17_30[,2])/est17_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[7,1])/b[7,1])


#r7 - 975
rb1_975 <- mean((est27_30[,2] - b[7,2])/abs(b[7,2]))
rb2_975 <- mean(abs(est27_30[,2] - b[7,2])/b[7,2])
MSE_975 <- mean((b[7,2] - est27_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est27_30[,2])/est27_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[7,2])/b[7,2])

r7 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)

#r8 - 25
rb1_25 <- mean((est18_30[,2] - b[8,1])/abs(b[8,1]))
rb2_25 <- mean(abs(est18_30[,2] - b[8,1])/b[8,1])
MSE_25 <- mean((b[8,1] - est18_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est18_30[,2])/est18_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[8,1])/b[8,1])


#r8 - 975
rb1_975 <- mean((est28_30[,2] - b[8,2])/abs(b[8,2]))
rb2_975 <- mean(abs(est28_30[,2] - b[8,2])/b[8,2])
MSE_975 <- mean((b[8,2] - est28_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est28_30[,2])/est28_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[8,2])/b[8,2])

r8 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)


#r9 - 25
rb1_25 <- mean((est19_30[,2] - b[9,1])/abs(b[9,1]))
rb2_25 <- mean(abs(est19_30[,2] - b[9,1])/b[9,1])
MSE_25 <- mean((b[9,1] - est19_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est19_30[,2])/est19_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[9,1])/b[9,1])


#r9 - 975
rb1_975 <- mean((est29_30[,2] - b[9,2])/abs(b[9,2]))
rb2_975 <- mean(abs(est29_30[,2] - b[9,2])/b[9,2])
MSE_975 <- mean((b[9,2] - est29_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est29_30[,2])/est29_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[9,2])/b[9,2])

r9 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)


#r10 - 25
rb1_25 <- mean((est110_30[,2] - b[10,1])/abs(b[10,1]))
rb2_25 <- mean(abs(est110_30[,2] - b[10,1])/b[10,1])
MSE_25 <- mean((b[10,1] - est110_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est110_30[,2])/est110_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[10,1])/b[10,1])


#r10 - 975
rb1_975 <- mean((est210_30[,2] - b[10,2])/abs(b[10,2]))
rb2_975 <- mean(abs(est210_30[,2] - b[10,2])/b[10,2])
MSE_975 <- mean((b[10,2] - est210_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est210_30[,2])/est210_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[10,2])/b[10,2])

r10 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)


#r11 - 25
rb1_25 <- mean((est111_30[,2] - b[11,1])/abs(b[11,1]))
rb2_25 <- mean(abs(est111_30[,2] - b[11,1])/b[11,1])
MSE_25 <- mean((b[11,1] - est111_30[,2])^2)
rMSE1_25 <- mean(abs(sqrt(MSE_25) - est111_30[,2])/est111_30[,2])
rMSE2_25 <- mean(abs(sqrt(MSE_25) - b[11,1])/b[11,1])


#r11 - 975
rb1_975 <- mean((est211_30[,2] - b[11,2])/abs(b[11,2]))
rb2_975 <- mean(abs(est211_30[,2] - b[11,2])/b[11,2])
MSE_975 <- mean((b[11,2] - est211_30[,2])^2)
rMSE1_975 <- mean(abs(sqrt(MSE_975) - est211_30[,2])/est211_30[,2])
rMSE2_975 <- mean(abs(sqrt(MSE_975) - b[11,2])/b[11,2])

r11 <- cbind.data.frame(MSE_25, MSE_975, rMSE1_25, rMSE1_975, rMSE2_25, rMSE2_975, rb1_25, rb1_975, rb2_25, rb2_975)

r <- rbind.data.frame(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11)
colnames(r) <- c("MSE_25", "MSE_975", "relMSE_25a", "relMSE_975a", "relMSE_25b", "relMSE_975b", "relbias_25a", "relbias_975a", "relbias_25b", "relbias_975b")






MSE <- mean((b[8,2] - est28_30[,2])^2)
rel.MSE1 <- mean(abs(sqrt(MSE) - est28_30[,2])/est28_30[,2])
rel.MSE2 <- mean(abs(sqrt(MSE) - b[8,2])/b[8,2])




rel.bias1 <- mean((hatb - b)/abs(b))
rel.bias2 <- mean(sum(abs(hatb - b)/abs(b)))
