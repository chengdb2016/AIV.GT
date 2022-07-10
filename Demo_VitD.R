### VITD dataset
#setwd("./AIV.GT-code")

source("./functions/AIV.GT.R")
source("./functions/PossCEs.R")

library(pcalg)
library(ivreg)

load("./Datasets/VitD.rData")

alpha <- 0.05
# data(VitD)
VitD$Vitd_std <- (VitD$vitd -20)/20

VitD <- VitD[,-3]
Y <- VitD$death
Tr <- VitD$Vitd_std
 
dat <- cbind.data.frame(VitD[,-c(4,5)],Tr,Y)
varnames <- colnames(dat)

##### AIV.GT ####
expdat <- dat
X_idx <- 4
Y_idx <- 5
K<-2

Poss_CEs<-PossCEs(expdat,X_idx, Y_idx,alpha)
est_AIV.GT <- AIV.GT(expdat, X_idx, Y_idx, 2, Poss_CEs)
hat.beta.wy <- est_AIV.GT$ace

print(hat.beta.wy)

### AIV.GT returns a pair of valid AIVs 
varnames[c(unlist(est_AIV.GT$S[which.max(est_AIV.GT$score_sel)]))]