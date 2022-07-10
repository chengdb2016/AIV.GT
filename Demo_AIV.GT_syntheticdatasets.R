### Synthetic datasets
#setwd("./AIV.GT-code") 
source("./functions/AIV.GT.R")
source("./functions/PossCEs.R")

library(pcalg)
library(ivreg)

#loading different synthetic dataset
load("./Datasets/Synthetic_a.Rdata")

alpha <- 0.05

##### AIV.GT ####
varnames<-colnames(expdat)

X_idx <- 3
Y_idx <- 4
K<-2

Poss_CEs<-PossCEs(expdat,X_idx, Y_idx,alpha)
est_AIV.GT <- AIV.GT(expdat, X_idx, Y_idx, 2, Poss_CEs)
hat.beta.wy <- est_AIV.GT$ace

print(hat.beta.wy)

Bias_AIV.GT <- (abs(hat.beta.wy-2)/2)*100
cat(sprintf("The estimated Bias %f \n", Bias_AIV.GT))
### AIV.GT returns a pair of valid AIVs 
varnames[c(unlist(est_AIV.GT$S[which.max(est_AIV.GT$score_sel)]))]

source("./functions/iv_tetrad_test.R")
n <- nrow(expdat)
dat_mat <- as.matrix(expdat)
W_scope <- setdiff(1:ncol(expdat), c(X_idx, Y_idx))

num_ivs<-3

tetrad_result <- iv_tetrad(dat_mat, W_scope, X_idx, Y_idx, num_ivs, blanket_prune = FALSE)
ACE.IV.tetrad <- tetrad_result$ace

Bias_IV.tetrad <- (abs(ACE.IV.tetrad-2)/2)*100
cat(sprintf("The estimated Bias %f \n", Bias_IV.tetrad))

#####sisVIVE
library(sisVIVE)
result = sisVIVE(Y = dat_mat[,Y_idx], D = dat_mat[,X_idx], Z = dat_mat[, W_scope])
ACE.sisVIVE <- mean(result$beta)
 
Bias_sisVIVE <- (abs(ACE.sisVIVE-2)/2)*100

cov_dat <- cov(expdat)
ACE.TSLS <- IV_TSLS(W_scope, c(), X_idx, Y_idx, cov_dat)
ACE.LSR <- solve(cov_dat[c(X_idx, W_scope), c(X_idx, W_scope)], cov_dat[c(X_idx, W_scope), Y_idx])[1]

Bias_TSLS <- (abs(ACE.TSLS-2)/2)*100
Bias_LSR <- (abs(ACE.LSR-2)/2)*100

Bias_all_estimators <- cbind.data.frame(Bias_LSR,Bias_TSLS,Bias_sisVIVE,Bias_IV.tetrad,Bias_AIV.GT)




