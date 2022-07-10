
source("./functions/AIV.GT.R")
source("./functions/PossCEs.R")
load("./Datasets/SchoolingReturns.rData")

############ AIV.GT #########
alpha <- 0.05
 
Varnames <- colnames(expdat)

X_idx <- 18
Y_idx <- 19

Poss_CEs<-PossCEs(expdat,X_idx, Y_idx,alpha)
est_AIV.GT <- AIV.GT(expdat, X_idx, Y_idx, 2, Poss_CEs)
hat.beta.wy <- est_AIV.GT$ace

##hat.beta.wy
print(hat.beta.wy)

### AIV.GT returns a pair of valid AIVs 
Varnames[c(unlist(est_AIV.GT$S[which.max(est_AIV.GT$score_sel)]))]

