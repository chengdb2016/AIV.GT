# ############################################################################################
# PossCEs.R
#
# Code by
#
#  - Anonymous Author(s)
#
# Current version: ***


#' @title Discovering ancestral instrumental variables for causal inference from observational data.
#'
#'
#' @description
#'
#'
#' @param train_dat training set.
#' @param x         index of treatment variable.
#' @param y         index of outcome variable.
#' @param alpha     the signficant level, default is 0.05
#' @param verbose   if TRUE, prints out information about the steps taken.

PossCEs<-function(traindat, x, y, alpha=0.05, verbose=FALSE){
  
  require(pcalg)
  
  m <- matrix(colMeans(traindat), nrow = 1)
  traindat <- expdat - m[rep(1, nrow(traindat)), ]
  train_corr <- cor(traindat)
  train_cov <- cov(traindat)
  
  p <- ncol(traindat)
  Varnames <- colnames(traindat)
  suffStat <- list(C=cor(traindat),n=nrow(traindat))
  fci.est <- rfci(suffStat, indepTest = gaussCItest, labels = colnames(traindat), alpha=alpha)
  
  AdjW <- setdiff(which(fci.est@amat[,x]!=0),y)
  AdjY <- setdiff(which(fci.est@amat[,y]!=0),x)
  
  CanS <- union(AdjW,AdjY)
  Z <- list()
  Alpah_hats<-c()
  est_sigma<-list()
  if(length(CanS)==0|length(CanS)==1){
    print("The asssumption of at least two valid IVs does not hold.")
  }else{
    for (i in 1:length(CanS)) {
      
      S <- CanS[i]
      Z_y <- setdiff(possAn(fci.est@amat,x=y, possible = TRUE, ds = FALSE,type = "pag"), c(x,y,S))
      Z_iv <- setdiff(possAn(fci.est@amat,x=S, possible = TRUE, ds = FALSE,type = "pag"),c(x,y,S))
      z_temp <-unique(c(Z_y,Z_iv))
      
      if (verbose) {
        cat(sprintf("--------- SEARCH for conditioning set Z that instrumentalises S to be an ancestral IV, \n", z_temp))
      }
      
      Z[[i]]<-z_temp
      
      Alpah_hat  <- IV_TSLS(S, z_temp,x,y,train_cov)
      est_sigma[[i]] <- sigam_est(S, z_temp,x,y,train_cov)
      
      Alpah_hats <-c(Alpah_hats,Alpah_hat)
    }
  }
  
  return(list(CandidateAIVs = CanS, ConditionalZ = Z, aces = Alpah_hats, sigma_est = est_sigma, fcI_est = fci.est))
  
}

#' @title Estimate differential causal effect using two-stage least squares.
#'
#' @description
#'
#' This method is the simplest variation of two-stage least squares.
#'
#' @param W         array containing indices of the instrumental variables.
#' @param Z         array containing indices of the conditioning variables.
#' @param x         index of treatment variable.
#' @param y         index of outcome variable.
#' @param cov_model covariance matrix of all variables
#'
#' @return Differential causal effect as implied by \code{cov_model}.

IV_TSLS <- function(W, Z, x, y, cov_model)
{
  WZ <- c(W, Z)
  XZ <- c(x, Z)
  
  # First stage
  
  beta_x <- solve(cov_model[WZ, WZ], cov_model[WZ, x])
  
  # Build covariance of X, Z, R and Y
  
  cov_XZ.R <- cov_model[c(x, Z), x] - colSums(beta_x * cov_model[WZ, XZ, drop = FALSE])
  cov_RR   <- cov_model[x, x] - sum(cov_model[x, WZ] * beta_x)
  cov_Y.R  <- cov_model[y, x] - sum(beta_x * cov_model[y, WZ])
  
  cov_XZR   <- rbind(cbind(cov_model[XZ, XZ], cov_XZ.R), c(cov_XZ.R, cov_RR))
  cov_XZR.Y <- c(cov_model[XZ, y], cov_Y.R)
  
  # Second stage
  
  beta_y <- solve(cov_XZR, cov_XZR.Y)
  
  # Return
  
  return(as.numeric(beta_y[1]))
}

standardize_data <- function(dat)
{
  N <- nrow(dat)
  sd_dat <- apply(dat, 2, sd)
  m_dat <- colMeans(dat)
  dat <- (dat - t(matrix(rep(m_dat, N), nrow = length(m_dat)))) / t(matrix(rep(sd_dat, N), nrow = length(sd_dat)))
  return(dat)
}

