############################################################################################
# AIV.GT
#
# The main functions for AIV.GT
#
# Code by
#
#  - Anonymous Author(s)
#
# Current version: ***
 

#' @title Discovering ancestral instrumental variables for causal inference from observational data.
#'
#' @description
#' 
#' @param traindat training set.
#' @param x         index of treatment variable.
#' @param y         index of outcome variable.
#' @param K         K=2 for searching for a pair of valid IVs.
#' @est.AIVs       Possible causal effects from all candidate IVs.

AIV.GT <- function(traindat, x, y, K=2, est.AIVs, verbose=FALSE)
{

  CanS <- est.AIVs$CandidateAIVs
  Z  <- est.AIVs$ConditionalZ
  Alpah_hats <- est.AIVs$aces
  est_sigma <- est.AIVs$sigma_est
  m <- matrix(colMeans(traindat), nrow = 1)
  traindat <- traindat - m[rep(1, nrow(traindat)), ]
  train_corr <- cor(traindat)
  train_cov <- cov(traindat)

  S <- list()
  Z1 <- list();Z2 <- list()
  tetrad_tests<-list()
  tetrad_scores <- c()
  aces<-c()
  Candidates <- combn(1:length(CanS),K)

  for (i in 1:ncol(Candidates)) {
    choice <- Candidates[,i]

    S[[i]] <-CanS[choice]
    Z1[[i]]<-unlist(Z[choice[1]])
    Z2[[i]]<-unlist(Z[choice[2]])

    sigma_1 <-est_sigma[[choice[1]]]$cov_swz * est_sigma[[choice[2]]]$cov_syz
    sigma_2 <-est_sigma[[choice[2]]]$cov_swz * est_sigma[[choice[1]]]$cov_syz
    ACE_1 <- Alpah_hats[choice][1]
    ACE_2 <- Alpah_hats[choice][2]
    tetrad_score <- -abs(abs(sigma_1 - sigma_2)-abs(ACE_1 -ACE_2))

    est_ace <- mean(Alpah_hats[choice])
    aces <- c(aces,est_ace)
    tetrad_scores <- c(tetrad_scores, tetrad_score)

  }

  tetrad_tests <- Tetrad_test(traindat, x, y, S, Z1, Z2, p_value = 0.05)
  score_sel<-tetrad_scores
  score_sel[tetrad_tests$passed == FALSE] <- -Inf
  ace <- aces[which.max(score_sel)]
  if (length(aces) == 0) ace <- 0

  return(list(ace = ace, S = S, Z1= Z1,Z2=Z2, aces = aces, tetrad_scores = tetrad_scores, score_sel = score_sel, tetrad_tests = tetrad_tests))
}


sigam_est <- function(W, Z, x, y, cov_model)
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


  return(list(cov_swz=cov_RR,cov_syz = cov_Y.R))
}


Tetrad_test <- function(dat, x, y, W, Z1, Z2, p_value, verbose = FALSE)
{

  num_solutions <- length(W)
  passed <- rep(FALSE, num_solutions)
  train_dat <- standardize_data(dat)
  train_cov <- cov(train_dat)
  N <- nrow(dat)
  score_pvalues <- rep(0, num_solutions)

  for (i in seq_len(num_solutions)) {

    if (verbose) cat(sprintf("Testing set %d out of %d\n", i, num_solutions))

    tetrad_test <- tetrad_test_temp(W[[i]], Z1[[i]], Z2[[i]], x, y, train_cov, N)
    pass_tetrad <- mean(tetrad_test > p_value, na.rm = TRUE) > 0.5
    score_pvalues[i] <- mean(log(tetrad_test))
    passed[i] <- pass_tetrad

  }

  return(list(passed = passed, score_pvalues = score_pvalues))
}


# tetrad_set_test:
#
# Assess whether tetrad constraints in pair (W, Z) hold.
#
# Input:
#
# - W: candidate instruments
# - Z: adjustment set
# - x, y: treatment and outcome
# - C: empirical covariance matrix
# - N: sample size
#
# Output:
#
# - pvalues: pvalues for all pairwise tests

tetrad_test_temp <- function(W, Z1, Z2, x, y, C, N)
{
  p <- length(W)
  if (p < 2) {
    return
  }
  w_count <- 0
  pvalues <- rep(0, p * (p - 1) / 2)

  for (w1i in 1:(p - 1)) {
    w1 <- W[w1i]
    for (w2i in (w1i + 1):p) {
      w2 <- W[w2i]
      w_count <- w_count + 1
      pvalues[w_count] <- tetrad_test_wishart(C, N, w1, w2, x, y, Z1, Z2)
    }
  }

  return(pvalues)
}

# tetrad_test_wishart:

tetrad_test_wishart <- function(C, N, S1, S2, x, y, Z1, Z2)
{
  xy  <- c(S1, S2, x, y)
  Z <- setdiff(unique(c(Z1,Z2)),xy)
  if (length(Z) > 0) {
    cxy <- C[xy, xy] - C[xy, Z, drop = FALSE] %*% solve(C[Z, Z, drop = FALSE], C[Z, xy, drop = FALSE])
  } else {
    cxy <- C[xy, xy]
  }
  D_12 <- cxy[1, 1] * cxy[2, 2] - cxy[1, 2]^2
  D_34  <- cxy[3, 3] * cxy[4, 4] - cxy[3, 4]^2
  D  <- det(cxy)
  std  <- sqrt(D_12  * D_34  * (N + 1) / ((N - 1) * (N - 2)) - D  / (N - 2))

  v <- cxy[1, 3] * cxy[2, 4] - cxy[1, 4] * cxy[2, 3]
  p_value <- 2 * pnorm(-abs(v), 0, std)

  return(p_value)
}

