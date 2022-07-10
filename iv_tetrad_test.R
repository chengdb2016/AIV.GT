############################################################################################
# iv_tetrad.R
#
# The main functions for IV-TETRAD
#
# Code by
#
#  - Ricardo Silva (ricardo@stats.ucl.ac.uk)
#
# Current version: 18/10/2016
# First version: 18/10/2016

#' @title Estimates causal effects by searching for plausible instruments.
#' 
#' @description
#' 
#' Generates synthetic examples and evaluation of different algorithms for estimating
#' causal effects in linear models using instrumental variables.
#' 
#' @param train_dat training set.
#' @param W_space   indices of the candidate instruments within the training set.
#' @param x         index of treatment variable.
#' @param y         index of outcome variable.
#' @param K_win     window of candidates as used by the algorithm.
#' @param verbose   if TRUE, prints out information about the steps taken.
#'
#'
#' @return A list containing several outcomes of interest:
#'   \item{\code{ace}}{a suggested average causal effect.}
#'   \item{\code{W}}{a list of candidate instruments.}
#'   \item{\code{Z}}{a list of corresponding candidate adjustment sets.}
#'   \item{\code{aces}}{a list of corresponding candidate causal effects.}
#'   \item{\code{tetrad_scores}}{a list of corresponding scores for each candidate causal effect.}   
#'   \item{\code{score_sel}}{the same scores, but with \code{-Inf} for those which do not pass a tetrad test.}
#'   \item{\code{tetrad_tests]}{the corresponding tetrad tests.}
#'   
#' @details
#' This provides a full example on how to run a simulation study using the main method
#' developed in this package, along with a comparison against standard procedures. See
#' the description of \code{run_experiment} for more details.
#'  
#' @export

iv_tetrad <- function(train_dat, W_space, x, y, K_win, verbose = FALSE, blanket_prune = TRUE)
{
  m <- matrix(colMeans(train_dat), nrow = 1)
  train_dat <- train_dat - m[rep(1, nrow(train_dat)), ]
  train_corr <- cor(train_dat)
  train_cov <- cov(train_dat)
  
  # Stage 1: remove weakly associated candidates
  
  if (blanket_prune) {
    W_remove <- yx_blanket(W_space, x, y, standardize_data(train_dat))
  } else {
    W_remove <- list(remove_x = c(), remove_y = c())
  }
  if (verbose) {
    cat(sprintf("Removing %d candidate instruments\n", length(unique(c(W_remove$remove_x, W_remove$remove_y)))))
  }
  if (K_win > length(W_space)) {
    K_win <- length(W_space)
    W <- W_space
    Z <- c()
  }
  
  # Stage 2: swipe W_space to maximize tetrads
  
  num_swipes <- 0
  W <- list()
  Z <- list()
  tetrad_scores <- c()
  alpha_hat <- get_alpha_hat(W_space, x, y, train_corr)
  candidates <- setdiff(W_space, c(W_remove$remove_x, W_remove$remove_y, W_space[which(is.infinite(alpha_hat))], W_space[which(is.nan(alpha_hat))]))
  
  while (length(candidates) >= K_win) {
    
    num_swipes <- num_swipes + 1
    
    if (verbose) {
      cat(sprintf("--------- TETRAD SEARCH, SWIPE #%d\n", num_swipes))
    }
    
    alpha_out <- alpha_search(K_win, W_space, x, y, train_corr, alpha_hat, candidates, TRUE)
    W[[num_swipes]] <- alpha_out$W
    Z[[num_swipes]] <- setdiff(W_space, W[[num_swipes]])
    tetrad_scores <- c(tetrad_scores, alpha_out$tetrad_score)
    candidates <- setdiff(candidates, alpha_out$W)
    
  }
  
  # Return
  
  aces <- rep(0, num_swipes)
  score_sel <- rep(0, num_swipes)
  tetrad_tests <- solution_test(train_dat, x, y, W, Z, p_value = 0.05)
  for (s in seq_len(num_swipes)) {
    aces[s] <- IV_TSLS(W[[s]], Z[[s]], x, y, train_cov)
    score_sel[s] <- tetrad_scores[s]
  }
  score_sel[tetrad_tests$passed == FALSE] <- -Inf
  ace <- aces[which.max(score_sel)]
  if (length(aces) == 0) ace <- 0
  
  return(list(ace = ace, W = W, Z = Z, aces = aces, tetrad_scores = tetrad_scores, score_sel = score_sel, tetrad_tests = tetrad_tests))  
}

# yx_blanket:
#
# Regression-based detection of Markov blanket of x and y among W_space.
#
# Input:
#
# - W_space: indices of candidate variables to include in the Markov blanket of x and y
# - x, y: index of treatment and outcome
# - dat: the data
#
# Output:
#
# - b_x, b_y: regression coefficients of the regression of x and y on W_space
# - remove_x, remove_y; index of variables to be removed from Markov blanket of x and y among W_space

yx_blanket <- function(W_space, x, y, dat)
{
  in_dat <- dat[, W_space]
  out_dat <- dat[, x]
  w_lars <- lars(in_dat, out_dat, type = "lasso", intercept = FALSE)
  w_cv <- cv.lars(in_dat, out_dat, plot.it = FALSE, intercept = FALSE)
  frac <- w_cv$index[which.min(w_cv$cv)]
  b_x <- predict(w_lars, type = "coefficients", mode = "fraction", s = frac)$coefficients
  
  out_dat <- dat[, y]
  w_lars <- lars(in_dat, out_dat, type = "lasso", intercept = FALSE)
  w_cv <- cv.lars(in_dat, out_dat, plot.it = FALSE, intercept = FALSE)
  frac <- w_cv$index[which.min(w_cv$cv)]
  b_y <- predict(w_lars, type = "coefficients", mode = "fraction", s = frac)$coefficients
  
  remove_x <- W_space[which(b_x == 0)]
  remove_y <- W_space[which(b_y == 0)]
  for (w in W_space) {
    if (cor.test(dat[, w], dat[, x])$p.value > 0.05) {
      remove_x <- c(remove_x, w)
    }
    if (cor.test(dat[, w], dat[, y])$p.value > 0.05) {
      remove_y <- c(remove_y, w)
    }
  }
  return(list(b_x = b_x, b_y = b_y, remove_x = unique(remove_x), remove_y = unique(remove_y)))  
}

# get_alpha_hat:
#
# Causal effect estimates for each variable assuming it to be an instrument, conditional
# on everybody else.
#
# Input:
#
# - W_space: indices of candidate variables to include in the Markov blanket of x and y
# - x, y: index of treatment and outcome
# - train_cov: covariance matrix of the data
#
# Output:
#
# - alpha_hat: the causal effect estimates.

get_alpha_hat <- function(W_space, x, y, train_cov)
{
  b_yx <- solve(train_cov[W_space, W_space], train_cov[W_space, c(y, x)])
  alpha_hat <- b_yx[, 1] / b_yx[, 2]
  return(alpha_hat)  
}

# alpha_search:
#
# Search for a set of variables that best correspond to an agreement on the estimated
# causal effect.
#
# Input:
#
# - K_win: initial size of the set to be found
# - W_space: indices of candidate variables to include in the Markov blanket of x and y
# - x, y: index of treatment and outcome
# - train_corr: correlation matrix of the data
# - pre_candidates: subset of W_space which can be considered. If NULL, this is the same as W_space
# - grow_it: if TRUE, continue to grow initial set of size K_win by extending it
#
# Output:
#
# - W: the subset found
# - var_score: variance of the set of causal effects corresponding to W
# - tetrad_score: tetrad score (simple_tetrad_score) corresponding to W

alpha_search <- function(K_win, W_space, x, y, train_corr, alpha_hat, pre_candidates = NULL, grow_it = TRUE)
{
  # Set filter
  
  if (is.null(pre_candidates)) {
    pre_candidates <- W_space
  }
  V_forbidden <- rep(1, ncol(train_corr))
  V_forbidden[pre_candidates] <- 0
  W_forbidden <- which(V_forbidden[W_space] == 1)
  
  # Swipe procedure
  
  alpha_sort <- sort(alpha_hat, index.return = TRUE)  
  candidates <- setdiff(alpha_sort$ix, W_forbidden)    
  
  best_tetrad_score <- -Inf
  best_choice <- c()
  
  for (i in 1:(length(candidates) - K_win + 1)) {
    choice <- candidates[i:(i + K_win - 1)]
    tetrad_score <- simple_tetrad_score(W_space[choice], W_space[-choice], x, y, alpha_hat[choice], train_corr)
    if (tetrad_score > best_tetrad_score) {
      best_tetrad_score <- tetrad_score
      best_choice <- choice
      start_choice <- i
      end_choice <- i + K_win - 1
    }
  }
  
  # Found best choice, now grow it
  
  if (grow_it) {
    
    while (TRUE) {
      
      best_score <- simple_tetrad_score(W_space[best_choice], W_space[-best_choice], x, y, alpha_hat[best_choice], train_corr)
      if (is.infinite((best_score))) {
        break
      }
      
      score_i <- -Inf
      if (start_choice > 1) {
        choice <- candidates[(start_choice - 1):end_choice]
        score_i <- simple_tetrad_score(W_space[choice], W_space[-choice], x, y, alpha_hat[choice], train_corr)
      }
      score_j <- -Inf
      if (end_choice < length(candidates)) {
        choice <- candidates[start_choice:(end_choice + 1)]
        score_j <- simple_tetrad_score(W_space[choice], W_space[-choice], x, y, alpha_hat[choice], train_corr)
      }
      
      if (best_score > max(c(score_i, score_j))) {
        break
      }
      
      if (score_i > score_j) {
        start_choice <- start_choice - 1
        best_score <- score_i
      } else {
        end_choice <- end_choice + 1
        best_score <- score_j
      }
      best_choice <- candidates[start_choice:end_choice]
      
    }
    
  }
  
  # Return
  
  W <- W_space[best_choice]  
  var_score <- var(alpha_hat[best_choice])
  
  return(list(W = W, var_score = var_score, tetrad_score = best_tetrad_score, alpha_hat = alpha_hat, alpha_idx = best_choice))
}

# simple_tetrad_score:
#
# A measure of agreement of how the causal effects implied by a set W conditioned on Z.
#
# Input:
# 
# - W: set of candidate IVs
# - Z: conditioning set
# - x, y: treatment and outcome
# - alpha_hat: candidate causal effects
# - train_cov: covariance matrix of the data
#
# Output:
#
# - tetrad_score: the measure of agreement

simple_tetrad_score <- function(W, Z, x, y, alpha_hat, train_cov)
{
  tetrad_score <- -median(abs(IV_TSLS(W, Z, x, y, train_cov) - alpha_hat))
  return(tetrad_score) 
}

# solution_test:
#
# Assess the fitness of the tetrad/residual independence constraints of a particular set.
#
# Input:
#
# - dat: training data
# - x, y: treatment and outcome
# - W: candidate instruments
# - Z: conditioning set
# - p_value: p value for assessing test
# - res_test_it: if TRUE, test residual independence; if FALSE, test tetrads
# - verbose: if TRUE, show steps of the testing
#
# Output:
#
# - passed: TRUE, if corresponding set (W, Z) passes the implied tetrad; FALSE otherwise
# - score_pvalues: corresponding pvaluess

solution_test <- function(dat, x, y, W, Z, p_value, res_test_it = FALSE, verbose = FALSE)
{
  num_solutions <- length(W)
  passed <- rep(FALSE, num_solutions)  
  train_dat <- standardize_data(dat)
  train_cov <- cov(train_dat)
  N <- nrow(dat)
  score_pvalues <- rep(0, num_solutions)
  
  for (i in seq_len(num_solutions)) {
    
    if (verbose) cat(sprintf("Testing set %d out of %d\n", i, num_solutions))
    
    tetrad_test <- tetrad_set_test(W[[i]], Z[[i]], x, y, train_cov, N)
    pass_tetrad <- mean(tetrad_test > p_value, na.rm = TRUE) > 0.5
    score_pvalues[i] <- mean(log(tetrad_test))
    passed[i] <- pass_tetrad
    if (!pass_tetrad) next
    
    if (res_test_it) {
      res_test <- residual_set_test(W[[i]], Z[[i]], x, train_dat, train_cov, type = 2, verbose = verbose)
      pass_lingam <- mean(res_test > p_value, na.rm = TRUE) > 0.5
      passed[i] <- pass_lingam
    }
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

tetrad_set_test <- function(W, Z, x, y, C, N)
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
      pvalues[w_count] <- test_tetrad_wishart(C, N, w1, w2, x, y, Z)
    }
  }
  
  return(pvalues)
}

# test_tetrad_wishart:

test_tetrad_wishart <- function(C, N, w1, w2, y1, y2, Z)
{  
  xy  <- c(w1, w2, y1, y2)
  if (length(Z) > 0) {
    cxy <- C[xy, xy] - C[xy, Z, drop = FALSE] %*% solve(C[Z, Z, drop = FALSE], C[Z, xy, drop = FALSE])
  } else {
    cxy <- C[xy, xy]
  }    
  
  D_12 <- cxy[1, 1] * cxy[2, 2] - cxy[1, 2]^2
  D_34 <- cxy[3, 3] * cxy[4, 4] - cxy[3, 4]^2
  D <- det(cxy)
  std <- sqrt(D_12 * D_34 * (N + 1) / ((N - 1) * (N - 2)) - D / (N - 2))
  v <- cxy[1, 3] * cxy[2, 4] - cxy[1, 4] * cxy[2, 3]
  p_value <- 2 * pnorm(-abs(v), 0, std)
  
  return(p_value)
}

# tetrad_value:

tetrad_value <- function(C, x1, x2, y1, y2, Z)
{  
  xy  <- c(x1, x2, y1, y2)
  if (length(Z) > 0) {
    cxy <- C[xy, xy] - C[xy, Z, drop = FALSE] %*% solve(C[Z, Z, drop = FALSE], C[Z, xy, drop = FALSE])
  } else {
    cxy <- C[xy, xy]
  }
  v <- cxy[1, 3] * cxy[2, 4] - cxy[1, 4] * cxy[2, 3]
  return(v)
}

# residual_set_test:

residual_set_test <- function(W, Z, x, dat, type = 1, verbose = FALSE)
{
  N <- min(5000, nrow(dat))
  XX <- cov(dat)
  
  if (length(Z) + length(W) > 1) {
    res_W <- matrix(ncol = length(W), nrow = nrow(dat))
    for (w in seq_along(W)) {
      ZW <- c(Z, W[-w])
      beta_W <- solve(XX[ZW, ZW, drop = FALSE], XX[ZW, W[w], drop = FALSE])      
      res_W[, w] <- dat[, W[w]] - dat[, ZW, drop = FALSE] %*% beta_W
    }
  } else {
    res_W <- dat[, W]
  }
  
  ZW <- c(Z, W)
  beta_x <- solve(XX[ZW, ZW, drop = FALSE], XX[ZW, x, drop = FALSE])
  res_x <- dat[, x] - dat[, ZW, drop = FALSE] %*% beta_x      
  test_result <- rep(0, length(W))
  
  for (w in seq_along(W)) {
    if (shapiro.test(res_W[1:N, w])$p.value > 0.05) {
      test_result[w] <- NA
    } else {
      if (type == 1) {
        test_result[w] <- hoeffd(res_W[, w], res_x)$P[2]
      } else if (type == 2) {
        hsic_result <- dhsic.test(res_W[, w], res_x)#, bandwidth_choice = 0.001)
        if (verbose) cat(sprintf("   HSIC RESULT = %1.2f\n", hsic_result$p.value))
        test_result[w] <- hsic_result$p.value
      } else if (type == 3) {
        ps <- rep(0, 8)
        ps[1] <- cor.test(tanh(res_x), res_W[, w])$p.value
        ps[2] <- cor.test(res_x, tanh(res_W[, w]))$p.value
        ps[3] <- cor.test(res_x^2, res_W[, w])$p.value
        ps[4] <- cor.test(res_x, res_W[, w]^2)$p.value
        ps[5] <- cor.test(tanh(res_x), tanh(res_W[, w]))$p.value
        ps[6] <- cor.test(tanh(res_x), res_W[, w]^2)$p.value
        ps[7] <- cor.test(res_x^2, tanh(res_W[, w]))$p.value
        ps[8] <- cor.test(res_x^2, res_W[, w]^2)$p.value
        test_result[w] <- min(ps)
      }
    }
  }
  
  return(test_result)
}

# lingam_iv_test:
#
# Assess whether tetrad constraints in pair (W, Z) hold.
#
# Input:
#
# - dat: the data
# - W: candidate instruments
# - Z: adjustment set
# - x, y: treatment and outcome
# - p_value: threshold under which we reject tests
# - verbose: if TRUE, show details of the testing
#
# Output:
#
# - residual_tests: individual test results
# - test_summary: proportion of independence hypothesis that are not rejected

lingam_iv_test <- function(dat, W, Z, x, y, p_value, verbose = TRUE)
{
  res_test <- residual_set_test(W, Z, x, dat, type = 2, verbose = verbose)
  test_summary <- mean(res_test > p_value, na.rm = TRUE)  
  return(list(residual_tests = res_test, test_summary = test_summary))
}

