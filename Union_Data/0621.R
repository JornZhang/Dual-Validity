#########################################################################
########### AR1 structured Xperp and eps #################### 
########### Add aggregate shock by loading1 and loading2 to eps #########
#########################################################################
rm(list = ls())




load("1108_PermSimu_Setting_c.Rdata")


# --- the function to get the diagonal entry in the projection matrix --- #
# Input: X is a matrix
# Output: the diagonal entries of the projection matrix (X%*%solve(t(X)%*%X)%*%t(X))
leverages_chol <- function(X) {
  R <- chol(crossprod(X))                 # p x p, upper-triangular
  Z <- backsolve(R, t(X), transpose = TRUE)  # solves t(R) Z = t(X); Z is p x n
  colSums(Z^2)                            # = rowSums((X R^{-T})^2) = diag(H)
}


hac_quad_bartlett <- function(g, Lag) {
  Tn <- length(g)
  Lag <- min(Lag, Tn - 1)
  
  out <- sum(g^2)
  
  if (Lag >= 1) {
    for (ell in 1:Lag) {
      w <- 1 - ell / (Lag + 1)
      out <- out + 2 * w * sum(g[(ell + 1):Tn] * g[1:(Tn - ell)])
    }
  }
  
  as.numeric(out)
}

hac_quad_prewhite_ar1 <- function(g, Lag) {
  g <- as.numeric(g)
  Tn <- length(g)
  
  rho_hat <- coef(lm(g[-1] ~ g[-Tn] - 1))[1]
  rho_hat <- max(min(rho_hat, 0.95), -0.95)
  
  v <- g[-1] - rho_hat * g[-Tn]
  
  hac_quad_bartlett(v, Lag = min(Lag, length(v) - 1)) / (1 - rho_hat)^2
}





# ---------- the function 
simulate_ar1 <- function(n_year, n_state = 47, phi = 0.8, sigma = 1) {
  # innovations: n x n_series, iid N(0, sigma^2)
  eps <- matrix(rnorm(n_year * n_state, sd = sigma), nrow = n_year, ncol = n_state)
  
  X <- matrix(NA_real_, nrow = n_year, ncol = n_state)
  
  # initialize at stationary distribution (optional)
  X[1, ] <- eps[1, ] / sqrt(1 - phi^2)
  
  # AR(1) recursion, vectorized over 47 series
  for (t in 2:n_year) {
    X[t, ] <- phi * X[t - 1, ] + eps[t, ]
  }
  
  X  # n x 47 matrix: each column is one AR(1) process
}









run = 10000

Lag = 20


shock.seq = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)





# ====================================================================== #
# ===== loading.strength will change across different files in HPC ===== #
# ====================================================================== #
loop.s = 3 # this value varies from 1 to 10, which is read from the sbatch file.
loading.strength = shock.seq[loop.s]
# ====================================================================== #



Z.eta = Z%*%eta
Z.gamm = Z%*%gamm


n_state = ncol(eps_perp_ori_wide)
n_year = nrow(eps_perp_ori_wide)

state.unq = unique(statefip_ipt)

sd.eps = sd(c(eps_perp_ori_wide))
sd.Xperp = sd(c(Xperp_ori_wide))

# loading1 = rep( c(rep(1, times=45), rep(0, times=32)), times= n_state) #1938-1982, pre-internet
# loading2 = rep( c(rep(0, times=45), rep(1, times=32)), times= n_state) #1983-2014, post-internet


loading1 = c(loading1_wide)
loading2 = c(loading2_wide)

# one.N = rep(1, times = (n_state*n_year) )

Wgt.state = matrix(0, n_year, n_year)
# Newey-West with finite sample adjustment
for (Wgtrow in 1:n_year) {
  for (Wgtcol in 1:n_year) {
    Wgt.state[Wgtrow, Wgtcol] = max( 1 - abs(Wgtrow-Wgtcol)/(Lag+1), 0  ) 
  }
}
rm(Wgtrow, Wgtcol)
Wgt_HAC = kronecker(diag(n_state), Wgt.state)

phi_eps = mean(res_c_acf_mat$ACF_mat[2,])
phi_Xperp = mean(Xperp_c_acf_mat$ACF_mat[2,])



  
## IMPORTANT: same seed for each shock level, as in your original code
set.seed(1130)
  
## The usual t and se2 under HAC without finite sample adjustment
t_HAC  <- numeric(run)
se2_HAC <- numeric(run)

# Satterthwaite adjusted df if we use critical value from t as the critical value
df_HAC <- numeric(run) 

# use prewhitened HAC on the score, if df_HAC above still low
t_HAC_PW = numeric(run)
se2_HAC_PW = numeric(run)


  
for (loop.run in 1:run) {
    
    ## ----- permute eps_perp via AR(1) ----- ##
    eps_ar1_wide <- simulate_ar1(
      n_year  = n_year,
      n_state = n_state,
      sigma   = sd.eps,
      phi     = phi_eps
    )
    eps_ar1 <- c(eps_ar1_wide)
    
    ## add aggregate shock
    f1 <- rnorm(1)
    f2 <- rnorm(1)
    
    eps <- sqrt(1 - loading.strength) * eps_ar1 +
      sqrt(loading.strength) * sd.eps * (f1 * loading1 + f2 * loading2)
    
    ## ----- permute Xperp via AR(1) ----- ##
    Xperp_ar1_wide <- simulate_ar1(
      n_year  = n_year,
      n_state = n_state,
      sigma   = sd.Xperp,
      phi     = phi_Xperp
    )
    Xperp <- c(Xperp_ar1_wide)
    
    ## construct X and y
    X <- Z.eta + Xperp
    y <- bet * X + Z.gamm + eps
    
    ## OLS fit
    M <- cbind(X, Z)
    colnames(M)[1] <- "X"
    
    fit.ols <- lm(y ~ M - 1)
    
    ## HAC via residuals of X on Z
    fit.X    <- lm(X ~ Z - 1)
    est.Xperp <- fit.X$residuals
    est.eps.hat <- fit.ols$residuals
    
    Hardmard <- t(est.Xperp) %*%
      ((est.eps.hat %*% t(est.eps.hat)) * Wgt_HAC) %*%
      est.Xperp
    
    se2_denom <- (sum(est.Xperp * est.Xperp))^2
    est.se2   <- as.numeric(Hardmard / se2_denom)
    
    se2_HAC[loop.run] <- est.se2
    t_HAC[loop.run]   <- (coef(fit.ols)[1] - bet) / sqrt(est.se2)
    
    
    
    
    # Calculate the Satterthwaite style degree of freedom
    score <- est.Xperp * est.eps.hat
    score_wide <- matrix(score, nrow = n_year, ncol = n_state)
    
    q_state <- colSums(score_wide * (Wgt.state %*% score_wide))
    
    df_eff <- 2 * (sum(q_state)^2) / sum(q_state^2)
    df_HAC[loop.run] <- max(min(df_eff, n_state - 1), 2)
    
    
    
    
    # prewhitened HAC on the scorescore <- est.Xperp * est.eps.hat
    Hardmard_PW <- sum(apply(score_wide, 2, hac_quad_prewhite_ar1, Lag = Lag))
    est.se2_PW <- as.numeric(Hardmard_PW / se2_denom)
    
    se2_HAC_PW[loop.run] = est.se2_PW
    t_HAC_PW[loop.run] = (coef(fit.ols)[1] - bet) / sqrt(est.se2_PW)
    
    
    
    print(loop.run)
    
}
  
## return a list for this shock level
reslt =  list(t = t_HAC, 
       se2 = se2_HAC,
       df = df_HAC,
       t_PW = t_HAC_PW,
       se2_PW = se2_HAC_PW)

hist(t_HAC[1:100])
hist(t_HAC_PW[1:100])
plot(t_HAC, t_HAC_PW)
abline(a=0,b=1,col=2)
plot(se2_HAC, se2_HAC_PW)
abline(a=0,b=1,col=2)

saveRDS(reslt, file = paste0("621_sameseed_epsShock=", loading.strength, ".rds"))




