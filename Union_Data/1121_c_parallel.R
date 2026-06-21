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









# ---------  cluster -------- #
library(doParallel)
library(foreach)

## use 10 cores (or min(10, detectCores()))
detectCores()
n_cores <- 10L
cl <- makeCluster(n_cores)
registerDoParallel(cl)




















# ----- start simulation --- #

run = 10000

Lag = 1

# shock.seq = seq(0,0.9, by=0.1)

shock.seq = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

Z.eta = Z%*%eta
Z.gamm = Z%*%gamm


n_state = ncol(eps_perp_ori_wide)
n_year = nrow(eps_perp_ori_wide)

state.unq = unique(statefip_ipt)

sd.eps = sd(c(eps_perp_ori_wide))
sd.Xperp = sd(c(Xperp_ori_wide))

loading1 = rep( c(rep(1, times=45), rep(0, times=32)), times= n_state) #1938-1982, pre-internet
loading2 = rep( c(rep(0, times=45), rep(1, times=32)), times= n_state) #1983-2014, post-internet

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


# ---------- parallel computing the simulation of 10 shock strengths --------------------------- #
res_list <- foreach(loop.s = seq_along(shock.seq), .export = c("simulate_ar1")) %dopar% {
  
  loading.strength <- shock.seq[loop.s]
  
  ## IMPORTANT: same seed for each shock level, as in your original code
  set.seed(1130)
  
  ## local storage for this shock level
  t_HAC_col  <- numeric(run)
  se2_HAC_col <- numeric(run)
  
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
    
    se2_HAC_col[loop.run] <- est.se2
    t_HAC_col[loop.run]   <- (coef(fit.ols)[1] - bet) / sqrt(est.se2)
    
    
  }
  
  ## return a list for this shock level
  list(t = t_HAC_col, se2 = se2_HAC_col)
}


## allocate final matrices
t_HAC  <- matrix(NA_real_, nrow = run, ncol = length(shock.seq))
se2_HAC <- matrix(NA_real_, nrow = run, ncol = length(shock.seq))

for (s in seq_along(shock.seq)) {
  t_HAC[,  s] <- res_list[[s]]$t
  se2_HAC[, s] <- res_list[[s]]$se2
}

stopCluster(cl)




save(list = c(
  "dat_target", "eps_perp_ori_wide", 
  "loading1_wide", "loading2_wide",
  "Q_p_res_c", "Q_p_Xperp_c",
  "Q_res_c", "Q_Xperp_c",
  "res_c_acf_mat", "res_c_pacf_mat",
  "se2_HAC", 
  "t_HAC", 
  "Xperp_c_acf_mat", "Xperp_c_pacf_mat",
  "Xperp_ori_wide",
  "Z", "Z.eta", "Z.gamm", 
  "bet", "eta", "gamm", "Lag",
  "loading1", "loading2",
  "n_state", "n_year",
  "northeast_states", "run", "shock.seq", "south_states",
  "state.unq", "statefip_ipt", "west_states", "year_ipt"
), 
file = "1121_c_AR1X&eps_epsShock_changeShock_L1_sameseed.Rdata"
)



# make plot about non-coverage prob #------------------
# --- HAC --- #
qtl = qnorm(c(0.15, 0.1,0.05, 0.025, 0.005), lower.tail = F) # non-coverage prob of 0.3, 0.2, 0.1, 0.05, 0.01

ncp = matrix(0, nrow = length(shock.seq), ncol = length(qtl))
rownames(ncp) = shock.seq
colnames(ncp) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.s in 1:length(shock.seq)) {
    
    ncp[loop.s, loop.qtl] = sum( abs(t_HAC[,loop.s]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.s, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=shock.seq, y=ncp[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp[,loop.qtl]),
       xlab = "shock strength",
       ylab = "Non-coverage probability HAC",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}






# --- qqplot --- #
prob = seq(0.01,0.99, by=0.01)

for (loop.s in 1:length(shock.seq)) {
  qt.emp = quantile(t_HAC[,loop.s], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("shock = ", shock.seq[loop.s]))
  abline(a=0, b=1, col=2, lwd=3)
}






# --- Check the size of standard error --- #
plot(x=shock.seq,
     y=apply(sqrt(se2_HAC), 2, mean),
     type = "b",
     pch=19,
     xlab = "shock",
     ylab = " ",
     main = "HAC Standard Error")



boxplot(t_HAC)

boxplot(se2_HAC)

boxplot(t_HAC*sqrt(se2_HAC))





# HC1 --------------
rm(list = ls())
load("1121_c_AR1X&eps_epsShock_changeShock_L20_sameseed.Rdata")
N = nrow(Z)
k = ncol(Z) + 1

t_HAC_HC1 = t_HAC * sqrt( (N - k) / N )

# make plot about non-coverage prob #------------------
# --- HAC --- #
qtl = qnorm(c(0.025), lower.tail = F) # non-coverage prob of 0.3, 0.2, 0.1, 0.05, 0.01

ncp = matrix(0, nrow = length(shock.seq), ncol = length(qtl))
rownames(ncp) = shock.seq
colnames(ncp) = c("nominal=0.05")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.s in 1:length(shock.seq)) {
    
    ncp[loop.s, loop.qtl] = sum( abs(t_HAC_HC1[,loop.s]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.s, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=shock.seq, y=ncp[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp[,loop.qtl]),
       xlab = "shock strength",
       ylab = "Non-coverage probability HAC",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}
