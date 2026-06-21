
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









# ----- start simulation --- #

run = 10000

Lag = 30

# shock.seq = seq(0,0.9, by=0.1)

shock.seq = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

t_OLS = matrix(0, nrow = run, ncol = length(shock.seq))

t_HAC = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC = matrix(0, nrow = run, ncol = length(shock.seq))

# t_HAC_HC2 = matrix(0, nrow = run, ncol = length(shock.seq))
# se2_HAC_HC2 = matrix(0, nrow = run, ncol = length(shock.seq))

# t_HAC_HC3 = matrix(0, nrow = run, ncol = length(shock.seq))
# se2_HAC_HC3 = matrix(0, nrow = run, ncol = length(shock.seq))


#Z = Z[,1:19]
#eta = eta[1:19]
#gamm = gamm[1:19]

# Z = as.matrix(rep(1, times = (77*47)))
# eta= eta[1]
# gamm = gamm[1]

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

set.seed(1130)
for (loop.s in 1:length(shock.seq)) {
  
  loading.strength = shock.seq[loop.s]
  
  
  
  for (loop.run in 1:run) {
    
    # permute epsperp #---------------------------------------
    eps_ar1_wide = simulate_ar1(n_year = 77, n_state = 47, sigma = sd.eps, phi = phi_eps)
    eps_ar1 = c(eps_ar1_wide)
    
    
    # add aggregate shock to generate eps
    f1 = rnorm(1)
    f2 = rnorm(1)
    
    eps =  sqrt(1 - loading.strength) * eps_ar1 + sqrt(loading.strength) * sd.eps * ( f1*loading1 + f2*loading2 ) 
    
    
    # f = rnorm(1)
    
    # eps =  sqrt(1 - loading.strength) * eps + sqrt(loading.strength) * sd.eps * ( f*one.N ) 
    
    # dim(eps_ar1_wide)
    # plot(eps_ar1_wide[,3], type = "b")
    # plot(eps) 
    # sd(eps); sd.eps
    
    
    # permute Xperp
    Xperp_ar1_wide = simulate_ar1(n_year = 77, n_state = 47, sigma = sd.Xperp, phi = phi_Xperp)
    
    # dim(Xperp_ar1_wide)
    # plot(Xperp_ar1_wide[,47], type = "b")
    
    Xperp = c(Xperp_ar1_wide)
    
    # plot(Xperp); sd(Xperp); sd.Xperp
    
    X  = Z.eta + Xperp 
    
    # plot(X)
    
    # generate y #----------------------------------------
    y = bet * X +  Z.gamm + eps
    
    
    M = cbind(X, Z)
    colnames(M)[1] = "X"
    
    fit.ols=lm(y~M-1)
    
    t_OLS[loop.run,loop.s] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / summary(fit.ols)$coefficients[1, "Std. Error"]
    
    # HAC using X groups (too slow) #-------------------------------------------------------
    # vcov.hac = solve(t(M)%*%M)%*%(t(M)%*%((fit.ols$residuals%*%t(fit.ols$residuals))*Wgt_HAC)%*%M)%*%solve(t(M)%*%M)
    # t_HAC[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(vcov.hac[1,1])
    # se2_HAC[loop.run, loop.L] = vcov.hac[1,1]
    #---------------------------------------------------------------------------
    
    
    # --- HAC --- #
    fit.X = lm(X~Z-1)
    
    est.Xperp = fit.X$residuals
    
    # max(abs(est.Xperp - Xperp))
    # plot(est.Xperp - Xperp)
    # plot(est.Xperp); plot(Xperp)
    
    est.eps.hat = fit.ols$residuals
    
    # sum(est.eps.hat^2)
    
    Hardmard = t(est.Xperp)%*% ( ( est.eps.hat%*%t(est.eps.hat) ) * Wgt_HAC ) %*% est.Xperp
    
    se2_denom = ( sum(est.Xperp*est.Xperp) )^2
    
    est.se2 = Hardmard / se2_denom
    
    # est.se2 - vcov.hac[1,1]
    
    se2_HAC[loop.run, loop.s] = est.se2
    t_HAC[loop.run,loop.s] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2)
    # -------------------------------------------------------------# 
    
    
    
    
    print(loop.run)
    
  }
  
  print(loop.s)
  
}


save(list = c(
  "dat_target", "eps_perp_ori_wide", 
  "loading1_wide", "loading2_wide",
  "Q_p_res_c", "Q_p_Xperp_c",
  "Q_res_c", "Q_Xperp_c",
  "res_c_acf_mat", "res_c_pacf_mat",
  "se2_HAC", 
  "t_HAC", "t_OLS",
  "Xperp_c_acf_mat", "Xperp_c_pacf_mat",
  "Xperp_ori_wide",
  "Z", "Z.eta", "Z.gamm", 
  "bet", "eta", "gamm", "Lag",
  "loading1", "loading2",
  "n_state", "n_year",
  "northeast_states", "run", "shock.seq", "south_states",
  "state.unq", "statefip_ipt", "west_states", "year_ipt"
), 
file = "1121_c_AR1X&eps_epsShock_changeShock_L30.Rdata"
)

load("1121_c_AR1X&eps_epsShock_changeShock_L30.Rdata")

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
load("1121_c_AR1X&eps_epsShock_changeShock_L5.Rdata")
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





# use ggplot with given nominal level #---------------------------
qtl = qnorm(c(0.05, 0.025, 0.005), lower.tail = F) # non-coverage prob of 0.1, 0.05, 0.01
N = nrow(Z)
k = ncol(Z) + 1
t_HAC_HC1 = t_HAC * sqrt( (N - k) / N )
ncp = list() # This is a list of matrices. Each matrix corresponds to a qtl. Given qtl, each col represents each method. Each row is the empirical prob of |t_j| > qtl for some rho

for (loop.qtl in 1:length(qtl)){
  
  ncp.qtl = matrix(0, nrow = length(shock.seq), ncol = 2)
  
  colnames(ncp.qtl) = c("OLS", "HAC" )
  rownames(ncp.qtl) = shock.seq
  
  for (loop.s in 1:length(shock.seq)) {
    
    ncp.qtl[loop.s, "OLS"] = sum(abs(t_OLS[,loop.s]) > qtl[loop.qtl]) / run
    ncp.qtl[loop.s, "HAC"] = sum(abs(t_HAC_HC1[,loop.s]) > qtl[loop.qtl]) / run
    
  }
  
  ncp[[loop.qtl]] = ncp.qtl
  
  rm(ncp.qtl)
  
}
names(ncp) = c("nominal=0.1", "nominal=0.05", "nominal=0.01")

ncp
ncp.0.05 = data.frame(ncp = c(ncp$`nominal=0.05`[,"HAC"], ncp$`nominal=0.05`[,"OLS"]),
                      Shock = c(shock.seq, shock.seq),
                      Method = c(rep("HAC", times = length(shock.seq)),
                                 rep("OLS", times = length(shock.seq))
                      ),
                      nominal = "nominal level: 0.05")
library(ggplot2)
library("ggpubr")
library(ggh4x)

# ncp.LZ_X$number_of_grp = factor(ncp.LZ_X$number_of_grp)

ggplot(ncp.0.05, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_color_manual(values = c("red", "blue"))  +
  labs(x=expression(Shock~strength~of~epsilon), y = " ") +
  scale_y_continuous(limits = c(0.05, 0.12),
                     breaks = c(0.05, 0.07, 0.09, 0.11)) +
  theme_bw()+
  theme(legend.position = "none")

















































































































































































































































##############################################
### Use do parallel to calculate the code above (shock on eps)
### Note that seeds are different for different shock strengths
#############################################



rm(list = ls())

library(doParallel)
library(foreach)

detectCores() # 12 cores
# Start parallel backend
n_cores <- 10  # you have 16; this uses 8
cl <- makeCluster(n_cores)
registerDoParallel(cl)




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









# ----- start simulation --- #

run = 10000

Lag = 1



shock.seq = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

t_HAC = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC = matrix(0, nrow = run, ncol = length(shock.seq))



Z.eta = Z%*%eta
Z.gamm = Z%*%gamm


n_state = ncol(eps_perp_ori_wide)
n_year = nrow(eps_perp_ori_wide)

state.unq = unique(statefip_ipt)

sd.eps = sd(c(eps_perp_ori_wide))
sd.Xperp = sd(c(Xperp_ori_wide))

loading1 = rep( c(rep(1, times=45), rep(0, times=32)), times= n_state) #1938-1982, pre-internet
loading2 = rep( c(rep(0, times=45), rep(1, times=32)), times= n_state) #1983-2014, post-internet



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

set.seed(715)



for (loop.s in 1:length(shock.seq)) {
  
  loading.strength <- shock.seq[loop.s]
  
  # parallel over loop.run
  res_mat <- foreach(loop.run = 1:run,
                     .combine = rbind,
                     .packages = "stats",
                     .export = c("simulate_ar1", "leverages_chol")) %dopar% {
                       
                       ## ---- everything inside your inner loop, but without assigning to
                       ## t_HAC[...] / se2_HAC[...] and without print(loop.run) ---- ##
                       
                       # permute epsperp
                       eps_ar1_wide <- simulate_ar1(n_year = 77, n_state = 47,
                                                    sigma = sd.eps, phi = phi_eps)
                       eps_ar1 <- c(eps_ar1_wide)
                       
                       # add aggregate shock to generate eps
                       f1 <- rnorm(1)
                       f2 <- rnorm(1)
                       
                       eps <- sqrt(1 - loading.strength) * eps_ar1 +
                         sqrt(loading.strength) * sd.eps * (f1 * loading1 + f2 * loading2)
                       
                       # permute Xperp
                       Xperp_ar1_wide <- simulate_ar1(n_year = 77, n_state = 47,
                                                      sigma = sd.Xperp, phi = phi_Xperp)
                       Xperp <- c(Xperp_ar1_wide)
                       
                       X <- Z.eta + Xperp
                       
                       # generate y
                       y <- bet * X + Z.gamm + eps
                       
                       M <- cbind(X, Z)
                       colnames(M)[1] <- "X"
                       
                       fit.ols <- lm(y ~ M - 1)
                       
                       # --- HAC --- #
                       fit.X <- lm(X ~ Z - 1)
                       est.Xperp <- fit.X$residuals
                       est.eps.hat <- fit.ols$residuals
                       
                       Hardmard <- t(est.Xperp) %*%
                         ((est.eps.hat %*% t(est.eps.hat)) * Wgt_HAC) %*%
                         est.Xperp
                       
                       se2_denom <- (sum(est.Xperp * est.Xperp))^2
                       est.se2 <- Hardmard / se2_denom
                       
                       t_val <- (coef(summary(fit.ols))[1, "Estimate"] - bet) / sqrt(est.se2)
                       
                       # return one row: c(se2, t)
                       c(se2_HAC = as.numeric(est.se2),
                         t_HAC   = as.numeric(t_val))
                     }
  
  # res_mat is run x 2: column 1 = se2, column 2 = t
  se2_HAC[, loop.s] <- res_mat[, "se2_HAC"]
  t_HAC[,  loop.s]  <- res_mat[, "t_HAC"]
  
  cat("Finished shock index", loop.s, "with loading =", loading.strength, "\n")
}

stopCluster(cl)
registerDoSEQ()   # go back to sequential


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
  "Z",  
  "bet", "eta", "gamm", "Lag",
  "loading1", "loading2",
  "n_state", "n_year",
  "northeast_states", "run", "shock.seq", "south_states", "sd.eps", "sd.Xperp",
  "state.unq", "statefip_ipt", "west_states", "year_ipt"
), 
file = "1121_c_AR1X&eps_epsShock_changeShock_L1.Rdata"
)
































































































































































































































##############################################
### Use do parallel to calculate the code above (shock on Xperp)
### Note that seeds are different for different shock strengths
#############################################



rm(list = ls())

library(doParallel)
library(foreach)

detectCores() # 12 cores
# Start parallel backend
n_cores <- 8  # you have 16; this uses 8
cl <- makeCluster(n_cores)
registerDoParallel(cl)




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









# ----- start simulation --- #

run = 10000

Lag = 5



shock.seq = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

t_HAC = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC = matrix(0, nrow = run, ncol = length(shock.seq))



Z.eta = Z%*%eta
Z.gamm = Z%*%gamm


n_state = ncol(eps_perp_ori_wide)
n_year = nrow(eps_perp_ori_wide)

state.unq = unique(statefip_ipt)

sd.eps = sd(c(eps_perp_ori_wide))
sd.Xperp = sd(c(Xperp_ori_wide))

loading1 = rep( c(rep(1, times=45), rep(0, times=32)), times= n_state) #1938-1982, pre-internet
loading2 = rep( c(rep(0, times=45), rep(1, times=32)), times= n_state) #1983-2014, post-internet



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

set.seed(715)



for (loop.s in 1:length(shock.seq)) {
  
  loading.strength <- shock.seq[loop.s]
  
  # parallel over loop.run
  res_mat <- foreach(loop.run = 1:run,
                     .combine = rbind,
                     .packages = "stats",
                     .export = c("simulate_ar1", "leverages_chol")) %dopar% {
                       
                       ## ---- everything inside your inner loop, but without assigning to
                       ## t_HAC[...] / se2_HAC[...] and without print(loop.run) ---- ##
                       
                       # permute epsperp
                       eps_ar1_wide <- simulate_ar1(n_year = 77, n_state = 47,
                                                    sigma = sd.eps, phi = phi_eps)
                       eps <- c(eps_ar1_wide)
                       
                       
                       
                       # permute Xperp
                       Xperp_ar1_wide <- simulate_ar1(n_year = 77, n_state = 47,
                                                      sigma = sd.Xperp, phi = phi_Xperp)
                       Xperp_ar1 <- c(Xperp_ar1_wide)
                       
                       # add aggregate shock to generate eps
                       f1 <- rnorm(1)
                       f2 <- rnorm(1)
                       
                       Xperp <- sqrt(1 - loading.strength) * Xperp_ar1 +
                         sqrt(loading.strength) * sd.Xperp * (f1 * loading1 + f2 * loading2)
                       
                       X <- Z.eta + Xperp
                       
                       # generate y
                       y <- bet * X + Z.gamm + eps
                       
                       M <- cbind(X, Z)
                       colnames(M)[1] <- "X"
                       
                       fit.ols <- lm(y ~ M - 1)
                       
                       # --- HAC --- #
                       fit.X <- lm(X ~ Z - 1)
                       est.Xperp <- fit.X$residuals
                       est.eps.hat <- fit.ols$residuals
                       
                       Hardmard <- t(est.Xperp) %*%
                         ((est.eps.hat %*% t(est.eps.hat)) * Wgt_HAC) %*%
                         est.Xperp
                       
                       se2_denom <- (sum(est.Xperp * est.Xperp))^2
                       est.se2 <- Hardmard / se2_denom
                       
                       t_val <- (coef(summary(fit.ols))[1, "Estimate"] - bet) / sqrt(est.se2)
                       
                       # return one row: c(se2, t)
                       c(se2_HAC = as.numeric(est.se2),
                         t_HAC   = as.numeric(t_val))
                     }
  
  # res_mat is run x 2: column 1 = se2, column 2 = t
  se2_HAC[, loop.s] <- res_mat[, "se2_HAC"]
  t_HAC[,  loop.s]  <- res_mat[, "t_HAC"]
  
  cat("Finished shock index", loop.s, "with loading =", loading.strength, "\n")
}

stopCluster(cl)
registerDoSEQ()   # go back to sequential


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
  "Z",  
  "bet", "eta", "gamm", "Lag",
  "loading1", "loading2",
  "n_state", "n_year",
  "northeast_states", "run", "shock.seq", "south_states", "sd.eps", "sd.Xperp",
  "state.unq", "statefip_ipt", "west_states", "year_ipt"
), 
file = "1121_c_AR1X&eps_XShock_changeShock_L5.Rdata"
)
load("1121_c_AR1X&eps_XShock_changeShock_L5.Rdata")

# The shock on Xperp seems to perform better than shock on eps
# This is weird.























































































































































































































































































































































































##############################################
### Use do parallel to calculate the code above (shock on epsperp)
### Use the same seed for different shock strengths
#############################################


