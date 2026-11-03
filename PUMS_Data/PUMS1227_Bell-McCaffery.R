#########################################################################
#### Bell-McCaffrey adjustment with Satterthwaite degrees of freedom ####
#########################################################################

rm(list=ls())
library(dplyr)
library(clubSandwich)

# folder = "D:\\Research\\20251214_PUMS"
folder = Sys.getenv("PUMS_DATA_DIR", unset = "data")
load(file.path(folder, "1222_PUMA_setting.Rdata"))

rm(list = setdiff(ls(), c("dat_tgt")))
# save(list = c("dat_tgt"), file = "dat_tgt.Rdata")
# load("dat_tgt.Rdata")


# The function used in the calculation of Bell-McCaffery Aepshat and the Satterthwaite df 
apply_Ai_vec <- function(Mi, v, R, tol_eig = 1e-12, tol_one = 1e-12) {
  # Mi: (n_i x p), v: (n_i)
  # R: (p x p) such that S = (M'M)^{-1} = R %*% t(R)
  Ui <- Mi %*% R                       # n_i x p
  Ki <- crossprod(Ui)                  # p x p
  
  eg <- eigen(Ki, symmetric = TRUE)
  d  <- eg$values
  Q  <- eg$vectors
  
  keep <- which(d > tol_eig)
  if (length(keep) == 0L) return(v)
  
  d  <- d[keep]
  Q  <- Q[, keep, drop = FALSE]
  
  # Orthonormal basis for col(Ui): Bi = Ui Q diag(1/sqrt(d))
  Bi <- (Ui %*% Q) %*% diag(1 / sqrt(d), nrow = length(d))
  
  r  <- as.numeric(crossprod(Bi, v))   # length = rank
  # guard against 1-d <= 0 due to numerical issues
  denom <- pmax(1 - d, tol_one)
  mult  <- (1 / sqrt(denom)) - 1
  
  v + as.numeric(Bi %*% (mult * r))
}




# cluster must be a factor
# M is the whole design matrix
# epshat is the residual
# the output is Aepshat using the BRL from Bell and 
BRL_Aepshat <- function(M, epshat, cluster) {
  stopifnot(is.matrix(M), length(epshat) == nrow(M), length(cluster) == nrow(M))
  
  N <- nrow(M); p <- ncol(M)
  
  # S = (M'M)^{-1} via chol
  XtX  <- crossprod(M)
  Rchol <- chol(XtX)                   # upper: XtX = t(Rchol) %*% Rchol
  Rinv  <- backsolve(Rchol, diag(p))   # Rchol^{-1}
  S     <- tcrossprod(Rinv)            # (M'M)^{-1}
  R     <- t(chol(S))                  # lower so that S = R %*% t(R)
  
  # sort by cluster to stream contiguous blocks
  ord <- order(cluster)
  cl  <- cluster[ord]
  runs <- rle(as.integer(cl))
  
  Aeps_sorted <- numeric(N)
  pos <- 1L
  
  for (k in seq_along(runs$lengths)) {
    ni <- runs$lengths[k]
    idx_sorted <- pos:(pos + ni - 1L)
    idx <- ord[idx_sorted]
    
    Mi <- M[idx, , drop = FALSE]
    ei <- epshat[idx]
    
    Aeps_sorted[idx_sorted] <- apply_Ai_vec(Mi, ei, R)
    
    pos <- pos + ni
  }
  
  # unsort back to original order
  Aepshat <- numeric(N)
  Aepshat[ord] <- Aeps_sorted
  Aepshat
}


# The function to calculate the Satterthwaite df 
# cluster must be a factor!
# M is the whole design matrix
# the output is the Satterthwaite df
BRL_df_e1 <- function(M, cluster, tol_eig = 1e-12) {
  # cluster <- as.factor(cluster)
  N <- nrow(M); p <- ncol(M)
  
  XtX  <- crossprod(M)
  Rchol <- chol(XtX)
  Rinv  <- backsolve(Rchol, diag(p))
  S     <- tcrossprod(Rinv)           # (M'M)^{-1}
  R     <- t(chol(S))                 # S = R %*% t(R)
  
  e1 <- c(1, rep(0, p - 1))
  b  <- as.numeric(S %*% e1)
  
  ord <- order(cluster)
  cl  <- cluster[ord]
  runs <- rle(as.integer(cl))
  G <- length(runs$lengths)
  
  s <- numeric(G)
  Tmat <- matrix(0, nrow = G, ncol = p)
  
  pos <- 1L
  for (g in seq_len(G)) {
    ni <- runs$lengths[g]
    idx_sorted <- pos:(pos + ni - 1L)
    idx <- ord[idx_sorted]
    
    Mi <- M[idx, , drop = FALSE]
    
    # u_i = A_i (M_i b)
    v  <- as.numeric(Mi %*% b)
    ui <- apply_Ai_vec(Mi, v, R, tol_eig = tol_eig)
    
    s[g] <- sum(ui^2)
    Tmat[g, ] <- as.numeric(crossprod(Mi, ui))
    
    pos <- pos + ni
  }
  
  # C = T S T'
  A <- Tmat %*% S                 # (G x p)
  C <- A %*% t(Tmat)              # (G x G) with C_ij = t_i' S t_j
  
  Gmat <- -C
  diag(Gmat) <- s - diag(C)
  
  trG  <- sum(diag(Gmat))
  trG2 <- sum(Gmat^2)             # Frobenius norm^2 = tr(G^2) since symmetric
  
  f <- (trG^2) / trG2
  f
}
















# Function: draw one value per group, then broadcast it to all rows in that group
# Input: grp must be a factor indicating groups.
# mean and sd are the rnorm arguments
# Output: A vector with the same length as grp. The indices in the same group get the same value from N(mean, sd^2).
rnorm_by_group_factor <- function(grp, mean = 0, sd = 1) {
  u <- rnorm(nlevels(grp), mean, sd)
  u[as.integer(grp)]
}



run = 1000

DVSN_index = factor(dat_tgt$DIVISION)
DVSN_index

ST_index = factor(dat_tgt$STATE)
ST_index

PUMA_index = dat_tgt$PUMA_full
PUMA_index # 2462 levels



# Parameter       Point.Est  Profile 2.5% Profile 97.5%  Profile 16% Profile 84%  
# \sigma_{P,eps}    0.18360     0.1781696     0.1892735    0.1808121   0.1864441
# \sigma_{S,eps}    0.19543     0.1567992     0.2498636    0.1742891   0.2205696
# \sigma_{D,eps}    0.25298     0.1550931     0.4541844    0.1964173   0.3342215
# \sigma_{e,eps}    0.86189     0.8606802     0.8630962    0.8612742   0.8625001

# Parameter     Point.Est  Profile 2.5% Profile 97.5%  Profile 16%  Profile 84%
# \sigma_{P,X}   0.137140    0.13293210     0.1415388   0.13498173    0.1393471 
# \sigma_{S,X}   0.092509    0.07187226     0.1210694   0.08124700    0.1057896
# \sigma_{D,X}   0.078655    0.03767062     0.1503290   0.05670803    0.1081511
# \sigma_{e,X}   0.789353    0.78824815     0.7904609   0.78879219    0.7899149


sigm_XP = 0.13714
sigm_XS = 0.09251
# sigm_XD.seq = c(0, 0.05, 0.1, 0.15, 0.2) 
sigm_XD.seq = c(0, 0.038, 0.057, 0.079, 0.108, 0.15)
sigm_Xe = 0.78935

sigm_epsP = 0.1836
sigm_epsS = 0.1954
sigm_epsD = 0.2530
sigm_epse = 0.8619

N = nrow(dat_tgt)

one.N = as.matrix(rep(1, times=N))

Z = cbind(one.N, 
          dat_tgt[,c("HHL_Spn", "HHL_Other", "RMSP", "YRBLT1940", "YRBLT1950",
                     "YRBLT1960", "YRBLT1970", "YRBLT1980", "YRBLT1990", "YRBLT2000",
                     "YRBLT2010", "YRBLT2020", "YRBLT2021", "YRBLT2022", "YRBLT2023",
                     "YRBLT2024", "BLD1U", "BLDMU", "AGEP", "SEX",
                     "MSP_sgl", "RAC1P2", "RAC1P3", "RAC1P4", "RAC1P5",
                     "RAC1P6", "RAC1P7", "RAC1P8", "RAC1P9", "SCHL_HSCOMP", "SCHL_COCOMP",
                     "DIS" )])
Z = as.matrix(Z)
colnames(Z)[1] = "one"

colnames(Z)

ZtZ_inv = solve(t(Z)%*%Z)

rm(dat_tgt)

#       
#           (Intercept)      HHL_Spn    HHL_Other         RMSP      YRBLT1940    YRBLT1950  
eta =   c(  10.202178,     0.009120,     0.058568,     0.058592,    -0.003536,     0.041945,  
            #           YRBLT1960      YRBLT1970    YRBLT1980     YRBLT1990      YRBLT2000    YRBLT2010  
            0.042215,     0.032049,     0.084491,     0.121783,     0.160453,     0.203288,  
            #           YRBLT2020      YRBLT2021     YRBLT2022     YRBLT2023    YRBLT2024        BLD1U  
            0.189724,     0.213712,     0.210512,     0.226751,     0.195161,     0.329970,  
            #             BLDMU         AGEP            SEX         MSP_sgl       RAC1P2        RAC1P3  
            0.459882,    -0.008931,     0.292543,    -0.283887,    -0.153351,    -0.160468,  
            #             RAC1P4       RAC1P5          RAC1P6       RAC1P7        RAC1P8        RAC1P9  
            -0.286003,    -0.149710,     0.078426,     0.007157,    -0.059985,    -0.022738,  
            #           SCHL_HSCOMP  SCHL_COCOMP          DIS  
            0.301461,     0.839699,    -0.257458)

bet = 0.27604

#      (Intercept)     HHL_Spn     HHL_Other         RMSP       YRBLT1940    YRBLT1950  
gamm = c(7.237291,    0.124433,     0.182241,     0.057493,     0.069301,     0.143184,  
         #       YRBLT1960     YRBLT1970     YRBLT1980     YRBLT1990     YRBLT2000     YRBLT2010     YRBLT2020  
         0.161374,     0.171950,     0.254317,     0.305908,     0.361102,     0.422672,     0.482539,  
         #       YRBLT2021     YRBLT2022     YRBLT2023    YRBLT2024        BLD1U        BLDMU         AGEP  
         0.477168,     0.519450,     0.53289,     0.526381,     1.331883,     1.465468,     0.003567,  
         #         SEX           MSP_sgl       RAC1P2       RAC1P3        RAC1P4        RAC1P5         RAC1P6  
         -0.003794,    -0.042660,    -0.190255,    -0.252357,    -0.489927,    -0.084555,     0.316191,  
         #         RAC1P7        RAC1P8        RAC1P9      SCHL_HSCOMP  SCHL_COCOMP        DIS  
         0.351471,     0.110249,     0.056949,     0.088430,    0.337412,    -0.086665 )

Z.eta = Z%*%eta
Z.gamm = Z%*%gamm




t_ols = matrix(0, nrow = run, length(sigm_XD.seq))
se2_ols = matrix(0, nrow = run, length(sigm_XD.seq))



t_ehw = matrix(0, nrow = run, length(sigm_XD.seq))
se2_ehw = matrix(0, nrow = run, length(sigm_XD.seq))



t_LZ.P = matrix(0, nrow = run, length(sigm_XD.seq))
se2_LZ.P = matrix(0, nrow = run, length(sigm_XD.seq))
t_LZ.P_A = matrix(0, nrow = run, length(sigm_XD.seq))
se2_LZ.P_A = matrix(0, nrow = run, length(sigm_XD.seq))
df_LZ.P_A = matrix(0, nrow = run, length(sigm_XD.seq))



t_LZ.S = matrix(0, nrow = run, length(sigm_XD.seq))
se2_LZ.S = matrix(0, nrow = run, length(sigm_XD.seq))
t_LZ.S_A = matrix(0, nrow = run, length(sigm_XD.seq))
se2_LZ.S_A = matrix(0, nrow = run, length(sigm_XD.seq))
df_LZ.S_A = matrix(0, nrow = run, length(sigm_XD.seq))



t_LZ.D = matrix(0, nrow = run, length(sigm_XD.seq))
se2_LZ.D = matrix(0, nrow = run, length(sigm_XD.seq))
t_LZ.D_A = matrix(0, nrow = run, length(sigm_XD.seq))
se2_LZ.D_A = matrix(0, nrow = run, length(sigm_XD.seq))
df_LZ.D_A = matrix(0, nrow = run, length(sigm_XD.seq))




for (loop.s in 1:length(sigm_XD.seq)) {
  
  sigm_XD = sigm_XD.seq[loop.s]
  
  set.seed(1125)
  
  for (loop.run in 1:run) {
    
    # generate X #---------------------------------------
    U_XD = rnorm_by_group_factor(DVSN_index, mean = 0, sd = sigm_XD)
    U_XS = rnorm_by_group_factor(ST_index, mean = 0, sd = sigm_XS)
    U_XP = rnorm_by_group_factor(PUMA_index, mean = 0, sd = sigm_XP)
    U_Xe = rnorm(N, mean = 0, sd = sigm_Xe)
    Xperp.tilde = U_XD + U_XS + U_XP + U_Xe
    
    X = Z.eta + Xperp.tilde
    
    # generate eps #-------------------------------------
    U_epsD = rnorm_by_group_factor(DVSN_index, mean = 0, sd = sigm_epsD)
    U_epsS = rnorm_by_group_factor(ST_index, mean = 0, sd = sigm_epsS)
    U_epsP = rnorm_by_group_factor(PUMA_index, mean = 0, sd = sigm_epsP)
    U_epse = rnorm(N, mean = 0, sd = sigm_epse)
    
    eps = U_epsD + U_epsS + U_epsP + U_epse 
    # plot(eps)
    
    # generate y #----------------------------------------
    y = bet * X +  Z.gamm + eps
    
    fit.X = lm(X~Z-1)
    Xperp = fit.X$residuals
    
    M = cbind(X, Z) 
    colnames(M)[1] = "X"
    
    fit.Y=lm(y~M-1)
    epshat = fit.Y$residuals
    
    bethat = summary(fit.Y)$coefficients[1, "Estimate"]
    
    # ols #-----------------------------------------------
    # t_ols[loop.run,loop.s] = ( sqrt(N) * sum(Xperp * eps) ) / ( sqrt(sum(Xperp^2)) * sqrt(sum(epshat^2)) )
    se2_ols[loop.run, loop.s] = (summary(fit.Y)$coefficients[1, "Std. Error"])^2
    t_ols[loop.run,loop.s] = ( bethat - bet ) / sqrt(se2_ols[loop.run, loop.s])
    #-----------------------------------------------------
    
    Xperp_norm2 = sum(Xperp^2)
    
    # EHW #-----------------------------------------------
    se2_ehw[loop.run, loop.s] = sum( Xperp*epshat*Xperp*epshat ) / (Xperp_norm2^2)
    t_ehw[loop.run,loop.s] = ( bethat - bet )/ sqrt(se2_ehw[loop.run, loop.s])
    #-----------------------------------------------------
    
    
    
    # LZ using PUMA groups #-------------------------------------------------------
    XGeps_P = rowsum(Xperp*epshat, PUMA_index, reorder = FALSE)
    se2_LZ.P[loop.run, loop.s] = sum(XGeps_P^2) / (Xperp_norm2^2)
    t_LZ.P[loop.run,loop.s] = ( bethat - bet ) / sqrt(se2_LZ.P[loop.run, loop.s])
    
    epshat_P_A = BRL_Aepshat(M=M, epshat = epshat, cluster = PUMA_index)   
    XGeps_P_A = rowsum(Xperp*epshat_P_A, PUMA_index, reorder = FALSE)
    se2_LZ.P_A[loop.run, loop.s] = sum(XGeps_P_A^2) / (Xperp_norm2^2)
    t_LZ.P_A[loop.run,loop.s] = ( bethat - bet ) / sqrt(se2_LZ.P_A[loop.run, loop.s])
    df_LZ.P_A[loop.run, loop.s] = BRL_df_e1(M=M, cluster = PUMA_index)
    #---------------------------------------------------------------------------
    
    
    
    
    # LZ using STATE groups #-------------------------------------------------------
    XGeps_S = rowsum(Xperp*epshat, ST_index, reorder = FALSE)
    se2_LZ.S[loop.run, loop.s] = sum(XGeps_S^2) / (Xperp_norm2^2)
    t_LZ.S[loop.run,loop.s] = ( bethat - bet ) / sqrt(se2_LZ.S[loop.run, loop.s])
    
    epshat_S_A = BRL_Aepshat(M=M, epshat = epshat, cluster = ST_index)
    XGeps_S_A = rowsum(Xperp*epshat_S_A, ST_index, reorder = FALSE)
    se2_LZ.S_A[loop.run, loop.s] = sum(XGeps_S_A^2) / (Xperp_norm2^2)
    t_LZ.S_A[loop.run,loop.s] = ( bethat - bet ) / sqrt(se2_LZ.S_A[loop.run, loop.s])
    df_LZ.S_A[loop.run, loop.s] = BRL_df_e1(M=M, cluster = ST_index)
    #---------------------------------------------------------------------------
    
    
    # LZ using DIVISION groups #-------------------------------------------------------
    XGeps_D = rowsum(Xperp*epshat, DVSN_index, reorder = FALSE)
    se2_LZ.D[loop.run, loop.s] = sum(XGeps_D^2) / (Xperp_norm2^2)
    t_LZ.D[loop.run,loop.s] = ( bethat - bet ) / sqrt(se2_LZ.D[loop.run, loop.s])
    
    epshat_D_A = BRL_Aepshat(M=M, epshat = epshat, cluster = DVSN_index)
    XGeps_D_A = rowsum(Xperp*epshat_D_A, DVSN_index, reorder = FALSE)
    se2_LZ.D_A[loop.run, loop.s] = sum(XGeps_D_A^2) / (Xperp_norm2^2)
    t_LZ.D_A[loop.run,loop.s] = ( bethat - bet ) / sqrt(se2_LZ.D_A[loop.run, loop.s])
    df_LZ.D_A[loop.run, loop.s] = BRL_df_e1(M=M, cluster = DVSN_index)
    #---------------------------------------------------------------------------
    
    
    print(loop.run)
    
  }
  
  print(loop.s)
  
}


save(list = c("df_LZ.D_A", "df_LZ.P_A", "df_LZ.S_A",
              "se2_ehw", "se2_LZ.D", "se2_LZ.D_A",
              "se2_LZ.P", "se2_LZ.P_A", "se2_LZ.S", "se2_LZ.S_A", "se2_ols",
              "t_ehw", "t_LZ.D", "t_LZ.D_A", "t_LZ.P", "t_LZ.P_A", "t_LZ.S", "t_LZ.S_A", "t_ols", 
              "Z",
              "bet", "DVSN_index", "eta", "gamm", "PUMA_index", "run", 
              "sigm_epsD", "sigm_epse", "sigm_epsP", "sigm_epsS",
              "sigm_XD.seq", "sigm_Xe", "sigm_XP", "sigm_XS",
              "ST_index"),
     file = "0103PUMS_change_sigm_XD_Bell_df.Rdata")
































































boxplot(cbind(se2_LZ.D, se2_LZ.D_A))
boxplot(cbind(se2_LZ.P, se2_LZ.P_A))
boxplot(cbind(se2_LZ.S, se2_LZ.S_A))

loop.s=4
se2_plot = cbind(se2_ehw[,loop.s], se2_LZ.P[,loop.s], se2_LZ.S[,loop.s], se2_LZ.D[,loop.s])
se2_plot = cbind(se2_ehw[,loop.s], se2_LZ.P_A[,loop.s], se2_LZ.S_A[,loop.s], se2_LZ.D_A[,loop.s])
colnames(se2_plot) = c("EHW", "PUMA", "STATE", "DIVISION")
boxplot(sqrt(se2_plot),
        main = expression(SE~when~sigma[XD]==0.079) )



plot(x = sqrt(se2_LZ.S)[,2], y = sqrt(se2_LZ.S_A)[,2])
abline(a=0,b=1,col=2,lwd=2)

plot(x = sqrt(se2_LZ.P)[,1], y = sqrt(se2_LZ.P_A)[,1])
abline(a=0,b=1,col=2,lwd=2)
max(abs(se2_LZ.P - se2_LZ.P_A)) # 3.63425e-08

plot(x = sqrt(se2_LZ.D)[,1], y = sqrt(se2_LZ.D_A)[,1])
abline(a=0,b=1,col=2,lwd=2)

boxplot(se2_LZ.S_A[,6] / se2_LZ.S[,6])

plot(x=abs(t_LZ.S[,1]), y = abs(t_LZ.S_A[,1]))
abline(a=0,b=1,col=2,lwd=2)


bias_ols = sqrt(se2_ols) * t_ols
bias_ehw = sqrt(se2_ehw) * t_ehw

bias_LZ.S = sqrt(se2_LZ.S) * t_LZ.S
bias_LZ.S_A = sqrt(se2_LZ.S_A) * t_LZ.S_A

bias_LZ.D = sqrt(se2_LZ.D) * t_LZ.D
bias_LZ.D_A = sqrt(se2_LZ.D_A) * t_LZ.D_A

bias_LZ.P = sqrt(se2_LZ.P) * t_LZ.P
bias_LZ.P_A = sqrt(se2_LZ.P_A) * t_LZ.P_A

max(abs(bias_ols - bias_ehw))
max(abs(bias_ols - bias_LZ.S))
max(abs(bias_ols - bias_LZ.S_A))
max(abs(bias_ols - bias_LZ.D))
max(abs(bias_ols - bias_LZ.D_A))
max(abs(bias_ols - bias_LZ.P))
max(abs(bias_ols - bias_LZ.P_A))

max(abs(t_LZ.S - bias_ols / sqrt(se2_LZ.S)))
t_LZ.S_A - bias_ols / sqrt(se2_LZ.S_A)
t_LZ.D - bias_ols / sqrt(se2_LZ.D)
t_LZ.D_A - bias_ols / sqrt(se2_LZ.D_A)
t_LZ.P - bias_ols / sqrt(se2_LZ.P)
t_LZ.P_A - bias_ols / sqrt(se2_LZ.P_A)

# draw qqplots with plot #------------------------------------------------------
prob = seq(0.01,0.99, by=0.01)

loop.s=1

qt.emp = quantile(t_ehw[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("EHW when ", " sigma_XD =", sigm_XD.seq[loop.s]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_ols[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("OLS when ", " sigma_XD =", sigm_XD.seq[loop.s]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.D[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("t_LZ.D when ", " sigma_XD = ", sigm_XD.seq[loop.s]))
abline(a=0, b=1, col=2, lwd=3)


qt.emp = quantile(t_LZ.S[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("t_LZ.S when ", " sigma_XD =", sigm_XD.seq[loop.s]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.P[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("t_LZ.P when ", " sigma_XD =", sigm_XD.seq[loop.s]))
abline(a=0, b=1, col=2, lwd=3)



qt.emp = quantile(t_LZ.D_A[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("t_LZ.D_A when ", " sigma_XD = ", sigm_XD.seq[loop.s]))
abline(a=0, b=1, col=2, lwd=3)


qt.emp = quantile(t_LZ.S_A[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("t_LZ.S_A when ", " sigma_XD =", sigm_XD.seq[loop.s]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.P_A[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("t_LZ.P_A when ", " sigma_XD =", sigm_XD.seq[loop.s]))
abline(a=0, b=1, col=2, lwd=3)



# Compare the qqplots of Liang-Zeger SE's based on three levels when \sigma_XD = 0:
loop.s=1
par(mfrow = c(2, 3))

#--- Not adjusted --- #
qt.emp = quantile(t_LZ.D[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("DIVISION, ", " sigma_XD = ", sigm_XD.seq[loop.s]),
     xlab = "Empirical", ylab = "Theoretical")
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.S[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("STATE, ", " sigma_XD =", sigm_XD.seq[loop.s]),
     xlab = "Empirical", ylab = "Theoretical")
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.P[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("PUMA, ", " sigma_XD =", sigm_XD.seq[loop.s]),
     xlab = "Empirical", ylab = "Theoretical")
abline(a=0, b=1, col=2, lwd=3)

# ----- adjusted ------ #
qt.emp = quantile(t_LZ.D_A[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("Adjusted DIVISION, ", " sigma_XD = ", sigm_XD.seq[loop.s]),
     xlab = "Empirical", ylab = "Theoretical")
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.S_A[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("Adjusted STATE, ", " sigma_XD =", sigm_XD.seq[loop.s]),
     xlab = "Empirical", ylab = "Theoretical")
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.P_A[,loop.s], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("Adjusted PUMA, ", " sigma_XD =", sigm_XD.seq[loop.s]),
     xlab = "Empirical", ylab = "Theoretical")
abline(a=0, b=1, col=2, lwd=3)
# ------------------------------------- #
par(mfrow = c(1, 1))












# make plot about non-coverage prob #------------------
nominal.seq = c(0.1, 0.05, 0.01)

ncp = list() # This is a list of matrices. 

for (loop.nominal in 1:length(nominal.seq)){
  
  nominal_level = nominal.seq[loop.nominal]
  
  ncp_nominal = matrix(0, nrow = length(sigm_XD.seq), ncol = 8)
  
  colnames(ncp_nominal) = c("OLS", "EHW", "LZ_P", "LZ_S", "LZ_D", "LZ_P_A", "LZ_S_A", "LZ_D_A")
  rownames(ncp_nominal) = sigm_XD.seq
  
  for (loop.s in 1:length(sigm_XD.seq)) {
    
    ncp_nominal[loop.s, "OLS"] = sum(abs(t_ols[,loop.s]) > qnorm( (nominal_level/2), lower.tail = F ) ) / run
    ncp_nominal[loop.s, "EHW"] = sum(abs(t_ehw[,loop.s]) > qnorm( (nominal_level/2), lower.tail = F ) ) / run
    ncp_nominal[loop.s, "LZ_P"] = sum(abs(t_LZ.P[,loop.s]) > qnorm( (nominal_level/2), lower.tail = F ) ) / run
    ncp_nominal[loop.s, "LZ_S"] = sum(abs(t_LZ.S[,loop.s]) > qnorm( (nominal_level/2), lower.tail = F ) ) / run
    ncp_nominal[loop.s, "LZ_D"] = sum(abs(t_LZ.D[,loop.s]) > qnorm( (nominal_level/2), lower.tail = F ) ) / run
    ncp_nominal[loop.s, "LZ_P_A"] = sum( abs(t_LZ.P_A[,loop.s]) >  qt( (nominal_level/2), df=df_LZ.P_A[,loop.s], lower.tail = F ) ) / run
    ncp_nominal[loop.s, "LZ_S_A"] = sum( abs(t_LZ.S_A[,loop.s]) >  qt( (nominal_level/2), df=df_LZ.S_A[,loop.s], lower.tail = F ) ) / run
    ncp_nominal[loop.s, "LZ_D_A"] = sum( abs(t_LZ.D_A[,loop.s]) >  qt( (nominal_level/2), df=df_LZ.D_A[,loop.s], lower.tail = F ) ) / run
    
  }
  
  ncp[[loop.nominal]] = ncp_nominal
  
  rm(ncp_nominal)
  
}
names(ncp) = c("nominal=0.1", "nominal=0.05", "nominal=0.01")

ncp



library(ggplot2)
library("ggpubr")
library(ggh4x)

## Adjusted LZ
ncp.0.05 = data.frame(ncp = c(ncp$`nominal=0.05`[,"LZ_P_A"],
                              ncp$`nominal=0.05`[,"LZ_D_A"],
                              ncp$`nominal=0.05`[,"LZ_S_A"]),
                      Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                      Method = c(rep("PUMA", times = length(sigm_XD.seq) ),
                                 rep("DIVISION", times = length(sigm_XD.seq) ),
                                 rep("STATE", times = length(sigm_XD.seq) )
                      ),
                      nominal = "nominal level: 0.05")

ncp.0.05$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.05, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_y_continuous(limits = c(-0.02, 1.03 ),
                     breaks = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1),
                     expand = expansion( mult = c(0, 0) ) ) +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) + 
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.05$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  +
  theme(
    ## Legend inside plot (bottom-right)
    legend.position = c(1, 0.4),
    legend.justification = c("right", "bottom"),
    #    legend.background = element_rect(fill = "white", color = "black"),
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    
    ## Axis labels
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    
    ## Axis tick labels
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    
    ## Facet strip text
    strip.text = element_text(size = 16)
  ) # + theme(legend.position = "none")



ncp.0.01 = data.frame(ncp = c(ncp$`nominal=0.01`[,"LZ_P_A"],
                              ncp$`nominal=0.01`[,"LZ_D_A"],
                              ncp$`nominal=0.01`[,"LZ_S_A"]),
                      Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                      Method = c(rep("PUMA", times = length(sigm_XD.seq)  ),
                                 rep("DIVISION", times = length(sigm_XD.seq)  ),
                                 rep("STATE", times = length(sigm_XD.seq)  )
                      ),
                      nominal = "nominal level: 0.01")

ncp.0.01$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.01, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.01, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_y_continuous(limits = c(-0.02, max(ncp$`nominal=0.01`[,c("LZ_P",  "LZ_S", "LZ_D")]) + 0.02),
                     breaks = c(0.01, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) +
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.01$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  # + theme(legend.position = "none")



ncp.0.1 = data.frame(ncp = c(ncp$`nominal=0.1`[,"LZ_P_A"],
                             ncp$`nominal=0.1`[,"LZ_D_A"],
                             ncp$`nominal=0.1`[,"LZ_S_A"]),
                     Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                     Method = c(rep("PUMA", times = length(sigm_XD.seq) ),
                                rep("DIVISION", times = length(sigm_XD.seq) ),
                                rep("STATE", times = length(sigm_XD.seq) )
                     ),
                     nominal = "nominal level: 0.1")

ncp.0.1$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.1, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.1, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) +
  scale_y_continuous(limits = c(-0.02, max(ncp$`nominal=0.1`[,c("LZ_P",  "LZ_S", "LZ_D")]) + 0.02),
                     breaks = c(0, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.1$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  # + theme(legend.position = "none")








## Not adjusted LZ
ncp.0.05 = data.frame(ncp = c(ncp$`nominal=0.05`[,"LZ_P"],
                              ncp$`nominal=0.05`[,"LZ_D"],
                              ncp$`nominal=0.05`[,"LZ_S"]),
                      Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                      Method = c(rep("PUMA", times = length(sigm_XD.seq) ),
                                 rep("DIVISION", times = length(sigm_XD.seq) ),
                                 rep("STATE", times = length(sigm_XD.seq) )
                      ),
                      nominal = "nominal level: 0.05")

ncp.0.05$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.05, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_y_continuous(limits = c(-0.02, (max(ncp$`nominal=0.05`[,c("LZ_P",  "LZ_S", "LZ_D")])+0.02) ),
                     breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) + 
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.05$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  # + theme(legend.position = "none")


ncp.0.01 = data.frame(ncp = c(ncp$`nominal=0.01`[,"LZ_P"],
                              ncp$`nominal=0.01`[,"LZ_D"],
                              ncp$`nominal=0.01`[,"LZ_S"]),
                      Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                      Method = c(rep("PUMA", times = length(sigm_XD.seq)  ),
                                 rep("DIVISION", times = length(sigm_XD.seq)  ),
                                 rep("STATE", times = length(sigm_XD.seq)  )
                      ),
                      nominal = "nominal level: 0.01")

ncp.0.01$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.01, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.01, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_y_continuous(limits = c(-0.02, max(ncp$`nominal=0.01`[,c("LZ_P",  "LZ_S", "LZ_D")]) + 0.02),
                     breaks = c(0.01, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) +
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.01$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  # + theme(legend.position = "none")


ncp.0.1 = data.frame(ncp = c(ncp$`nominal=0.1`[,"LZ_P"],
                             ncp$`nominal=0.1`[,"LZ_D"],
                             ncp$`nominal=0.1`[,"LZ_S"]),
                     Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                     Method = c(rep("PUMA", times = length(sigm_XD.seq) ),
                                rep("DIVISION", times = length(sigm_XD.seq) ),
                                rep("STATE", times = length(sigm_XD.seq) )
                     ),
                     nominal = "nominal level: 0.1")

ncp.0.1$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.1, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.1, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) +
  scale_y_continuous(limits = c(-0.02, max(ncp$`nominal=0.1`[,c("LZ_P",  "LZ_S", "LZ_D")]) + 0.02),
                     breaks = c(0, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.1$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  # + theme(legend.position = "none")












## The critical value of BM adjustment uses t distribution with complicated df
## I use number of cluster - 1 first
rm(list=ls())
load("0103PUMS_change_sigm_XD_Bell_df.Rdata")
# make plot about non-coverage prob #------------------
ncp = list() # This is a list of matrices. Each matrix corresponds to a qtl. Given qtl, each col represents each method. Each row is the empirical prob of |t_j| > qtl for some rho
nominal.seq = c(0.1, 0.05, 0.01)
for (loop.nominal in 1:length(nominal.seq)){
  
  nominal_level = nominal.seq[loop.nominal]
  
  ncp.loop = matrix(0, nrow = length(sigm_XD.seq), ncol = 3)
  
  colnames(ncp.loop) = c("LZ_P_A", "LZ_S_A", "LZ_D_A")
  rownames(ncp.loop) = sigm_XD.seq
  
  for (loop.s in 1:length(sigm_XD.seq)) {
    
    ncp.loop[loop.s, "LZ_P_A"] = sum(abs(t_LZ.P_A[,loop.s]) > qt( (1-(nominal_level/2)), df = (nlevels(PUMA_index) - 1) ) ) / run
    ncp.loop[loop.s, "LZ_S_A"] = sum(abs(t_LZ.S_A[,loop.s]) > qt( (1-(nominal_level/2)), df = (nlevels(ST_index) - 1) ) ) / run
    ncp.loop[loop.s, "LZ_D_A"] = sum(abs(t_LZ.D_A[,loop.s]) > qt( (1-(nominal_level/2)), df = (nlevels(DVSN_index) - 1) ) ) / run
    
  }
  
  ncp[[loop.nominal]] = ncp.loop
  
  rm(ncp.loop)
  
}
names(ncp) = c("nominal=0.1", "nominal=0.05", "nominal=0.01")

ncp

library(ggplot2)
library("ggpubr")
library(ggh4x)

## Adjusted LZ
ncp.0.05 = data.frame(ncp = c(ncp$`nominal=0.05`[,"LZ_P_A"],
                              ncp$`nominal=0.05`[,"LZ_D_A"],
                              ncp$`nominal=0.05`[,"LZ_S_A"]),
                      Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                      Method = c(rep("PUMA", times = length(sigm_XD.seq) ),
                                 rep("DIVISION", times = length(sigm_XD.seq) ),
                                 rep("STATE", times = length(sigm_XD.seq) )
                      ),
                      nominal = "nominal level: 0.05")

ncp.0.05$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.05, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_y_continuous(limits = c(-0.02, (max(ncp$`nominal=0.05`)+0.02) ),
                     breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) + 
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.05$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  # + theme(legend.position = "none")



ncp.0.01 = data.frame(ncp = c(ncp$`nominal=0.01`[,"LZ_P_A"],
                              ncp$`nominal=0.01`[,"LZ_D_A"],
                              ncp$`nominal=0.01`[,"LZ_S_A"]),
                      Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                      Method = c(rep("PUMA", times = length(sigm_XD.seq)  ),
                                 rep("DIVISION", times = length(sigm_XD.seq)  ),
                                 rep("STATE", times = length(sigm_XD.seq)  )
                      ),
                      nominal = "nominal level: 0.01")

ncp.0.01$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.01, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.01, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_y_continuous(limits = c(-0.02, max(ncp$`nominal=0.01`) + 0.02),
                     breaks = c(0.01, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) +
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.01$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  # + theme(legend.position = "none")



ncp.0.1 = data.frame(ncp = c(ncp$`nominal=0.1`[,"LZ_P_A"],
                             ncp$`nominal=0.1`[,"LZ_D_A"],
                             ncp$`nominal=0.1`[,"LZ_S_A"]),
                     Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
                     Method = c(rep("PUMA", times = length(sigm_XD.seq) ),
                                rep("DIVISION", times = length(sigm_XD.seq) ),
                                rep("STATE", times = length(sigm_XD.seq) )
                     ),
                     nominal = "nominal level: 0.1")

ncp.0.1$Method = factor(ncp.0.05$Method, levels = c("DIVISION", "STATE", "PUMA"))

ggplot(ncp.0.1, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.1, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) +
  scale_y_continuous(limits = c(-0.02, max(ncp$`nominal=0.1`) + 0.02),
                     breaks = c(0, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.1$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  # + theme(legend.position = "none")

# It seems that using nubmer of cluster -1 as the df of t-distribution, the division stays stable at the nominal level of ncp, which is best. State only achieve the nominal at \sigma_XD = 0.


# ---------- draw the plot of 0.05 level ncp -------------- #
ggplot(ncp.0.05, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_y_continuous(limits = c(-0.02, (max(ncp$`nominal=0.05`)+0.02) ),
                     breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion( mult = c(0, 0) ) ) +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) + 
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(ncp.0.05$nominal)),
    aes(x = sigm_XD_real, y = ncp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  +
  theme(
    ## Legend inside plot (bottom-right)
    legend.position = c(1, 0.5),
    legend.justification = c("right", "bottom"),
    #    legend.background = element_rect(fill = "white", color = "black"),
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    
    ## Axis labels
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    
    ## Axis tick labels
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    
    ## Facet strip text
    strip.text = element_text(size = 16)
  ) # + theme(legend.position = "none")

# ------ coverage prob plot ------------ #
cp.frame = ncp.0.05
cp.frame$cp = 1 - ncp.0.05$ncp
cp.frame$nominal = "nominal level: 0.95"

ggplot(cp.frame, aes(x = Shock, y = cp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.95, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ") +
  scale_y_continuous(limits = c(0.1, 1.03 ),
                     breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1),
                     expand = expansion( mult = c(0, 0) ) ) +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) + 
  theme_bw() +
  geom_point(  # add a point to indicate the point estimate of \sigm_XD from the real data
    data = data.frame(sigm_XD_real = 0.079, cp = 0.1, nominal = unique(cp.frame$nominal)),
    aes(x = sigm_XD_real, y = cp),
    inherit.aes = FALSE,
    color = "red",
    shape = 17,   # triangle
    size = 3
  )  +
  theme(
    ## Legend inside plot (bottom-right)
    legend.position = c(0.95, 0.4),
    legend.justification = c("right", "bottom"),
    #    legend.background = element_rect(fill = "white", color = "black"),
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    
    ## Axis labels
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    
    ## Axis tick labels
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    
    ## Facet strip text
    strip.text = element_text(size = 16)
  ) # + theme(legend.position = "none")




































































