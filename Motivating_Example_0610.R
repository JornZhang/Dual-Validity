#########################################
###  EXP1 in Motivation of the Intro ####
#####  For the validity of HAC ##########
#########################################

rm(list = ls())


library("sandwich")
library("lmtest")
library(Matrix)

# The function to get V0.5 given V.#--------------------------------------------
# The input V must be a NND matrix.#
Vhalf = function(V){
  aa = eigen(V)
  #Lambda=diag(c(aa$values))
  Q=aa$vectors
  Lambda0.5=diag(sqrt(aa$values))
  V0.5=Q%*%Lambda0.5%*%t(Q)
  return(list(V0.5 = V0.5, error = max(abs(V0.5%*%V0.5-V)) ))
}

# If n is too large, eigen() will produce error message #-----------------------
# if V is NND, the SVD is the same as eigen #
Vhalf = function(V){
  aa = svd(V)
  Lambda=diag(c(aa$d)) # max(abs((aa$u - aa$v)) is small because aa$u and aa$v are both eigenvectors of V, if V is NND.
  Q=aa$u  # max(abs(Q%*%Lambda%*%t(Q) - V)) is small, so Q%*%Lambda%*%t(Q) is just the eigen decomposition of V
  Lambda0.5=diag(sqrt(aa$d))
  V0.5=Q%*%Lambda0.5%*%t(Q)
  return(list(V0.5 = V0.5, error = max(abs(V0.5%*%V0.5-V)) ))
}
#-----------------------------------------------------------------------------


rho.seq = seq(-0.9, 0.9, by=0.15)

n=10000

d = 1

beta = 1

gamma = 1 # the first one is the target in the outcome model

run = 1000


t_ols = matrix(0, nrow = run, ncol = length(rho.seq))

t_ehw = matrix(0, nrow = run, ncol = length(rho.seq))

t_hac = matrix(0, nrow = run, ncol = length(rho.seq))

one.n = as.matrix(rep(1, times=n))

# since eps has two factors, loading_eps is the loading matrix for these two factors
loading_eps = cbind( c(rep(1, times = (n/2)), rep(0, times = (n/2))),
                     c(rep(0, times = (n/2)), rep(1, times = (n/2))) )


# HAC weight matrix #-----------------------------------------------
Lag = floor( 100 * ( (n/100)^(2/9) ) )
Wgt = matrix(0, n, n)
# Newey-West with finite sample adjustment
for (Wgtrow in 1:n) {
  for (Wgtcol in 1:n) {
    Wgt[Wgtrow, Wgtcol] = max( 1 - abs(Wgtrow-Wgtcol)/(Lag+1), 0  ) 
  }
}
rm(Wgtrow, Wgtcol)
#--------------------------------------------------


# start the simulation #---------------------------------------------
set.seed(715)
for (loop.rho in 1:length(rho.seq)) {
  
  rho = rho.seq[loop.rho]
  
  # V = matrix(0, n, n)
  
  # if (rho==0) {
  #   V = diag(n)
  # } else {
  #   for (Vrow in 1:n) {
  #     for (Vcol in 1:n) {
  #       V[Vrow, Vcol] = rho^(abs(Vrow-Vcol))
  #     }
  #   }
  #  rm(Vrow, Vcol)
  # }
  
  # tmp = Vhalf(V)
  # V0.5 = tmp$V0.5
  # print(tmp$error)
  # rm(tmp)
  
  for (loop.run in 1:run) {
    
    
    
    # generate Z #---------------------------------------
    # this Z has intercept
    Z = one.n
    
    # generate X #---------------------------------------
    # w_X =  as.matrix(rnorm(n))
    # X = V0.5%*%w_X # X has AR(1) structure
    
    # another way to generate X #------------------------
    X = as.matrix(arima.sim(n=n, 
                            model = list(ar = rho),
                            sd    = 1     # standard deviation of white noise
    ) # the default distribution of innovation is normal
    )
    
    M = cbind(X, Z) 
    
    # generate eps #-------------------------------------
    f_eps = rnorm(2)
    # w_eps = runif(n, min = -sqrt(3), max = sqrt(3))
    # w_eps = 2*rbinom(n,1,0.5)-1
    # w_eps = rgamma(n, shape = 2, scale = (1/sqrt(2)) ) - sqrt(2)
    # w_eps = rgamma(n, shape = 1, scale = 1) - 1
    w_eps =  rnorm(n)
    eps = sqrt(0.9) * loading_eps%*%f_eps + sqrt(0.1) * w_eps
    
    
    # generate y #----------------------------------------
    y = beta * X +  gamma * Z + eps
    
    
    # ols #-----------------------------------------------
    fit.ols=lm(y~M-1)
    t_ols[loop.run,loop.rho] = ( summary(fit.ols)$coefficients[1, "Estimate"] - beta ) / summary(fit.ols)$coefficients[1, "Std. Error"]
    #-----------------------------------------------------
    
    
    
    # EHW #-----------------------------------------------
    vcov.ehw=solve(t(M)%*%M) %*% (t(M) %*% diag(fit.ols$residuals)%*%diag(fit.ols$residuals) %*% M) %*% solve(t(M)%*%M)
    t_ehw[loop.run,loop.rho] = ( summary(fit.ols)$coefficients[1, "Estimate"] - beta ) / sqrt(vcov.ehw[1,1])
    #-----------------------------------------------------
    
    
    # HAC #-----------------------------------------------
    vcov.hac = solve(t(M)%*%M)%*%(t(M)%*%((fit.ols$residuals%*%t(fit.ols$residuals))*Wgt)%*%M)%*%solve(t(M)%*%M)
    t_hac[loop.run,loop.rho] = ( summary(fit.ols)$coefficients[1, "Estimate"] - beta ) / sqrt(vcov.hac[1,1])
    #-------------------------------------------------------
    
    
    print(loop.run)
    
  }
  
  print(loop.rho)
  
}


save(list = c("loading_eps", "t_ehw", "t_hac", "t_ols", "Lag",
              "beta", "d", "gamma", "n", "rho.seq", "run"),
     file = "HAC_X_vary_rho_n=5000.Rdata")



par(mfrow = c(1, 1))
plot(eps, type="b")
plot(X, type="b")

# draw qqplots with plot #------------------------------------------------------
prob = seq(0.01,0.99, by=0.01)

loop.rho=13

qt.emp = quantile(t_ehw[,loop.rho], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("EHW", " rho=", rho.seq[loop.rho]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_ols[,loop.rho], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("OLS", " rho=", rho.seq[loop.rho]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_hac[,loop.rho], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("HAC", " rho=", rho.seq[loop.rho]))
abline(a=0, b=1, col=2, lwd=3)




# make plot about non-coverage prob #------------------
qtl = qnorm(c(0.05, 0.025, 0.005), lower.tail = F) # non-coverage prob of 0.1, 0.05, 0.01

ncp = list() # This is a list of matrices. Each matrix corresponds to a qtl. Given qtl, each col represents each method. Each row is the empirical prob of |t_j| > qtl for some rho

for (loop.qtl in 1:length(qtl)){
  
  ncp.qtl = matrix(0, nrow = length(rho.seq), ncol = 3)
  
  colnames(ncp.qtl) = c("HAC", "OLS", "EHW")
  rownames(ncp.qtl) = rho.seq
  
  for (loop.rho in 1:length(rho.seq)) {
    
    ncp.qtl[loop.rho, "HAC"] = sum(abs(t_hac[,loop.rho]) > qtl[loop.qtl]) / run
    ncp.qtl[loop.rho, "OLS"] = sum(abs(t_ols[,loop.rho]) > qtl[loop.qtl]) / run
    ncp.qtl[loop.rho, "EHW"] = sum(abs(t_ehw[,loop.rho]) > qtl[loop.qtl]) / run
    
  }
  
  ncp[[loop.qtl]] = ncp.qtl
  
  rm(ncp.qtl)
  
}
names(ncp) = c("nominal=0.1", "nominal=0.05", "nominal=0.01")

plot(x=rho.seq, y=ncp$`nominal=0.1`[,"HAC"], 
     type = "b", lwd=2, pch=19,
     col=2,
     ylim = c(min(ncp$`nominal=0.1`), max(ncp$`nominal=0.1`)),
     xlab = expression(phi),
     ylab = "Non-coverage probability",
     main = "nominal level = 0.1")
abline(h=0.1, lty=2, col=1)
lines(x=rho.seq, y=ncp$`nominal=0.1`[,"OLS"], 
      type = "b", lwd=2, pch=5,
      col="magenta")
lines(x=rho.seq, y=ncp$`nominal=0.1`[,"EHW"], 
      type = "b", lwd=2, pch=5,
      col="brown")


plot(x=rho.seq, y=ncp$`nominal=0.05`[,"HAC"], 
     type = "b", lwd=2, pch=19,
     col=2,
     ylim = c(min(ncp$`nominal=0.05`), max(ncp$`nominal=0.05`)),
     xlab = expression(phi),
     ylab = "Non-coverage probability",
     main = "nominal level = 0.05")
abline(h=0.05, lty=2, col=1)
lines(x=rho.seq, y=ncp$`nominal=0.05`[,"OLS"], 
      type = "b", lwd=2, pch=5,
      col="magenta")
lines(x=rho.seq, y=ncp$`nominal=0.05`[,"EHW"], 
      type = "b", lwd=2, pch=5,
      col="brown")

plot(x=rho.seq, y=ncp$`nominal=0.01`[,"HAC"], 
     type = "b", lwd=2, pch=19,
     col=2,
     ylim = c(min(ncp$`nominal=0.01`), max(ncp$`nominal=0.01`)),
     xlab = expression(phi),
     ylab = "Non-coverage probability",
     main = "nominal level = 0.01")
abline(h=0.01, lty=2, col=1)
lines(x=rho.seq, y=ncp$`nominal=0.01`[,"OLS"], 
      type = "b", lwd=2, pch=5,
      col="magenta")
lines(x=rho.seq, y=ncp$`nominal=0.01`[,"EHW"], 
      type = "b", lwd=2, pch=5,
      col="brown")

# use ggplot with given nominal level #---------------------------
ncp.HAC = data.frame(ncp = ncp$`nominal=0.05`[,"HAC"],
                     rho = rho.seq,
                     nominal = "nominal level: 0.05")
library(ggplot2)
library("ggpubr")
library(ggh4x)

ggplot(ncp.HAC, aes(x = rho, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, color = 2, pch = 19) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  labs(x = expression(phi), y = " ") +
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = seq(0, 0.1, by = 0.025)) +
  theme_bw()

ncp.frame = data.frame(
  ncp = c( ncp$`nominal=0.05`[,"HAC"], ncp$`nominal=0.05`[,"EHW"]  ),
  rho = c(rho.seq, rho.seq),
  SE = rep(c("HAC SE", "robust SE"), each = length(rho.seq)),
  nominal = "nominal level: 0.05"
)

ggplot(ncp.frame, aes(x = rho, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(aes(color = SE, shape = SE, group = SE),
                 size = 2, linewidth = 0.8) + 
  scale_color_manual(values = c("red", "black")) +
  scale_shape_manual(values = c(17, 19)) + 
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  labs(x = expression(phi), y = " ")+ 
  scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  theme_bw() +
  theme(legend.position = "none") # hide the legend


ggplot(ncp.frame, aes(x = rho, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) + 
  geom_point(aes(color = SE)) +
  geom_line(aes(color = SE, linetype=SE), linewidth=0.8) + 
  geom_hline(yintercept=0.05,color = "black", linetype=2)+
  scale_y_continuous(breaks=c(0.05, 0.15, 0.3, 0.45, 0.6)) +
  labs(x=expression(phi), y=" ") +
  scale_color_manual(values = c("red", "black")) + 
  scale_linetype_manual(values = c(1,2)) +
  theme_bw() + 
  theme(legend.position = "none") # hide the legend




cp.frame = data.frame(
  cp = 1-c( ncp$`nominal=0.05`[,"HAC"], ncp$`nominal=0.05`[,"EHW"]  ),
  rho = c(rho.seq, rho.seq),
  SE = rep(c("HAC SE", "robust SE"), each = length(rho.seq)),
  nominal = "nominal level: 0.95"
)


ggplot(cp.frame, aes(x = rho, y = cp)) +
  facet_wrap(~nominal, ncol = 1) + 
  geom_point(aes(color = SE)) +
  geom_line(aes(color = SE, linetype=SE), linewidth=0.8) + 
  geom_hline(yintercept=0.95,color = "black", linetype=2)+
  scale_y_continuous(breaks=c(0.1, 0.3, 0.5, 0.75, 0.9, 0.95)) +
  labs(x=expression(phi), y=" ") +
  scale_color_manual(values = c("red", "black")) + 
  scale_linetype_manual(values = c(1,2)) +
  theme_bw() + 
  theme(legend.position = "none") # hide the legend
















































































































































###########################################
####  EXP2 in Motivation of the Intro #####
#####  For the validity model LZ in X #####
###########################################
rm(list = ls())


library("sandwich")
library("lmtest")
library(Matrix)

grp.size_eps.seq = c(5, seq(100,1000,by=100))

N=5000

d = 1

beta = 1

gamma = 1 # the first one is the target in the outcome model

run = 1000


t_ols = matrix(0, nrow = run, ncol = length(grp.size_eps.seq))

t_ehw = matrix(0, nrow = run, ncol = length(grp.size_eps.seq))

t_LZ.X = matrix(0, nrow = run, ncol = length(grp.size_eps.seq))
se2_LZ.X = matrix(0, nrow = run, ncol = length(grp.size_eps.seq))

t_LZ.eps = matrix(0, nrow = run, ncol = length(grp.size_eps.seq))
se2_LZ.eps = matrix(0, nrow = run, ncol = length(grp.size_eps.seq))

one.N = as.matrix(rep(1, times=N))



grp.size_X = 5
n.X = N/grp.size_X




blk_loading_X = replicate(n.X, matrix(1, nrow = grp.size_X, ncol=1), simplify = F)
loading_X = as.matrix(bdiag(blk_loading_X))






set.seed(715)
for (loop.epsgrp in 1:length(grp.size_eps.seq)) {
  
  grp.size_eps = grp.size_eps.seq[loop.epsgrp]
  n.eps = N/grp.size_eps
  
  if (floor(n.eps) == n.eps) {
    blk_loading_eps = replicate(n.eps, matrix(1, nrow = grp.size_eps, ncol=1), simplify = F)
    loading_eps = as.matrix(bdiag(blk_loading_eps))
  } else {
    blk_loading_eps = replicate(floor(n.eps), matrix(1, nrow = grp.size_eps, ncol=1), simplify = F)
    blk_loading_eps[[(floor(n.eps)+1)]] = matrix(1, nrow = (N-floor(n.eps)*grp.size_eps), ncol = 1 )
    loading_eps = as.matrix(bdiag(blk_loading_eps))
  }
  
  I.grp.size_X <- matrix(1, nrow = grp.size_X, ncol = grp.size_X)
  Wgt_LZ_X = kronecker(diag(n.X), I.grp.size_X) # the weight matrix of LZ SE if we use X clusters.
  
  I.grp.size_eps <- matrix(1, nrow = grp.size_eps, ncol = grp.size_eps)
  if (floor(n.eps) == n.eps) {
    Wgt_LZ_eps = kronecker(diag(n.eps), I.grp.size_eps) # the weight matrix of LZ SE if we use eps clusters
  } else {
    Wgt_LZ_eps = as.matrix( 
      bdiag(
        kronecker(diag(floor(n.eps)), I.grp.size_eps), 
        matrix(1, nrow = (N-floor(n.eps)*grp.size_eps), ncol = (N-floor(n.eps)*grp.size_eps) )
      )
    )
  }
  
  for (loop.run in 1:run) {
    
    
    
    # generate Z #---------------------------------------
    # this Z has intercept
    Z = one.N
    
    # generate X #---------------------------------------
    w_X = as.matrix(rnorm(N))
    f_X = as.matrix(rnorm(n.X))
    X = sqrt(0.9)*(loading_X%*%f_X) + sqrt(0.1) * w_X
    
    M = cbind(X, Z) 
    
    # generate eps #-------------------------------------
    f_eps = rnorm(ncol(loading_eps))
    # w_eps = runif(n, min = -sqrt(3), max = sqrt(3))
    # w_eps = 2*rbinom(n,1,0.5)-1
    # w_eps = rgamma(n, shape = 2, scale = (1/sqrt(2)) ) - sqrt(2)
    # w_eps = rgamma(n, shape = 1, scale = 1) - 1
    w_eps =  rnorm(N)
    eps = sqrt(0.9) * loading_eps%*%f_eps + sqrt(0.1) * w_eps
    
    
    # generate y #----------------------------------------
    y = beta * X +  gamma * Z + eps
    
    
    # ols #-----------------------------------------------
    fit.ols=lm(y~M-1)
    t_ols[loop.run,loop.epsgrp] = ( summary(fit.ols)$coefficients[1, "Estimate"] - beta ) / summary(fit.ols)$coefficients[1, "Std. Error"]
    #-----------------------------------------------------
    
    
    
    # EHW #-----------------------------------------------
    vcov.ehw=solve(t(M)%*%M) %*% (t(M) %*% diag(fit.ols$residuals)%*%diag(fit.ols$residuals) %*% M) %*% solve(t(M)%*%M)
    t_ehw[loop.run,loop.epsgrp] = ( summary(fit.ols)$coefficients[1, "Estimate"] - beta ) / sqrt(vcov.ehw[1,1])
    #-----------------------------------------------------
    
    
    
    # LZ using X groups #-------------------------------------------------------
    vcov.lz_X = solve(t(M)%*%M)%*%(t(M)%*%((fit.ols$residuals%*%t(fit.ols$residuals))*Wgt_LZ_X)%*%M)%*%solve(t(M)%*%M)
    t_LZ.X[loop.run,loop.epsgrp] = ( summary(fit.ols)$coefficients[1, "Estimate"] - beta ) / sqrt(vcov.lz_X[1,1])
    se2_LZ.X[loop.run, loop.epsgrp] = vcov.lz_X[1,1]
    #---------------------------------------------------------------------------
    
    
    
    
    
    
    
    # LZ using eps groups #-------------------------------------------------------
    vcov.lz_eps = solve(t(M)%*%M)%*%(t(M)%*%((fit.ols$residuals%*%t(fit.ols$residuals))*Wgt_LZ_eps)%*%M)%*%solve(t(M)%*%M)
    t_LZ.eps[loop.run,loop.epsgrp] = ( summary(fit.ols)$coefficients[1, "Estimate"] - beta ) / sqrt(vcov.lz_eps[1,1])
    se2_LZ.eps[loop.run, loop.epsgrp] = vcov.lz_eps[1,1]
    #---------------------------------------------------------------------------
    
    
    
    print(loop.run)
    
  }
  
  print(loop.epsgrp)
  
}


save(list = c("se2_LZ.eps", "se2_LZ.X", "t_ehw", "t_LZ.eps", "t_LZ.X", "t_ols",
              "beta", "d", "gamma", "grp.size_eps.seq", "grp.size_X", "N", "n.X", "run"),
     file = "LZ_X_vary_epsgrp.Rdata")


rm(loading_eps, blk_loading_eps, loading_X, blk_loading_X)


par(mfrow = c(1, 1))
plot(eps, type="b")
plot(X, type="b")

# draw qqplots with plot #------------------------------------------------------
prob = seq(0.01,0.99, by=0.01)

loop.epsgrp=6

qt.emp = quantile(t_ehw[,loop.epsgrp], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("EHW", " epsilon group size=", grp.size_eps.seq[loop.epsgrp]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_ols[,loop.epsgrp], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("OLS", " epsilon group size=", grp.size_eps.seq[loop.epsgrp]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.X[,loop.epsgrp], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("LZ_X", " epsilon group size=", grp.size_eps.seq[loop.epsgrp]))
abline(a=0, b=1, col=2, lwd=3)

qt.emp = quantile(t_LZ.eps[,loop.epsgrp], probs = prob)
plot(x = qt.emp, y = qnorm(prob), main = paste0("LZ_eps", " epsilon group size=", grp.size_eps.seq[loop.epsgrp]))
abline(a=0, b=1, col=2, lwd=3)


# make plot about non-coverage prob #------------------
qtl = qnorm(c(0.05, 0.025, 0.005), lower.tail = F) # non-coverage prob of 0.1, 0.05, 0.01

ncp = list() # This is a list of matrices. Each matrix corresponds to a qtl. Given qtl, each col represents each method. Each row is the empirical prob of |t_j| > qtl for some rho

for (loop.qtl in 1:length(qtl)){
  
  ncp.qtl = matrix(0, nrow = length(grp.size_eps.seq), ncol = 4)
  
  colnames(ncp.qtl) = c("OLS", "EHW", "LZ_X", "LZ_eps")
  rownames(ncp.qtl) = grp.size_eps.seq
  
  for (loop.epsgrp in 1:length(grp.size_eps.seq)) {
    
    ncp.qtl[loop.epsgrp, "OLS"] = sum(abs(t_ols[,loop.epsgrp]) > qtl[loop.qtl]) / run
    ncp.qtl[loop.epsgrp, "EHW"] = sum(abs(t_ehw[,loop.epsgrp]) > qtl[loop.qtl]) / run
    ncp.qtl[loop.epsgrp, "LZ_X"] = sum(abs(t_LZ.X[,loop.epsgrp]) > qtl[loop.qtl]) / run
    ncp.qtl[loop.epsgrp, "LZ_eps"] = sum(abs(t_LZ.eps[,loop.epsgrp]) > qtl[loop.qtl]) / run
    
  }
  
  ncp[[loop.qtl]] = ncp.qtl
  
  rm(ncp.qtl)
  
}
names(ncp) = c("nominal=0.1", "nominal=0.05", "nominal=0.01")

plot(x=grp.size_eps.seq, y=ncp$`nominal=0.1`[,"LZ_X"], 
     type = "b", lwd=2, pch=19,
     col=2,
     ylim = c(min(ncp$`nominal=0.1`), max(ncp$`nominal=0.1`)),
     xlab = expression(group~size~of~epsilon),
     ylab = "Non-coverage probability",
     main = "nominal level = 0.1")
abline(h=0.1, lty=2, col=1)
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.1`[,"LZ_eps"], 
      type = "b", lwd=2, pch=5,
      col=4)
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.1`[,"OLS"], 
      type = "b", lwd=2, pch=5,
      col="magenta")
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.1`[,"EHW"], 
      type = "b", lwd=2, pch=5,
      col="brown")


plot(x=grp.size_eps.seq, y=ncp$`nominal=0.05`[,"LZ_X"], 
     type = "b", lwd=2, pch=19,
     col=2,
     ylim = c(min(ncp$`nominal=0.05`), max(ncp$`nominal=0.05`)),
     xlab = expression(group~size~of~epsilon),
     ylab = "Non-coverage probability",
     main = "nominal level = 0.05")
abline(h=0.05, lty=2, col=1)
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.05`[,"LZ_eps"], 
      type = "b", lwd=2, pch=5,
      col=4)
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.05`[,"OLS"], 
      type = "b", lwd=2, pch=5,
      col="magenta")
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.05`[,"EHW"], 
      type = "b", lwd=2, pch=5,
      col="brown")

plot(x=grp.size_eps.seq, y=ncp$`nominal=0.01`[,"LZ_X"], 
     type = "b", lwd=2, pch=19,
     col=2,
     ylim = c(min(ncp$`nominal=0.01`), max(ncp$`nominal=0.01`)),
     xlab = expression(group~size~of~epsilon),
     ylab = "Non-coverage probability",
     main = "nominal level = 0.01")
abline(h=0.01, lty=2, col=1)
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.01`[,"LZ_eps"], 
      type = "b", lwd=2, pch=5,
      col=4)
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.01`[,"OLS"], 
      type = "b", lwd=2, pch=5,
      col="magenta")
lines(x=grp.size_eps.seq, y=ncp$`nominal=0.01`[,"EHW"], 
      type = "b", lwd=2, pch=5,
      col="brown")

# use ggplot with given nominal level #---------------------------
ncp.LZ_X = data.frame(ncp = ncp$`nominal=0.05`[,"LZ_X"],
                      grp.size = grp.size_eps.seq,
                      number_of_grp = ceiling(N/grp.size_eps.seq),
                      nominal = "nominal level: 0.05")
library(ggplot2)
library("ggpubr")
library(ggh4x)

ncp.LZ_X$number_of_grp = factor(ncp.LZ_X$number_of_grp)

ggplot(ncp.LZ_X, aes(x = grp.size, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, color = 2, pch = 19) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  labs(x=expression(group~size~of~epsilon), y = " ") +
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = seq(0, 0.1, by = 0.025)) +
  theme_bw()


ggplot(ncp.LZ_X, aes(x = number_of_grp, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, color = 2, pch = 19) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  labs(x=expression(number~of~groups~of~epsilon), y = " ") +
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = seq(0, 0.1, by = 0.025)) +
  theme_bw()

# check which SE is smaller #--------------------------------------
plot(x=grp.size_eps.seq, y = apply(se2_LZ.eps, 2, mean),
     type = "b", col=4,
     ylab = "average of SE2",
     xlab = "group size in eps")
lines(x=grp.size_eps.seq, y = apply(se2_LZ.X,2,mean),
      type = "b", col=2)
boxplot(cbind(se2_LZ.eps, se2_LZ.X))


ncp.frame = data.frame(
  ncp = c( ncp$`nominal=0.05`[,"LZ_X"], ncp$`nominal=0.05`[,"OLS"]  ),
  clstr.size = c(grp.size_eps.seq, grp.size_eps.seq),
  number_of_clstr = c(ceiling(N/grp.size_eps.seq), ceiling(N/grp.size_eps.seq)  ),
  SE = rep(c("cluster-robust", "non-robust"), each = length(grp.size_eps.seq)),
  nominal = "nominal level: 0.05"
)

ncp.frame$number_of_clstr = factor(ncp.frame$number_of_clstr
)

ggplot(ncp.frame, aes(x = clstr.size, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(aes(color = SE, shape = SE, group = SE),
                 size = 2, linewidth = 0.8) + 
  scale_color_manual(values = c("blue", "black")) +
  scale_shape_manual(values = c(15, 19)) + 
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.2, 0.3)) + 
  labs(x = expression(cluster~size~of~epsilon), y = " ") +
  theme_bw() +
  theme(legend.position = "none") # hide the legend

ggplot(ncp.frame, aes(x = number_of_clstr, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(aes(color = SE, shape = SE, group = SE),
                 linewidth = 0.8) + 
  scale_color_manual(values = c("blue", "black")) +
  scale_shape_manual(values = c(15, 19)) + 
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.2, 0.3)) + 
  labs(x = expression(number~of~clusters~of~epsilon), y = " ") +
  theme_bw() +
  theme(legend.position = "none") # hide the legend


ggplot(ncp.frame, aes(x = clstr.size, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) + 
  geom_point(aes(color = SE)) +
  geom_line(aes(color = SE, linetype=SE), linewidth=0.8) + 
  geom_hline(yintercept=0.05,color = "black", linetype=2)+
  scale_y_continuous(breaks=c(0.05, 0.15, 0.3)) +
  labs(x=expression(cluster~size~of~epsilon), y=" ") +
  scale_color_manual(values = c("blue", "black")) + 
  scale_linetype_manual(values = c(1,2)) +
  theme_bw() + 
  theme(legend.position = "none") # hide the legend


ggplot(ncp.frame, aes(x = number_of_clstr, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_point(aes(color = SE)) +
  geom_line(aes(color = SE, linetype=SE, group = SE), linewidth=0.8) + 
  scale_color_manual(values = c("blue", "black")) +
  scale_shape_manual(values = c(15, 19)) + 
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.2, 0.3)) + 
  labs(x = expression(number~of~clusters~of~epsilon), y = " ") +
  theme_bw() +
  theme(legend.position = "none") # hide the legend




cp.frame = data.frame(
  cp = 1-c( ncp$`nominal=0.05`[,"LZ_X"], ncp$`nominal=0.05`[,"EHW"]  ),
  clstr.size = c(grp.size_eps.seq, grp.size_eps.seq),
  number_of_clstr = c(ceiling(N/grp.size_eps.seq), ceiling(N/grp.size_eps.seq)  ),
  SE = rep(c("cluster-robust", "robust"), each = length(grp.size_eps.seq)),
  nominal = "nominal level: 0.95"
)

cp.frame$number_of_clstr = factor(ncp.frame$number_of_clstr
)

ggplot(cp.frame, aes(x = number_of_clstr, y = cp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_point(aes(color = SE)) +
  geom_line(aes(color = SE, linetype=SE, group = SE), linewidth=0.8) + 
  scale_color_manual(values = c("blue", "black")) +
  scale_shape_manual(values = c(15, 19)) + 
  geom_hline(yintercept = 0.95, color = "black", linetype = 2) +
  scale_y_continuous(breaks=c(0.6, 0.7, 0.8, 0.9, 0.95)) + 
  labs(x = expression(number~"of"~clusters~"in"~epsilon), y = " ") +
  theme_bw() +
  theme(legend.position = "none") # hide the legend









