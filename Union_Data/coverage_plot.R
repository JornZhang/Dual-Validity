#########################################################################
########### AR1 structured Xperp and eps #################### 
########### Add aggregate shock by loading1 and loading2 to eps #########
#########################################################################
rm(list = ls())

shock.seq = c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9")
task.seq = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10")

## allocate final matrices
t_HAC_eps_shock  <- matrix(NA_real_, nrow = 10000, ncol = length(shock.seq))
t_HAC_HC1_eps_shock <- matrix(NA_real_, nrow = 10000, ncol = length(shock.seq))
df_HAC_eps_shock <- matrix(NA_real_, nrow = 10000, ncol = length(shock.seq))

t_HAC_Xperp_shock  <- matrix(NA_real_, nrow = 10000, ncol = length(shock.seq))
t_HAC_HC1_Xperp_shock <- matrix(NA_real_, nrow = 10000, ncol = length(shock.seq))
df_HAC_Xperp_shock <- matrix(NA_real_, nrow = 10000, ncol = length(shock.seq))

#######################################################################
Lag_eps = 20
Lag_X = 20
#######################################################################

for (s in seq_along(shock.seq)) {
   
  res_eps_shock = readRDS(paste0("results_eps_shock", 
                            "/Shock_on_eps_L", Lag_eps, "_task", task.seq[s], "_kappa", shock.seq[s],".rds"))
  
  t_HAC_eps_shock[,  s] <- res_eps_shock$t
  t_HAC_HC1_eps_shock[,s]  = res_eps_shock$t_HC1
  df_HAC_eps_shock[, s] <- res_eps_shock$df
  
  
  res_Xperp_shock = readRDS(paste0("results_Xperp_shock", 
                            "/Shock_on_Xperp_L", Lag_X, "_task", task.seq[s], "_kappa", shock.seq[s],".rds"))
  
  t_HAC_Xperp_shock[,  s] <- res_Xperp_shock$t
  t_HAC_HC1_Xperp_shock[,s]  = res_Xperp_shock$t_HC1
  df_HAC_Xperp_shock[, s] <- res_Xperp_shock$df
  
}


nominal = 0.05

cp_HAC_eps_shock = colSums(abs(t_HAC_eps_shock) < qnorm(1 - nominal/2)) / nrow(t_HAC_eps_shock)

cp_HAC_HC1_eps_shock = colSums(abs(t_HAC_HC1_eps_shock) < qnorm(1 - nominal/2)) / nrow(t_HAC_eps_shock)

cp_HAC_HC1_Satterthwaite_eps_shock = colSums(abs(t_HAC_HC1_eps_shock) < qt( (1 - nominal/2), df = df_HAC_eps_shock)) / nrow(t_HAC_eps_shock)




cp_HAC_Xperp_shock = colSums(abs(t_HAC_Xperp_shock) < qnorm(1 - nominal/2)) / nrow(t_HAC_Xperp_shock)

cp_HAC_HC1_Xperp_shock = colSums(abs(t_HAC_HC1_Xperp_shock) < qnorm(1 - nominal/2)) / nrow(t_HAC_Xperp_shock)

cp_HAC_HC1_Satterthwaite_Xperp_shock = colSums(abs(t_HAC_HC1_Xperp_shock) < qt( (1 - nominal/2), df = df_HAC_Xperp_shock)) / nrow(t_HAC_Xperp_shock)


apply(df_HAC_eps_shock, 2, mean)
apply(df_HAC_Xperp_shock, 2, mean)


shock.seq = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

plot(x = shock.seq,
     y = cp_HAC_eps_shock,
     type = "b",
     col=2,
     pch=19,
     xlab= "kappa",
     ylab= "coverage probability",
     ylim=c(0.9, 0.95))
lines(x = shock.seq,
       y = cp_HAC_Xperp_shock,
       type = "b",
       col=4,
       pch=19)
abline(h=0.95, lty=2)






plot(x = shock.seq,
     y = cp_HAC_HC1_eps_shock,
     type = "b",
     col=2,
     pch=19,
     xlab= "kappa",
     ylab= "coverage probability",
     ylim=c(0.9, 0.95))
lines(x = shock.seq,
      y = cp_HAC_HC1_Xperp_shock,
      type = "b",
      col=4,
      pch=19)
abline(h=0.95, lty=2)





plot(x = shock.seq,
     y = cp_HAC_HC1_Satterthwaite_eps_shock,
     type = "b",
     col=2,
     pch=19,
     xlab= "kappa",
     ylab= "coverage probability",
     ylim=c(0.8, 0.95))
lines(x = shock.seq,
      y = cp_HAC_HC1_Satterthwaite_Xperp_shock,
      type = "b",
      col=4,
      pch=19)
abline(h=0.95, lty=2)








df_cp = data.frame(
  shock.seq = rep(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                  times = 2),
  cp = c(cp_HAC_HC1_Satterthwaite_eps_shock,
              cp_HAC_HC1_Satterthwaite_Xperp_shock),
  nominal = rep("nominal level: 0.95", times = (2*length(shock.seq)) ),
  type = rep(c("shock on eps",
               "shock on X"), each = length(cp_HAC_HC1_Satterthwaite_eps_shock))
)

library(ggplot2)
library("ggpubr")
library(ggh4x)
ggplot(df_cp, aes(x = shock.seq, y = cp, color = type)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, 
                 linewidth = 0.8, 
                 shape = 19) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.85, 1))




ggplot(df_cp, aes(x = shock.seq, y = cp, 
                  color = type, 
                  linetype = type)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, shape = 19) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.85, 1))


ggplot(df_cp, aes(x = shock.seq, y = cp, 
                  color = type, 
                  linetype = type)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, shape = 19) +
  scale_linetype_manual(values = c("solid", 
                                   "dashed")) + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.85, 1))


ggplot(df_cp, aes(x = shock.seq, y = cp, 
                  color = type, 
                  linetype = type)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_line(linewidth = 0.8) +
  scale_linetype_manual(values = c("solid", 
                                   "dashed")) + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.85, 1))



