#############################################################
#### This file uses the Bell-McCaffery correction ###########
#############################################################

rm(list=ls())
library(dplyr)
library(clubSandwich)
library(fwildclusterboot)


# folder = "D:\\Research\\20251214_PUMS"
folder = Sys.getenv("PUMS_DATA_DIR", unset = "data")
load(file.path(folder, "1222_PUMA_setting.Rdata"))


rm(list = setdiff(ls(), c("dat_tgt", "folder")))
# save(list = c("dat_tgt"), file = "dat_tgt.Rdata")
# load("dat_tgt.Rdata")

##------------------------------------------------#
# Function: draw one value per group, then broadcast it to all rows in that group
# Input: grp must be a factor indicating groups.
# mean and sd are the rnorm arguments
# Output: A vector with the same length as grp. The indices in the same group get the same value from N(mean, sd^2).
rnorm_by_group_factor <- function(grp, mean = 0, sd = 1) {
  u <- rnorm(nlevels(grp), mean, sd)
  u[as.integer(grp)]
}
# --------------------------------------------# 


run = 100

B.boot = 100 # This is the number of bootstrap samples



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
# sigm_XD.seq = c(0, 0.05, 0.1, 0.15, 0.2, 0.4) 
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

colnames(Z)[1] = "one"

colnames(Z)

rm(dat_tgt)

Z.matrix = as.matrix(Z)



#       
#           (Intercept)      HHL_Spn    HHL_Other         RMSP      YRBLT1940    YRBLT1950  
eta =   c(  10.202178,     0.009120,     0.058568,     0.058592,    -0.003536,     0.041945,  
            #YRBLT1960      YRBLT1970    YRBLT1980     YRBLT1990      YRBLT2000    YRBLT2010  
            0.042215,     0.032049,     0.084491,     0.121783,     0.160453,     0.203288,  
            #YRBLT2020      YRBLT2021     YRBLT2022     YRBLT2023    YRBLT2024        BLD1U  
            0.189724,     0.213712,     0.210512,     0.226751,     0.195161,     0.329970,  
            #BLDMU         AGEP            SEX         MSP_sgl       RAC1P2        RAC1P3  
            0.459882,    -0.008931,     0.292543,    -0.283887,    -0.153351,    -0.160468,  
            #RAC1P4       RAC1P5          RAC1P6       RAC1P7        RAC1P8        RAC1P9  
            -0.286003,    -0.149710,     0.078426,     0.007157,    -0.059985,    -0.022738,  
            #SCHL_HSCOMP  SCHL_COCOMP          DIS  
            0.301461,     0.839699,    -0.257458)

bet = 0.27604

#      (Intercept)     HHL_Spn     HHL_Other         RMSP       YRBLT1940    YRBLT1950  
gamm = c(7.237291,    0.124433,     0.182241,     0.057493,     0.069301,     0.143184,  
         #YRBLT1960     YRBLT1970     YRBLT1980     YRBLT1990     YRBLT2000     YRBLT2010     YRBLT2020  
         0.161374,     0.171950,     0.254317,     0.305908,     0.361102,     0.422672,     0.482539,  
         #YRBLT2021     YRBLT2022     YRBLT2023    YRBLT2024        BLD1U        BLDMU         AGEP  
         0.477168,     0.519450,     0.53289,     0.526381,     1.331883,     1.465468,     0.003567,  
         #SEX           MSP_sgl       RAC1P2       RAC1P3        RAC1P4        RAC1P5         RAC1P6  
         -0.003794,    -0.042660,    -0.190255,    -0.252357,    -0.489927,    -0.084555,     0.316191,  
         #RAC1P7        RAC1P8        RAC1P9      SCHL_HSCOMP  SCHL_COCOMP        DIS  
         0.351471,     0.110249,     0.056949,     0.088430,    0.337412,    -0.086665 )

Z.eta = Z.matrix%*%eta
Z.gamm = Z.matrix%*%gamm

rm(Z.matrix)


# --- PUMA level results --- #
t_LZ.P = matrix(0, nrow = run, length(sigm_XD.seq)) # original data t-statistic using adjusted LZ standard error, where the adjustment is the LZ SE times sqrt of (G / (G-1)) * ((N-1)/(N-p)), G is the number of cluster (PUMA), and N is the number of observations.
CI_L_LZ.P = matrix(0, nrow = run, length(sigm_XD.seq)) # 95% bootstrap confidence set, left end point
CI_R_LZ.P = matrix(0, nrow = run, length(sigm_XD.seq)) # 95% bootstrap confidence set, right end point
p_LZ.P = matrix(0, nrow = run, length(sigm_XD.seq)) # p-value from the bootstrap sample

# --- STATE level results --- #
t_LZ.S = matrix(0, nrow = run, length(sigm_XD.seq))
CI_L_LZ.S = matrix(0, nrow = run, length(sigm_XD.seq))
CI_R_LZ.S = matrix(0, nrow = run, length(sigm_XD.seq))
p_LZ.S = matrix(0, nrow = run, length(sigm_XD.seq))

# --- DIVISION level results --- #
t_LZ.D = matrix(0, nrow = run, length(sigm_XD.seq))
CI_L_LZ.D = matrix(0, nrow = run, length(sigm_XD.seq))
CI_R_LZ.D = matrix(0, nrow = run, length(sigm_XD.seq))
p_LZ.D = matrix(0, nrow = run, length(sigm_XD.seq))






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
    
    # fit #
    dat_sim = cbind(y, X, Z, DVSN_index, ST_index, PUMA_index)
    fit.Y = lm(y ~  0 + X + one + HHL_Spn + HHL_Other + RMSP + YRBLT1940 + YRBLT1950 + YRBLT1960 +
               YRBLT1970 + YRBLT1980 + YRBLT1990 + YRBLT2000 + YRBLT2010 + YRBLT2020 + YRBLT2021 +
               YRBLT2022 + YRBLT2023 + YRBLT2024 + BLD1U + BLDMU + AGEP + SEX +        
               MSP_sgl + RAC1P2 + RAC1P3 + RAC1P4 + RAC1P5 + RAC1P6 + RAC1P7 +   
               RAC1P8 + RAC1P9 + SCHL_HSCOMP + SCHL_COCOMP + DIS,
               data = dat_sim)
    
    # point est of beta
    bethat = summary(fit.Y)$coefficients["X", "Estimate"]
    
    
    # LZ using PUMA groups #-------------------------------------------------------
    boot_P = boottest(
      fit.Y, 
      clustid = "PUMA_index",
      R = c(1), 
      r=bet, # the null value
      param = "X",    # H0: t(R)%*%param = r, where r must be a scalar.
      B = B.boot # number of bootstrap samples
    )
    
    # summary(boot_P)
    t_LZ.P[loop.run,loop.s] = boot_P$t_stat
    CI_L_LZ.P[loop.run, loop.s] = boot_P$conf_int[1]
    CI_R_LZ.P[loop.run, loop.s] = boot_P$conf_int[2]
    p_LZ.P[loop.run, loop.s] = boot_P$p_val
    #---------------------------------------------------------------------------
    
    
    
    
    # LZ using STATE groups #-------------------------------------------------------
    boot_S = boottest(
      fit.Y,
      clustid = "ST_index",
      R = c(1),
      r=bet,
      param = "X",
      B= B.boot
    )
    
    # summary(boot_S)
    t_LZ.S[loop.run, loop.s] = boot_S$t_stat
    CI_L_LZ.S[loop.run, loop.s] = boot_S$conf_int[1]
    CI_R_LZ.S[loop.run, loop.s] = boot_S$conf_int[2]
    p_LZ.S[loop.run, loop.s] = boot_S$p_val
    #---------------------------------------------------------------------------
    
    
    # LZ using DIVISION groups #-------------------------------------------------------
    boot_D = boottest(
      fit.Y,
      clustid = "DVSN_index",
      R = c(1),
      r = bet,
      param = "X",
      B = B.boot # , type    = "webb"
    )
    
    # summary(boot_D)
    t_LZ.D[loop.run, loop.s] = boot_D$t_stat
    CI_L_LZ.D[loop.run, loop.s] = boot_D$conf_int[1]
    CI_R_LZ.D[loop.run, loop.s] = boot_D$conf_int[2]
    p_LZ.D[loop.run, loop.s] = boot_D$p_val
    #---------------------------------------------------------------------------
    
    
    print(loop.run)
    
  }
  
  print(loop.s)
  
}


save(list = c("CI_L_LZ.D", "CI_L_LZ.P", "CI_L_LZ.S", "CI_R_LZ.D", "CI_R_LZ.P", "CI_R_LZ.S",
              "p_LZ.D", "p_LZ.P", "p_LZ.S", "t_LZ.D", "t_LZ.P", "t_LZ.S",
              "Z", "B.boot", 
              "bet", "DVSN_index", "eta", "gamm", "PUMA_index", "run", 
              "sigm_epsD", "sigm_epse", "sigm_epsP", "sigm_epsS",
              "sigm_XD.seq", "sigm_Xe", "sigm_XP", "sigm_XS",
              "ST_index"),
     file = "0107PUMS_change_sigm_XD_boot_auto.Rdata")





# --------- draw the type I error plot using cluster wild bootstrap p values -------------- #
rm(list = ls())
load("0107PUMS_change_sigm_XD_boot_auto.Rdata")
nominal_level = 0.05

TypeI_boot.D = 1 * (p_LZ.D < nominal_level) # division level SE type I error indicator
TypeI_boot.S = 1 * (p_LZ.S < nominal_level) # state level SE type I error indicator
TypeI_boot.P = 1 * (p_LZ.P < nominal_level) # puma level SE type I error indicator

TypeI_rate = data.frame(
  rate = c(
    apply(TypeI_boot.D, 2, mean),
    apply(TypeI_boot.S, 2, mean),
    apply(TypeI_boot.P, 2, mean)
  ),
  Shock = c(sigm_XD.seq, sigm_XD.seq, sigm_XD.seq),
  Method = c(rep("DIVISION", times = length(sigm_XD.seq) ),
             rep("STATE", times = length(sigm_XD.seq) ),
             rep("PUMA", times = length(sigm_XD.seq) )
  ),
  nominal = paste0("nominal level = ", nominal_level)
)

TypeI_rate

TypeI_rate$Method = factor(TypeI_rate$Method, levels = c("DIVISION", "STATE", "PUMA"))

library(ggplot2)
library("ggpubr")
library(ggh4x)


ggplot(TypeI_rate, aes(x = Shock, y = rate)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(size = 2, linewidth = 0.8, pch = 19, 
                 aes(color=Method, group=Method, linetype = Method)) +
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_color_manual(values = c("orange",  "magenta", "brown"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) + 
  labs(x=expression(sigma[XD]), y = " ")  +
  scale_x_continuous(limits = range(sigm_XD.seq),
                     breaks = sigm_XD.seq) + 
  theme_bw() 


ggplot(TypeI_rate, aes(x = Shock, y = rate)) +
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
    data = data.frame(sigm_XD_real = 0.079, ncp = -0.02, nominal = unique(TypeI_rate$nominal)),
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


# ------ coverage prob ----------------- #
cp.frame = TypeI_rate
cp.frame$cp = 1 - TypeI_rate$rate
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


