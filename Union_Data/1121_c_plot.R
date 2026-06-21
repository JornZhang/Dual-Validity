##################################################################
################ Draw some plot ##################################
##################################################################









# ------------------------- Draw the U-shape union density vs social inequality ---------------------#
rm(list = ls())

# --- 0.1) Libraries --- #
library(haven)
library(dplyr)
library(stringr)
library(tidyr)
library(fixest)




# --- 1.1) Gallup + trim --- #
gallup <- read_dta("composite_gallup_boot.dta") %>%
  select(statefip, year, union, wgt) %>%
  mutate(across(c(statefip, year), as.integer),
         across(c(union, wgt), as.numeric)) %>%
  filter(!is.na(statefip), !is.na(year), !is.na(union), !is.na(wgt), wgt > 0)

# the original "composite_gallup_boot.dta" has 8 columns, with one column named "dataset" with lables: 
# 0   ANES
# 1    CPS
# 2 GALLUP
# 3 BIKINI
# 4   1936
# 5   PSID
# And the "dataset" variable in "composite_gallup_boot.dta" only contains value "2".

# --- 1.2) CPS + trim --- #
cps <- read_dta("cps_for_union_boot.dta") %>%
  select(statefip, year, union_hh, wgt) %>%
  transmute(
    statefip = as.integer(statefip),
    year     = as.integer(year),
    union    = as.numeric(union_hh),   
    wgt      = as.numeric(wgt)
  ) %>%
  filter(!is.na(statefip), !is.na(year), !is.na(union), !is.na(wgt), wgt > 0)

# --- 1.3) State-year outcome/control panel --- #
panel <- read_dta("state_year_union_ineq.dta") %>%
  mutate(across(c(statefip, year), as.integer))

# --- 2) Construct state-year union density, without sample splitting --- #
# --- 2.1) Union density for gallup data --- #
g_g = gallup %>%
  group_by(statefip, year) %>%
  summarise(
    n = n(),
    union_density = stats::weighted.mean(union, wgt, na.rm = TRUE)
  )

# --- 2.2) Union density for cps data --- #
g_c = cps %>%
  group_by(statefip, year) %>%
  summarise(
    n = n(),
    union_density = stats::weighted.mean(union, wgt, na.rm = TRUE)
  )


# --- 2.3) Merge gallup and cps state-year union density data --- #
gallup_cps_union <-  full_join(g_g, g_c,
                               by = c("statefip","year"),
                               suffix = c(".g",".c")) %>%   # # Combine gallup_split and cps_split on the keys statefip and year. The function full_join keeps all (state, year) pairs that appear in either dataset (the union of keys). The argument suffix = c(".g",".c") disambiguates same-named columns: columns coming from Gallup get .g, from CPS get .c.
  transmute( # transmute() builds a new data frame with just the variables we specify.
    statefip, year, # We keep statefip, year, and create “spliced” columns.
    # Means: prefer CPS when available, else Gallup
    union_st_density = if_else(!is.na(union_density.c), union_density.c, union_density.g),
    # Cell sizes aligned with the chosen source for each half
    obs_union_density       = case_when(
      !is.na(union_density.c) ~ as.double(n.c),  # If the CPS split-mean exists (not NA), prefer CPS and take the CPS cell size n.c.
      !is.na(union_density.g) ~ as.double(n.g),  # Otherwise, if the Gallup split-mean exists, use the Gallup cell size n.g.
      TRUE ~ NA_real_    # If neither CPS nor Gallup has a valid split-mean (both missing), set the cell size to a typed numeric missing (NA_real_).
    )
  )


# --- 3) Define the southern states --- #
#                 DE,DC,FL,GA,MD,NC,VA,WV, AL,KY,MS,TN,AR,LA,OK,TX
south_states <- c(10,11,12,13,24,37,51,54, 1, 21,28,47,5, 22,40,48)


# --- 4) Merge the panel data containing state-year dependent variables with galluup_sps_union, which contains the state-year union density.
dat0 <- panel %>%
  inner_join(gallup_cps_union, by = c("statefip","year")) %>% # Keep only rows whose (statefip, year) exist in both panel and iv_df. inner_join keeps only overlaps.
  mutate(
    south    = as.integer(statefip %in% south_states), # Create a 0/1 indicator for “South” based on whether statefip is in the vector south_states.
    statefip = as.integer(statefip),
    year     = as.integer(year)
  )



# --- 4.1) Draw the plot of union density v.s. top10 --- #
# The plot shows the annual union density and top10 share, averaged across states
# So we first construct the data frame averaging over states
union_top10_ave <- dat0 %>% group_by(year) %>% summarise(
  ave_union = mean(union_st_density),
  ave_top10 = mean(top10)
) 
print(union_top10_ave, n=100)


# draw the plot of union density v.s. top10 #
plot(x=union_top10_ave$year, 
     y=union_top10_ave$ave_top10, type = "b", col=2,
     pch = 19,
     ylim = range( 
       c(
         range(union_top10_ave$ave_union*100), 
         range(union_top10_ave$ave_top10[!is.na(union_top10_ave$ave_top10)])
       ) 
     ),
     xlab = " ",
     ylab = "Percent",
     main = "Union Density v.s. Top-10 Share"
)
lines(x=union_top10_ave$year, 
      y=union_top10_ave$ave_union*100, type = "b", col=4,
      pch = 15)
legend("right",
       legend = c("Top-10 Share", "Union Density"),
       col    = c(2, 4),
       lwd    = c(1,1), 
       pch = c(19, 15), 
       cex = 0.9)


union_top10_ave_gg = data.frame(
  value = c(
    (100*union_top10_ave$ave_union[1:79]),
    union_top10_ave$ave_top10[1:79]
  ),
  year = c(1937:2015, 1937:2015),
  label = c(
    rep("Union Density", times = 79),
    rep("Top-10% Share", times = 79)
  )
)

library(ggplot2)
library(dplyr)

ggplot(union_top10_ave_gg, aes(x = year, y = value, color = label, group = label)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Union Density" = "red",
                                "Top-10% Share" = "blue")) +
  labs(x = "Year", y = NULL, color = NULL) +
  theme_bw() +
  theme(
    legend.position = c(0.98, 0.65),            # inside: right/top-ish
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                     colour = "grey80"),
    legend.key = element_blank(),
    legend.text = element_text(size = 10),  # make text larger
    legend.key.size = unit(1, "lines"),                  # enlarge color lines
    legend.spacing.y = unit(0.1, "cm")                     # vertical spacing
  )

rm(list = ls())























library(viridis)

# ------------ draw the ACF plot --------------------- #
rm(list = ls())
load("1031_PermSimu_Setting_c.Rdata")

library(latex2exp)

lag_max = 40
off = 0.12

# abs(ACF) AVERAGED ACROSS STATES
plot( x = ( (1:lag_max) - off ), 
      y =  apply( ( abs(res_c_acf_mat$ACF_mat[-1,]) ), 1,mean ) ,
      type = "h", lwd = 3, lend = "butt",
      xlab = "Lag", ylab = " ",
      ylim = c(-0.1,0.75), 
      xlim = c((0-3*off), lag_max+2*off),
      xaxs = "i", 
      yaxs = "i",
      col  = "#1f77b4",
      main = " "
)
lines(x = ( (1:lag_max) + off ),
      y =  apply( ( abs(Xperp_c_acf_mat$ACF_mat[-1,]) ), 1, mean) , 
      type = "h", lwd = 3, lend = "butt", col = "#d62728")
points(x=( (1:lag_max) - off ), 
       apply( ( abs(res_c_acf_mat$ACF_mat[-1,]) ), 1,mean ), 
       pch = 16, col = "#1f77b4")
points(x = ( (1:lag_max) + off ), 
       apply( ( abs(Xperp_c_acf_mat$ACF_mat[-1,]) ), 1, mean) ,
       pch = 16, col = "#d62728")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
legend("topright",
       legend = TeX(c("$\\widehat{X}_{\\perp}$", "$\\widehat{\\epsilon}$")),
       col    = c("#d62728", "#1f77b4"),
       lwd    = 3, 
       pch = 16, 
       bty = "n", 
       cex    = 1.2,      # ↑ text size (1.0 is default)
       pt.cex = 0.8       # symbol size in legend (independent of text)
)



library(ggplot2)

acf_df = data.frame(
  acf_value = c(
    apply( ( abs(res_c_acf_mat$ACF_mat[-1,]) ), 1,mean ),
    apply( ( abs(Xperp_c_acf_mat$ACF_mat[-1,]) ), 1, mean)
  ),
  label = factor(c(rep("hat(epsilon)", lag_max), rep("Xperp", lag_max)),
                 levels = c("Xperp", "hat(epsilon)")),   # red first, blue second
  lag = rep(1:lag_max, 2)
)
# Dotted-line threshold
thr <- qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used))

pos <- position_dodge(width = 0.35)

ggplot(acf_df, aes(x = lag, y = acf_value, color = label)) +
  geom_segment(aes(xend = lag, y = 0, yend = acf_value),
               position = pos, linewidth = 0.7) +
  geom_point(position = pos, size = 2) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "grey40") +
  scale_color_manual(
    name = NULL,
    values = c("Xperp" = "red", "hat(epsilon)" = "blue"),
    labels = c(expression(X["\u22A5"]), expression(hat(epsilon)))
  ) +
  labs(x = "Lag", y = NULL) +
  theme_bw() +
  theme(
    legend.position = c(0.99, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                     colour = "grey80"),
    legend.key = element_blank(),
    legend.text = element_text(size = 12),  # larger legend text
    legend.key.size = unit(1.5, "lines"),                  # larger color line
    legend.spacing.y = unit(0.4, "cm")                     # more vertical spacing
  )













# ----------- Q statistic ------------------------- #
h.seq_Q = c(10:20)


# The raw statistic
colnames(Q_res_c)= h.seq_Q
colnames(Q_Xperp_c) = h.seq_Q

loop.h = 1
plot(x=Q_Xperp_c[,loop.h], y = Q_res_c[,loop.h], 
     xlab = expression(X_(perp)), ylab = expression(hat(epsilon)),
     main = paste0("Ljung-Box Q statsitic with h=",h.seq_Q[loop.h])
)
abline(a=0,b=1,col=2)


library(ggplot2)

df_tmp <- data.frame(x = Q_Xperp_c[,loop.h], y = Q_res_c[,loop.h])

ggplot(df_tmp, aes(x, y)) +
  geom_point(alpha = 0.7, size = 2, col="blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1, linetype = "solid")  + 
  labs(
    x = expression(Q[X,i]),          # ⊥ as Unicode subscript
    y = expression(Q[epsilon,i]),          # \widehat{ε}
    title = paste0("Ljung-Box Q statsitic with h=",h.seq_Q[loop.h])
  )  +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),  # axis label size
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text     = element_text(size = 10)                 # tick labels (optional)
  )

ggplot(df_tmp, aes(x, y)) +
  geom_point(alpha = 0.7, size = 2, col="blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1, linetype = "solid")  + 
  labs(
    x = expression(Q[X,i]),          # ⊥ as Unicode subscript
    y = expression(Q[epsilon,i]),          # \widehat{ε}
  )  +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),  # axis label size
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text     = element_text(size = 10)                 # tick labels (optional)
  )


# the p values
colnames(Q_p_res_c)= h.seq_Q
colnames(Q_p_Xperp_c) = h.seq_Q

plot(x=Q_p_Xperp_c[,loop.h], y = Q_p_res_c[,loop.h], 
     xlab = expression(X_(perp)), ylab = expression(hat(epsilon)),
     main = paste0("p-values of Ljung-Box Q statsitic with h=",h.seq_Q[loop.h])
)
abline(a=0,b=1,col=2)


df_tmp <- data.frame(x = -log(Q_p_Xperp_c[,loop.h]),
                     y = -log(Q_p_res_c[,loop.h])
)
df_tmp <- data.frame(x = (Q_p_Xperp_c[,loop.h]),
                     y = (Q_p_res_c[,loop.h])
)

ggplot(df_tmp, aes(x, y)) +
  geom_point(alpha = 0.7, size = 2, col="blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1, linetype = "solid")  + 
  labs(
    x = expression(hat(X)["\u22A5"]),          # ⊥ as Unicode subscript
    y = expression(hat(epsilon)),          # \widehat{ε}
    title = paste0("p values of Ljung-Box Q statsitic with h=",h.seq_Q[loop.h])
  )  +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),  # axis label size
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10)                 # tick labels (optional)
  ) +
  geom_vline(xintercept = (0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = (0.05), linetype = "dashed", color = "black", linewidth = 0.8) + 
  scale_x_continuous(
    breaks = c(0, 0.05, 0.25, 0.5, 0.75),
    labels = c("0", "0.05", "0.25", "0.5", "0.75")
  ) + 
  scale_y_continuous(
    breaks = c(0, 0.05, 0.1, 0.2, 0.3),
    labels = c("0", "0.05", "0.1", "0.2", "0.3")
  ) 

ggplot(df_tmp, aes(x, y)) +
  geom_point(alpha = 0.7, size = 2, col="blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1, linetype = "solid")  + 
  labs(
    x = expression(p[X,i]),          # ⊥ as Unicode subscript
    y = expression(p[epsilon,i]),          # \widehat{ε}
  )  +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),  # axis label size
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10)                 # tick labels (optional)
  ) +
  geom_vline(xintercept = (0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = (0.05), linetype = "dashed", color = "black", linewidth = 0.8) + 
  scale_x_continuous(
    breaks = c(0, 0.05, 0.25, 0.5, 0.75),
    labels = c("0", "0.05", "0.25", "0.5", "0.75")
  ) + 
  scale_y_continuous(
    breaks = c(0, 0.05, 0.1, 0.2, 0.3),
    labels = c("0", "0.05", "0.1", "0.2", "0.3")
  ) 













































































































































# ---------- Only draw the non-covarage prob plot of epsShock same seed in the simulation ---------------- #
# take lag 1, 5, 20
rm(list = ls())

# --- only draw the plot for nominal level 0.05 ---- #
nominal = 0.05



# -- epsShock, with Lag = 1 -- # 
load("1121_c_AR1X&eps_epsShock_changeShock_L1_sameseed.Rdata")
# make plot about non-coverage prob #------------------

N = nrow(Z)

k = ncol(Z) + 1

ncp_epsShk_L1 = matrix(0, nrow = length(shock.seq), ncol = 2 )
rownames(ncp_epsShk_L1) = shock.seq
colnames(ncp_epsShk_L1) = c("HAC", "HAC HC1")

for (loop.S in 1:length(shock.seq)) {
  ncp_epsShk_L1[loop.S, "HAC"] = sum( abs(t_HAC[,loop.S]) > qnorm( (nominal/2), lower.tail = F) ) / run
}

for (loop.S in 1:length(shock.seq)) {
  ncp_epsShk_L1[loop.S, "HAC HC1"] = sum( abs( t_HAC[,loop.S] * sqrt( (N - k) / N ) ) > qnorm( (nominal/2), lower.tail = F) ) / run
}


Shock.seq_eps_L1 = shock.seq

rm(list = setdiff(ls(), 
                  c("nominal", "ncp_epsShk_L1", "Shock.seq_eps_L1", "N", "k")
)
)


# -- epsShock, with Lag = 5 -- # 
load("1121_c_AR1X&eps_epsShock_changeShock_L5_sameseed.Rdata")
# make plot about non-coverage prob #------------------

ncp_epsShk_L5 = matrix(0, nrow = length(shock.seq), ncol = 2 )
rownames(ncp_epsShk_L5) = shock.seq
colnames(ncp_epsShk_L5) = c("HAC", "HAC HC1")

for (loop.S in 1:length(shock.seq)) {
  ncp_epsShk_L5[loop.S, "HAC"] = sum( abs(t_HAC[,loop.S]) > qnorm( (nominal/2), lower.tail = F) ) / run
}

for (loop.S in 1:length(shock.seq)) {
  ncp_epsShk_L5[loop.S, "HAC HC1"] = sum( abs( t_HAC[,loop.S] * sqrt( (N - k) / N ) ) > qnorm( (nominal/2), lower.tail = F) ) / run
}


Shock.seq_eps_L5 = shock.seq

rm(list = setdiff(ls(), 
                  c("nominal", "ncp_epsShk_L1", "Shock.seq_eps_L1",
                    "ncp_epsShk_L5", "Shock.seq_eps_L5", "N", "k")
)
)











# -- epsShock, with Lag = 20 -- # 
load("1121_c_AR1X&eps_epsShock_changeShock_L20_sameseed.Rdata")
# make plot about non-coverage prob #------------------

ncp_epsShk_L20 = matrix(0, nrow = length(shock.seq), ncol = 2 )
rownames(ncp_epsShk_L20) = shock.seq
colnames(ncp_epsShk_L20) = c("HAC", "HAC HC1")

for (loop.S in 1:length(shock.seq)) {
  ncp_epsShk_L20[loop.S, "HAC"] = sum( abs(t_HAC[,loop.S]) > qnorm( (nominal/2), lower.tail = F) ) / run
}

for (loop.S in 1:length(shock.seq)) {
  ncp_epsShk_L20[loop.S, "HAC HC1"] = sum( abs( t_HAC[,loop.S] * sqrt( (N - k) / N ) ) > qnorm( (nominal/2), lower.tail = F) ) / run
}




Shock.seq_eps_L20 = shock.seq

rm(list = setdiff(ls(), 
                  c("nominal", "ncp_epsShk_L1", "Shock.seq_eps_L1",
                    "ncp_epsShk_L5", "Shock.seq_eps_L5",
                    "ncp_epsShk_L20", "Shock.seq_eps_L20", "N", "k"
                  )
)
)


ncp.frame_HAC = data.frame(
  # if you want to see HAC HC3, just change ncp column here and other parts do not need to be changed
  ncp = c( 
    ncp_epsShk_L1[,"HAC"], ncp_epsShk_L5[,"HAC"],    ncp_epsShk_L20[,"HAC"]
  ),
  
  Shock = c( 
    Shock.seq_eps_L1, Shock.seq_eps_L5,  Shock.seq_eps_L20 
  ),
  Lag = c(
    rep(1, times = length(Shock.seq_eps_L1)),
    rep(5, times = length(Shock.seq_eps_L5)),
    rep(20, times = length(Shock.seq_eps_L20))
  ),
  nominal = "nominal level: 0.05"
  )







ncp.frame_HACHC1 = data.frame(
  # if you want to see HAC HC3, just change ncp column here and other parts do not need to be changed
  ncp = c( 
    ncp_epsShk_L1[,"HAC HC1"], 
    ncp_epsShk_L5[,"HAC HC1"], 
    ncp_epsShk_L20[,"HAC HC1"]
  ),
  
  Shock = c( 
    Shock.seq_eps_L1, Shock.seq_eps_L5, Shock.seq_eps_L20 
  ),
  Lag = c(
    rep(1, times = length(Shock.seq_eps_L1)),
    rep(5, times = length(Shock.seq_eps_L5)),
    rep(20, times = length(Shock.seq_eps_L20))
  ),
  nominal = "nominal level: 0.05"
)




library(ggplot2)
library("ggpubr")
library(ggh4x)

ggplot(ncp.frame_HAC, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1)  +
  geom_pointpath(aes(color = factor(Lag),
                     linetype = factor(Lag) ),   # <- map linetype to Lag),
                 size = 2, linewidth = 0.8, shape = 19
  ) +
  scale_color_manual(values = c("orange", "brown", "magenta"  ))  +
  scale_linetype_manual(values = c("solid","solid","solid")) +  # 3 lags
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_y_continuous(limits = c(0.05, 0.11),
                     breaks = c(0.05, 0.07, 0.09)) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(legend.position = "none")




ggplot(ncp.frame_HACHC1, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(aes(color = factor(Lag),
                     linetype = factor(Lag) ),   # <- map linetype to Lag),
                 size = 2, linewidth = 0.8, shape = 19
  ) +
  scale_color_manual(values = c("orange", "brown", "magenta"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) +  # 3 lags
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_y_continuous(limits = c(0.05, 0.09),
                     breaks = c(0.05, 0.07, 0.09)) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(legend.position = "none")






ggplot(ncp.frame_HACHC1, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(aes(color = factor(Lag),
                     linetype = factor(Lag) ),   # <- map linetype to Lag),
                 size = 2, linewidth = 0.8, shape = 19
  ) +
  scale_color_manual(values = c("orange", "brown", "magenta"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) +  # 3 lags
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  scale_y_continuous(limits = c(0., 0.2),
                     breaks = c(0, 0.05, 0.1, 0.15)) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(legend.position = "none")


# ---------- non-coverage probability with larger legend -------------------- #
ncp.frame_HACHC1$Lag = paste0("Lag = ", ncp.frame_HACHC1$Lag)
ncp.frame_HACHC1$Lag = factor(ncp.frame_HACHC1$Lag, levels = c("Lag = 1", "Lag = 5", "Lag = 20" ))

ggplot(ncp.frame_HACHC1, aes(x = Shock, y = ncp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(aes(color = Lag,
                     linetype = Lag ),   # <- map linetype to Lag),
                 size = 2, linewidth = 0.8, shape = 19
  ) +
  scale_color_manual(values = c("orange", "brown", "magenta"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) +  # 3 lags
  geom_hline(yintercept = 0.05, color = "black", linetype = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = 2) + 
  scale_y_continuous(limits = c(0., 0.2),
                     breaks = c(0, 0.05, 0.1, 0.15)) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(
    ## Legend inside plot (bottom-right)
    legend.position = c(1, 0.6),
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


cp.frame_HACHC1 = ncp.frame_HACHC1
cp.frame_HACHC1$cp = 1 - cp.frame_HACHC1$ncp
cp.frame_HACHC1$nominal = "nominal level: 0.95"



# -------- draw the coverage probability ------------- #
ggplot(cp.frame_HACHC1, aes(x = Shock, y = cp)) +
  facet_wrap(~nominal, ncol = 1) +
  geom_pointpath(aes(color = Lag,
                     linetype = Lag ),   # <- map linetype to Lag),
                 size = 2, linewidth = 0.8, shape = 19
  ) +
  scale_color_manual(values = c("orange", "brown", "magenta"))  +
  scale_linetype_manual(values = c("solid","solid","solid")) +  # 3 lags
  geom_hline(yintercept = 0.95, color = "black", linetype = 2) +
  scale_y_continuous(limits = c(0.8, 0.96),
                     breaks = c(0.85, 0.9, 0.95, 1)) +
  labs(x = expression(shock~strength~kappa), y = "") +
  theme_bw() +
  theme(
    ## Legend inside plot (bottom-right)
    legend.position = c(0.96, 0.1),
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




