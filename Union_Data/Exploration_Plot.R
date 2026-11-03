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
load("Data_Setting.Rdata")

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

