##################################################################################################################
######################### No sample splitting, IV, or resampling, just analyze the original data #################
##################################################################################################################
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
    legend.key.size = unit(1.2, "lines"),                  # enlarge color lines
    legend.spacing.y = unit(0.3, "cm")                     # vertical spacing
  )

rm(union_top10_ave, union_top10_ave_gg)

# ---  5) Choose the response and predictor variables in Table III --- #
## dependent variable
# depvar = "gini" 
depvar = "top10"

# I will not consider the 1929 case in col (13) and (16)
if (depvar %in% c("lshare_nni", "top10")) {
  skill_var <- "ln_chsr_eu_ipolate_age1"
} else {
  skill_var <- "ln_chsr_eu_age1" # skill share variable
}



# --- 6) Choose variables --- # 
# --- 6.1) Full controls (right columns in Table III): income/tax/policy + industry shares (drop one as baseline) --- #
controls <- c("lnincpc", "share_return", "max_mw_all",
              "pollib_median", "state_top_mtr", "democrat_gov",  #
              "emp_share_mining", "emp_share_mining_intp",
              "emp_share_construct", "emp_share_construct_intp",
              "emp_share_manuf", "emp_share_manuf_intp",
              "emp_share_transport", "emp_share_transport_intp",
              "emp_share_trade", "emp_share_trade_intp",
              "emp_share_finance", "emp_share_finance_intp",
              "emp_share_services", "emp_share_services_intp",
              "emp_share_gov", "emp_share_gov_intp")

controls

# --- 6.2) fixed effects --- #
fe <- " | statefip + year^south"




# --- 6.3) Time window --- #
# I will not consider the 1929 case in col (13) and (16)
if (depvar %in% c("lshare_nni", "top10")) {
  min_year <- 1938L # it should be 1937, but after remove the missing values in the control variable, there is only one observation in 1937, so we remove this year.
} else {
  min_year <- 1940L # this is from "if year>=1940" in the boot_union_program.do file
}

max_year <- 2014L # This is from that pollib_median only shows non NA in the dataframe panel until 2014.

# --- 6.4) Trim dat0 into the selected time window --- #
dat0 <- dat0 %>% filter(year >= min_year, year <= max_year) # 3604 obs now.



# ---  7) Build the actual model frames (what feeds eq. 7 and eq. 6) --- #
# --- 7.1) Variables needed by each column --- #
# need9 is a character vector listing every variable your column (9) regression uses
need9  <- c(depvar, skill_var, "union_st_density") 
# need10 adds the extra controls used in column (10).
need10 <- c(need9, controls) 


# --- 7.2) Force a common sample across cols 9 & 10 --- # 
# complete.cases() returns a logical vector (TRUE/FALSE) marking rows with no NAs in any of those variables.
keep <- complete.cases(dat0[, need10, drop = FALSE])   # dat0[, need10, drop = FALSE] subsets dat0 by columns to just those needed for column (10). drop = FALSE keeps this as a data frame even if the subset ends up with a single column (avoids accidental vector conversion).
# dat is the data frame containing all the used variables without NA's.
dat  <- dat0[keep, , drop = FALSE]  # Filters dat0 down to only those rows with all required variables present.
nrow(dat) # 3550 state-year units remains for the following regression of eqn (7) and (6).





### --- 7.3) the final target dataset to be used in the regression --- #
dat_target = dat # 3550 state-year units



# --- 8) Eq. (6) without IV: depvar on union_st_density + skill_var (+ controls) | FE --- #
# --- Without control --- #
fml6_nc <- as.formula(paste0( 
  depvar, " ~ ", "union_st_density", " + ", skill_var,
  fe
))
fit_nc <- feols(fml6_nc, 
                data = dat_target, # If you want to use the non-bootstrapped data, set "data = dat,"
                cluster = ~ statefip)
summary(fit_nc) # union_st_density  has estimated coefficient: -2.35375, significant


# --- With control --- #
fml6_c <- as.formula(paste0(
  depvar, " ~ ", paste(c("union_st_density", skill_var, controls), collapse = " + "),
  fe
))
fit_c <- feols(fml6_c, 
               data = dat_target, # If you want to use the non-bootstrapped data, set "data = dat,"
               cluster = ~ statefip)
summary(fit_c) # union_st_density has estimated coefficient: -2.052551

fit_c$collin.var # the removed variables by collinearity

# --- 9) ACF and PACF --- #
# --- 9.0) Define some functions --- #

### Function to transform a data frame with (state, year) rows, to a wide data frame for a variable of interest ###
# Input: dat_sy is the data frame, such as dat_target, which we used in the feols; target_var is the variable of interest (residual) of the feols
# length(target_var) = nrow(dat_sy), each row in dat_state_year is (state year).
# Output: a data frame, whose first column is "year", and the rest columns are the target_var values for each state.
transwide = function(dat_sy, target_var){
  
  output <- dat_sy %>%
    mutate(target = target_var) %>%
    select(year, statefip, target) %>%
    arrange(year, statefip) %>%
    pivot_wider(names_from = statefip, values_from = target, values_fill = NA)
  
  return(output)
}

### Function to calculate PACF mat ###
## Input: a wide data frame, whose first column is "year", and the rest columns are residuals for each state
## Output: a matrix pacf for each state, a vector which are n_used's for the pacf across states to construct the dashed lines of white noise.
PACF_mat = function(data_wide, lag_max=60){
  
  state_cols <- setdiff(names(data_wide), "year")
  
  PACF_mat = matrix(0, nrow = (lag_max), ncol=length(state_cols) )
  
  n_used = c()
  
  for (loop.st in 1:length(state_cols)) {
    x <- data_wide[[state_cols[loop.st]]]
    x <- x - mean(x, na.rm = TRUE)      # center
    x <- x[is.finite(x)]                # drop NAs
    if (length(x) < 2) {
      PACF_mat[,loop.st] = rep(NA_real_, lag_max)
      n_used[loop.st] = NA
    } else {
      pacf.x = stats::pacf(x, lag.max = lag_max, plot = FALSE)
      v <- as.numeric(pacf.x$acf)
      n_used[loop.st] = as.numeric(pacf.x$n.used)
      if (length(v) < (lag_max)) {
        PACF_mat[,loop.st] = c(v, rep(NA_real_, ( lag_max - length(v) ) ) )
      } else {
        PACF_mat[,loop.st] = v
      }
    }
  }
  
  return(list(PACF_mat = PACF_mat, n_used = n_used))
}

### Function to calculate ACF mat (This function is defined above, but I need n_used, so re-define it.)###
## Input: a wide data frame, whose first column is "year", and the rest columns are residuals for each state
## Output: a matrix acf for each state, a vector which are n_used's for acf across states to construct the dashed lines of white noise.
ACF_mat = function(data_wide, lag_max=60){
  
  state_cols <- setdiff(names(data_wide), "year")
  
  ACF_mat = matrix(0, nrow = (lag_max + 1), ncol=length(state_cols) )
  n_used = c()
  
  for (loop.st in 1:length(state_cols)) {
    x <- data_wide[[state_cols[loop.st]]]
    x <- x - mean(x, na.rm = TRUE)      # center
    x <- x[is.finite(x)]                # drop NAs
    if (length(x) < 2) {
      ACF_mat[,loop.st] = rep(NA_real_, lag_max + 1)
      n_used[loop.st] = NA
    } else {
      acf.x = stats::acf(x, lag.max = lag_max, plot = FALSE)
      v <- as.numeric(acf.x$acf)
      n_used[loop.st] = acf.x$n.used
      if (length(v) < (lag_max + 1)) {
        ACF_mat[,loop.st] = c(v, rep(NA_real_, (lag_max + 1) - length(v)))
      } else {
        ACF_mat[,loop.st] = v
      }
    }
  }
  
  return(list(ACF_mat = ACF_mat, n_used = n_used))
}




























# lag_max = length(unique(dat_target$year))

lag_max = 40

# ---  10.1) Without Control --- #

# --- 10.1.1 ACF ---#
# ACF across states for fit_nc residuals 
res_nc_wide = transwide(dat_sy = dat_target, target_var = fit_nc$residuals)
res_nc_acf_mat = ACF_mat(data_wide = res_nc_wide, lag_max = lag_max)

plot(x = 0:(lag_max),  y=apply(res_nc_acf_mat$ACF_mat, 1,mean), 
     type = "h", xlab = "Lag", ylab = "Average ACF")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )

# ACF across states for X_perp
fml_aux_nc <- as.formula(paste0(
  "union_st_density ~ ", skill_var,             # same control as iv_nc
  fe                                   # same FE as iv_nc
))
aux_nc <- feols(fml_aux_nc, 
                data = dat_target, 
                cluster = ~ statefip)
Xperp_nc_wide = transwide(dat_sy = dat_target, target_var = aux_nc$residuals)
Xperp_nc_acf_mat = ACF_mat(data_wide = Xperp_nc_wide, lag_max = lag_max)

plot(x = 1:(lag_max),  y=apply(Xperp_nc_acf_mat$ACF_mat[-1,], 1,mean), 
     type = "h", xlab = "Lag", ylab = "Average ACF")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(Xperp_nc_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(Xperp_nc_acf_mat$n_used) ) ), lty = 3 )




off  <- 0.12

# ACF for each state:
for (loop.s in 1:(ncol(res_nc_acf_mat$ACF_mat)) ) {
  plot(x = ( (1:lag_max) - off ),  
       y=res_nc_acf_mat$ACF_mat[-1,loop.s], 
       type = "h", lwd = 3, lend = "butt",
       xlab = "Lag", ylab = "Average ACF",
       ylim = c(-0.3,1.1), 
       xlim = c((0-3*off), lag_max),
       xaxs = "i", 
       yaxs = "i",
       col  = "#1f77b4",
       main = paste0("State ", loop.s, " ACF w/o Control") )
  lines(x = ( (1:lag_max) + off ),
        y=Xperp_nc_acf_mat$ACF_mat[-1, loop.s], 
        type = "h", lwd = 3, lend = "butt", col = "#d62728")
  points(x=( (1:lag_max) - off ), 
         y=res_nc_acf_mat$ACF_mat[-1,loop.s], 
         pch = 16, col = "#1f77b4")
  points(x = ( (1:lag_max) + off ), 
         y=Xperp_nc_acf_mat$ACF_mat[-1,loop.s],
         pch = 16, col = "#d62728")
  abline(h = 0, lty = 2)
  abline(h = ( qnorm(0.975) / sqrt(res_nc_acf_mat$n_used[loop.s]  ) ) , lty = 3 )
  abline(h = -( qnorm(0.975) / sqrt(res_nc_acf_mat$n_used[loop.s] ) ), lty = 3 )
  legend("topright",
         legend = c("X_perp ACF", "epsilon_hat ACF"),
         col    = c("#d62728", "#1f77b4"),
         lwd    = 3, pch = 16, bty = "n", cex = 0.9)
  
}


# ACF AVERAGED ACROSS STATES
plot( x = ( (1:lag_max) - off ), 
      y = apply(res_nc_acf_mat$ACF_mat[-1,], 1,mean),
      type = "h", lwd = 3, lend = "butt",
      xlab = "Lag", ylab = "Average ACF",
      ylim = c(-0.3,0.8), 
      xlim = c((0-3*off), lag_max),
      xaxs = "i", 
      yaxs = "i",
      col  = "#1f77b4",
      main = paste0("ACF ave across 47 states without control for ", depvar)
)
lines(x = ( (1:lag_max) + off ),
      y=apply(Xperp_nc_acf_mat$ACF_mat[-1,], 1,mean), 
      type = "h", lwd = 3, lend = "butt", col = "#d62728")
points(x=( (1:lag_max) - off ), 
       apply(res_nc_acf_mat$ACF_mat[-1,], 1,mean), 
       pch = 16, col = "#1f77b4")
points(x = ( (1:lag_max) + off ), 
       apply(Xperp_nc_acf_mat$ACF_mat[-1,], 1,mean),
       pch = 16, col = "#d62728")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )
legend("topright",
       legend = c("X_perp ACF", "epsilon_hat ACF"),
       col    = c("#d62728", "#1f77b4"),
       lwd    = 3, pch = 16, bty = "n", cex = 0.9)





# abs(ACF) AVERAGED ACROSS STATES
plot( x = ( (1:lag_max) - off ), 
      y = apply(abs(res_nc_acf_mat$ACF_mat[-1,]), 1,mean),
      type = "h", lwd = 3, lend = "butt",
      xlab = "Lag", ylab = "Average ACF",
      ylim = c(-0.3,0.9), 
      xlim = c((0-3*off), lag_max),
      xaxs = "i", 
      yaxs = "i",
      col  = "#1f77b4",
      main = paste0("ACF ave across 47 states without control for ", depvar)
)
lines(x = ( (1:lag_max) + off ),
      y=apply(abs(Xperp_nc_acf_mat$ACF_mat[-1,]), 1,mean), 
      type = "h", lwd = 3, lend = "butt", col = "#d62728")
points(x=( (1:lag_max) - off ), 
       apply(abs(res_nc_acf_mat$ACF_mat[-1,]), 1,mean), 
       pch = 16, col = "#1f77b4")
points(x = ( (1:lag_max) + off ), 
       apply(abs(Xperp_nc_acf_mat$ACF_mat[-1,]), 1,mean),
       pch = 16, col = "#d62728")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )
legend("topright",
       legend = c("X_perp ACF", "epsilon_hat ACF"),
       col    = c("#d62728", "#1f77b4"),
       lwd    = 3, pch = 16, bty = "n", cex = 0.9)



# --- 10.1.2 PACF --- #
# PACF across states for fit_nc residuals 
res_nc_wide = transwide(dat_sy = dat_target, target_var = fit_nc$residuals)
res_nc_pacf_mat = PACF_mat(data_wide = res_nc_wide, lag_max = lag_max)

plot(x = 1:(lag_max),  y=apply(res_nc_pacf_mat$PACF_mat, 1,mean), 
     type = "h", xlab = "Lag", ylab = "Average PACF")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )

# PACF across states for X_perp
Xperp_nc_wide = transwide(dat_sy = dat_target, target_var = aux_nc$residuals)
Xperp_nc_pacf_mat = PACF_mat(data_wide = Xperp_nc_wide, lag_max = lag_max)

plot(x = 1:(lag_max),  y=apply(Xperp_nc_pacf_mat$PACF_mat, 1,mean), 
     type = "h", xlab = "Lag", ylab = "Average ACF")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(Xperp_nc_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(Xperp_nc_acf_mat$n_used) ) ), lty = 3 )




off  <- 0.12

# PACF for each state:
for (loop.s in 1:(ncol(res_nc_pacf_mat$PACF_mat)) ) {
  plot(x = ( (1:lag_max) - off ),  
       y=res_nc_pacf_mat$PACF_mat[,loop.s], 
       type = "h", lwd = 3, lend = "butt",
       xlab = "Lag", ylab = "Average ACF",
       ylim = c(-0.3,1.1), 
       xlim = c((0-3*off), lag_max),
       xaxs = "i", 
       yaxs = "i",
       col  = "#1f77b4",
       main = paste0("State ", loop.s, " PACF w/o Control") )
  lines(x = ( (1:lag_max) + off ),
        y=Xperp_nc_pacf_mat$PACF_mat[, loop.s], 
        type = "h", lwd = 3, lend = "butt", col = "#d62728")
  points(x=( (1:lag_max) - off ), 
         y=res_nc_pacf_mat$PACF_mat[,loop.s], 
         pch = 16, col = "#1f77b4")
  points(x = ( (1:lag_max) + off ), 
         y=Xperp_nc_pacf_mat$PACF_mat[,loop.s],
         pch = 16, col = "#d62728")
  abline(h = 0, lty = 2)
  abline(h = ( qnorm(0.975) / sqrt(res_nc_acf_mat$n_used[loop.s]  ) ) , lty = 3 )
  abline(h = -( qnorm(0.975) / sqrt(res_nc_acf_mat$n_used[loop.s] ) ), lty = 3 )
  legend("topright",
         legend = c("X_perp PACF", "epsilon_hat PACF"),
         col    = c("#d62728", "#1f77b4"),
         lwd    = 3, pch = 16, bty = "n", cex = 0.9)
  
}


# PACF AVERAGED ACROSS STATES
plot( x = ( (1:lag_max) - off ), 
      y = apply(res_nc_pacf_mat$PACF_mat, 1,mean),
      type = "h", lwd = 3, lend = "butt",
      xlab = "Lag", ylab = "Average ACF",
      ylim = c(-0.3,1.1), 
      xlim = c((0-3*off), lag_max),
      xaxs = "i", 
      yaxs = "i",
      col  = "#1f77b4",
      main = paste0("PACF ave across 47 states w/o control for ", depvar)
)
lines(x = ( (1:lag_max) + off ),
      y=apply(Xperp_nc_pacf_mat$PACF_mat, 1,mean), 
      type = "h", lwd = 3, lend = "butt", col = "#d62728")
points(x=( (1:lag_max) - off ), 
       apply(res_nc_pacf_mat$PACF_mat, 1,mean), 
       pch = 16, col = "#1f77b4")
points(x = ( (1:lag_max) + off ), 
       apply(Xperp_nc_pacf_mat$PACF_mat, 1,mean),
       pch = 16, col = "#d62728")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_nc_acf_mat$n_used) ) ), lty = 3 )
legend("topright",
       legend = c("X_perp PACF", "epsilon_hat PACF"),
       col    = c("#d62728", "#1f77b4"),
       lwd    = 3, pch = 16, bty = "n", cex = 0.9)



















# --- 10.2) With Control --- #


# --- 10.1.1 ACF ---#
# ACF across states for fit_nc residuals 
res_c_wide = transwide(dat_sy = dat_target, target_var = fit_c$residuals)
res_c_acf_mat = ACF_mat(data_wide = res_c_wide, lag_max = lag_max)

plot(x = 1:(lag_max),  y=apply(res_c_acf_mat$ACF_mat[-1,], 1,mean), 
     type = "h", xlab = "Lag", ylab = "Average ACF")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )

# ACF across states for X_perp
fml_aux_c <- as.formula(paste0(
  "union_st_density ~ ", paste(c(skill_var, controls), collapse = " + "),            # same control as iv_nc
  fe                                   # same FE as iv_nc
))
aux_c <- feols(fml_aux_c, 
               data = dat_target, 
               cluster = ~ statefip)
aux_c$collin.var

aux_c$collin.var == fit_c$collin.var # same as fit_c multicollinearity

Xperp_c_wide = transwide(dat_sy = dat_target, target_var = aux_c$residuals)
Xperp_c_acf_mat = ACF_mat(data_wide = Xperp_c_wide, lag_max = lag_max)

plot(x = 1:(lag_max),  y=apply(Xperp_c_acf_mat$ACF_mat[-1,], 1,mean), 
     type = "h", xlab = "Lag", ylab = "Average ACF")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(Xperp_c_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(Xperp_c_acf_mat$n_used) ) ), lty = 3 )



off  <- 0.12

# ACF for each state:
for (loop.s in 1:(ncol(res_c_acf_mat$ACF_mat)) ) {
  plot(x = ( (1:lag_max) - off ),  
       y=res_c_acf_mat$ACF_mat[-1,loop.s], 
       type = "h", lwd = 3, lend = "butt",
       xlab = "Lag", ylab = "Average ACF",
       ylim = c(-0.5,0.9), 
       xlim = c((0-3*off), lag_max),
       xaxs = "i", 
       yaxs = "i",
       col  = "#1f77b4",
       main = paste0("State ", loop.s, " ACF w/ Control") )
  lines(x = ( (1:lag_max) + off ),
        y=Xperp_c_acf_mat$ACF_mat[-1, loop.s], 
        type = "h", lwd = 3, lend = "butt", col = "#d62728")
  points(x=( (1:lag_max) - off ), 
         y=res_c_acf_mat$ACF_mat[-1,loop.s], 
         pch = 16, col = "#1f77b4")
  points(x = ( (1:lag_max) + off ), 
         y=Xperp_c_acf_mat$ACF_mat[-1,loop.s],
         pch = 16, col = "#d62728")
  abline(h = 0, lty = 2)
  abline(h = ( qnorm(0.975) / sqrt(res_c_acf_mat$n_used[loop.s]  ) ) , lty = 3 )
  abline(h = -( qnorm(0.975) / sqrt(res_c_acf_mat$n_used[loop.s] ) ), lty = 3 )
  legend("topright",
         legend = c("X_perp ACF", "epsilon_hat ACF"),
         col    = c("#d62728", "#1f77b4"),
         lwd    = 3, pch = 16, bty = "n", cex = 0.9)
  
}




# some funciton of (ACF) AVERAGED ACROSS STATES

# Averaged ACF across states
plot( x = ( (1:lag_max) - off ), 
      y = sqrt( apply( ( (res_c_acf_mat$ACF_mat[-1,])^2 ), 1,mean ) ),
      type = "h", lwd = 3, lend = "butt",
      xlab = "Lag", ylab = " ",
      ylim = c(-0.1,0.9), 
      xlim = c((0-3*off), lag_max+2*off),
      xaxs = "i", 
      yaxs = "i",
      col  = "#1f77b4",
      main = paste0("sqrt[ ave (ACF^2) across 47 states] w/ control for ", depvar)
)
lines(x = ( (1:lag_max) + off ),
      y = sqrt( apply( ( (Xperp_c_acf_mat$ACF_mat[-1,])^2), 1, mean) ), 
      type = "h", lwd = 3, lend = "butt", col = "#d62728")
points(x=( (1:lag_max) - off ), 
       sqrt( apply( ( ( res_c_acf_mat$ACF_mat[-1,])^2 ), 1,mean) ), 
       pch = 16, col = "#1f77b4")
points(x = ( (1:lag_max) + off ), 
       sqrt( apply( ( (Xperp_c_acf_mat$ACF_mat[-1,])^2 ), 1,mean ) ),
       pch = 16, col = "#d62728")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
legend("topright",
       legend = c("X_perp ACF", "epsilon_hat ACF"),
       col    = c("#d62728", "#1f77b4"),
       lwd    = 3, pch = 16, bty = "n", cex = 0.9)


# abs(ACF) AVERAGED ACROSS STATES
plot( x = ( (1:lag_max) - off ), 
      y =  apply( ( abs(res_c_acf_mat$ACF_mat[-1,]) ), 1,mean ) ,
      type = "h", lwd = 3, lend = "butt",
      xlab = "Lag", ylab = " ",
      ylim = c(-0.1,0.9), 
      xlim = c((0-3*off), lag_max+2*off),
      xaxs = "i", 
      yaxs = "i",
      col  = "#1f77b4",
      main = paste0("abs(ACF) ave across 47 states w/ control")
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
       legend = c("X_perp ACF", "epsilon_hat ACF"),
       col    = c("#d62728", "#1f77b4"),
       lwd    = 3, pch = 16, bty = "n", cex = 0.9)




# ACF AVERAGED ACROSS STATES
plot( x = ( (1:lag_max) - off ), 
      y = apply(res_c_acf_mat$ACF_mat[-1,], 1,mean),
      type = "h", lwd = 3, lend = "butt",
      xlab = "Lag", ylab = "Average ACF",
      ylim = c(-0.3,0.9), 
      xlim = c((0-3*off), lag_max),
      xaxs = "i", 
      yaxs = "i",
      col  = "#1f77b4",
      main = paste0("ACF ave across 47 states w/ control for ", depvar)
)
lines(x = ( (1:lag_max) + off ),
      y=apply(Xperp_c_acf_mat$ACF_mat[-1,], 1,mean), 
      type = "h", lwd = 3, lend = "butt", col = "#d62728")
points(x=( (1:lag_max) - off ), 
       apply(res_c_acf_mat$ACF_mat[-1,], 1,mean), 
       pch = 16, col = "#1f77b4")
points(x = ( (1:lag_max) + off ), 
       apply(Xperp_c_acf_mat$ACF_mat[-1,], 1,mean),
       pch = 16, col = "#d62728")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
legend("topright",
       legend = c("X_perp ACF", "epsilon_hat ACF"),
       col    = c("#d62728", "#1f77b4"),
       lwd    = 3, pch = 16, bty = "n", cex = 0.9)





# --- 10.1.2 PACF --- #
# PACF across states for fit_nc residuals 
res_c_wide = transwide(dat_sy = dat_target, target_var = fit_c$residuals)
res_c_pacf_mat = PACF_mat(data_wide = res_c_wide, lag_max = lag_max)

plot(x = 1:(lag_max),  y=apply(res_c_pacf_mat$PACF_mat, 1,mean), 
     type = "h", xlab = "Lag", ylab = "Average PACF")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )

# PACF across states for X_perp
Xperp_c_wide = transwide(dat_sy = dat_target, target_var = aux_c$residuals)
Xperp_c_pacf_mat = PACF_mat(data_wide = Xperp_c_wide, lag_max = lag_max)

plot(x = 1:(lag_max),  y=apply(Xperp_c_pacf_mat$PACF_mat, 1,mean), 
     type = "h", xlab = "Lag", ylab = "Average ACF")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(Xperp_c_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(Xperp_c_acf_mat$n_used) ) ), lty = 3 )


off  <- 0.12

# PACF for each state:
for (loop.s in 1:(ncol(res_c_pacf_mat$PACF_mat)) ) {
  plot(x = ( (1:lag_max) - off ),  
       y=res_c_pacf_mat$PACF_mat[,loop.s], 
       type = "h", lwd = 3, lend = "butt",
       xlab = "Lag", ylab = "Average PACF",
       ylim = c(-0.3,1.1), 
       xlim = c((0-3*off), lag_max),
       xaxs = "i", 
       yaxs = "i",
       col  = "#1f77b4",
       main = paste0("State ", loop.s, " PACF w/ Control") )
  lines(x = ( (1:lag_max) + off ),
        y=Xperp_c_pacf_mat$PACF_mat[, loop.s], 
        type = "h", lwd = 3, lend = "butt", col = "#d62728")
  points(x=( (1:lag_max) - off ), 
         y=res_c_pacf_mat$PACF_mat[,loop.s], 
         pch = 16, col = "#1f77b4")
  points(x = ( (1:lag_max) + off ), 
         y=Xperp_c_pacf_mat$PACF_mat[,loop.s],
         pch = 16, col = "#d62728")
  abline(h = 0, lty = 2)
  abline(h = ( qnorm(0.975) / sqrt(res_c_acf_mat$n_used[loop.s]  ) ) , lty = 3 )
  abline(h = -( qnorm(0.975) / sqrt(res_c_acf_mat$n_used[loop.s] ) ), lty = 3 )
  legend("topright",
         legend = c("X_perp PACF", "epsilon_hat PACF"),
         col    = c("#d62728", "#1f77b4"),
         lwd    = 3, pch = 16, bty = "n", cex = 0.9)
  
}


# PACF AVERAGED ACROSS STATES
plot( x = ( (1:lag_max) - off ), 
      y = apply(res_c_pacf_mat$PACF_mat, 1,mean),
      type = "h", lwd = 3, lend = "butt",
      xlab = "Lag", ylab = "Average ACF",
      ylim = c(-0.3,1.1), 
      xlim = c((0-3*off), lag_max),
      xaxs = "i", 
      yaxs = "i",
      col  = "#1f77b4",
      main = paste0("PACF ave across 47 states w/ control for ", depvar)
)
lines(x = ( (1:lag_max) + off ),
      y=apply(Xperp_c_pacf_mat$PACF_mat, 1,mean), 
      type = "h", lwd = 3, lend = "butt", col = "#d62728")
points(x=( (1:lag_max) - off ), 
       apply(res_c_pacf_mat$PACF_mat, 1,mean), 
       pch = 16, col = "#1f77b4")
points(x = ( (1:lag_max) + off ), 
       apply(Xperp_c_pacf_mat$PACF_mat, 1,mean),
       pch = 16, col = "#d62728")
abline(h = 0, lty = 2)
abline(h = ( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
abline(h = -( qnorm(0.975) / sqrt(mean(res_c_acf_mat$n_used) ) ), lty = 3 )
legend("topright",
       legend = c("X_perp PACF", "epsilon_hat PACF"),
       col    = c("#d62728", "#1f77b4"),
       lwd    = 3, pch = 16, bty = "n", cex = 0.9)



# 18 remained controls
controls.remain = setdiff(controls, fit_c$collin.var)

boxplot(dat_target[,c(controls.remain, skill_var)])

# additonally, besides the 18 remained controls, we need to consider skill_share
skill_var
hist(dat_target[[skill_var]])
summary(dat_target$ln_chsr_eu_ipolate_age1) # mean = -0.62323, sd = 0.7501657
sd((dat_target$ln_chsr_eu_ipolate_age1))

apply(dat_target[,controls.remain], 2, mean)
apply(dat_target[,controls.remain], 2, sd)
sapply(dat_target[,controls.remain], summary)


hist(dat_target[["lnincpc"]]) # kind of symmetric
mean(dat_target[["lnincpc"]]) # 1.546928
sd(dat_target[["lnincpc"]]) # 1.315307


hist(dat_target[["share_return"]]) # kind of symmetric
mean(dat_target[["share_return"]]) # 0.8598822
sd(dat_target[["share_return"]]) # 0.205121





hist(dat_target[["max_mw_all"]]) # right skewed
hist(log(dat_target[["max_mw_all"]])) # flat, like uniform
mean(log(dat_target[["max_mw_all"]])) # 0.6523264
sd(log(dat_target[["max_mw_all"]])) # 1.018366



hist(dat_target[["pollib_median"]]) # symmetric
mean(dat_target[["pollib_median"]]) # -0.01295696
sd(dat_target[["pollib_median"]]) # 0.9774815


hist(dat_target[["state_top_mtr"]]) # right skewed with negative values, so no log
mean((dat_target[["state_top_mtr"]])) # 5.551241
sd((dat_target[["state_top_mtr"]])) # 4.528699



hist(dat_target[["democrat_gov"]]) # binary
mean((dat_target[["democrat_gov"]])) # 0.5497043




hist(dat_target[["emp_share_mining"]]) # right skewed with negative values, so no log
mean((dat_target[["emp_share_mining"]])) # 0.01776567
sd((dat_target[["emp_share_mining"]])) # 0.031457



hist(dat_target[["emp_share_mining_intp"]], breaks = 10) # binary
mean((dat_target[["emp_share_mining_intp"]])) # 0.8203323


hist(dat_target[["emp_share_construct"]]) # symmetric
mean((dat_target[["emp_share_construct"]])) # 0.05266156
sd((dat_target[["emp_share_construct"]])) # 0.01587423



hist(dat_target[["emp_share_construct_intp"]]) # binary
mean((dat_target[["emp_share_construct_intp"]])) # 0.8251197


hist(dat_target[["emp_share_manuf"]]) # right skewed, all positive
mean((dat_target[["emp_share_manuf"]])) # 0.2188569
sd((dat_target[["emp_share_manuf"]])) # 0.1171774


hist(dat_target[["emp_share_manuf_intp"]]) # binary
mean((dat_target[["emp_share_manuf_intp"]])) # 0.8268093



hist(dat_target[["emp_share_transport"]]) # right skewed, all positive
mean((dat_target[["emp_share_transport"]])) # 0.06734218
sd((dat_target[["emp_share_transport"]])) # 0.02619768



hist(dat_target[["emp_share_trade"]]) # symmetric
mean((dat_target[["emp_share_trade"]])) # 0.2217757
sd((dat_target[["emp_share_trade"]])) # 0.02900262


hist(dat_target[["emp_share_finance"]]) # slightly right skewed
mean((dat_target[["emp_share_finance"]])) # 0.04724264
sd((dat_target[["emp_share_finance"]])) # 0.01461992


hist(dat_target[["emp_share_services"]]) # right skewed, all positive
mean((dat_target[["emp_share_services"]])) # 0.1987146
sd((dat_target[["emp_share_services"]])) # 0.08344819

hist(dat_target[["emp_share_services_intp"]]) # binary
mean((dat_target[["emp_share_services_intp"]])) # 0.818361

hist(dat_target[["emp_share_gov"]]) # symmetric
mean((dat_target[["emp_share_gov"]])) # 0.1755461
sd((dat_target[["emp_share_gov"]])) # 0.03978153










# ------------- Using lm() to Recover fitted values ---------------- #

# --- if fe <- " | statefip + year^south", it is easy to recover it using lm() --- #

# --- With controls as an example --- #

fe  

fml6_c 
# recall fit_c
fit_c

fit_c$coefficients
summary(fit_c)

# fixed effect matrix
M_fe <- model.matrix(
  ~ 0 +
    factor(statefip) +
    factor(year) + factor(south) + factor(year):factor(south),
  data = dat_target
)

colnames(M_fe)

nrow(dat_target) # 3550
length(unique(dat_target$statefip)) # 47 states
length(unique(dat_target$year)) # 77 years

dim(M_fe) # 3550  200, where nrow(X_fe) = nrow(dat_target); ncol(X_fe) = 47(states) + 76(number of years-1) + 1(south) + 76(for south:year, which is number of years-1)

rhs <- c("union_st_density", "ln_chsr_eu_ipolate_age1", controls)

M_main = dat_target[,c("union_st_density",skill_var, controls)]

head(M_main)

y <- dat_target[[depvar]] 

# 4) Fit OLS with explicit FE columns
M_full <- cbind(M_main, M_fe)

M_full = as.matrix(M_full)

# remove those multicollinear columns
M_full = M_full[,
                !( colnames(M_full) 
                   %in% 
                     c("emp_share_transport_intp", "emp_share_trade_intp", "emp_share_finance_intp", "emp_share_gov_intp", "factor(south)1" )
                   )
                ]

min(eigen(t(M_full)%*%M_full)$values)

dim(M_full) # 3550  219

fit_c_lm = lm(y~M_full-1)

sum(is.na(fit_c_lm))

fit_c$coefficients
fit_c_lm$coefficients # coef are the same

max(abs(fit_c$residuals - fit_c_lm$residuals)) # same residuals
max(abs(fit_c$fitted.values - fit_c_lm$fitted.values)) # same fitted values




# construct (\beta, \gamma)
fit_c_lm_coef = fit_c_lm$coefficients
sum(is.na(fit_c_lm_coef) )  # construct the coef vector fitting y~M_full, and fill those NA's with zeros.

max(abs(fit_c_lm$fitted.values - M_full  %*% fit_c_lm_coef )) # the new constructed coef filled with zeros for NA's is valid
max(abs(y - fit_c_lm$fitted.values - fit_c_lm$residuals))


# --- check X_perp --- #
# recall aux_c
fml_aux_c 
aux_c 
aux_c$collin.var

colnames(M_full)
aux_c_lm = lm(M_full[,1 ]~M_full[,-1]-1)

aux_c$coefficients
aux_c_lm$coefficients # coeffficients are the same

sum(is.na(aux_c_lm$coefficients)) # no missing (no collinearity)

max(abs(aux_c$residuals - aux_c_lm$residuals))
max(abs(aux_c$fitted.values - aux_c_lm$fitted.values))


# construct \gamma
aux_c_lm_coef = aux_c_lm$coefficients
sum(is.na(aux_c_lm_coef) )

max(abs(aux_c_lm$fitted.values - M_full[,-1]  %*% aux_c_lm_coef )) # the new constructed coef filled with zeros for NA's is valid
max(abs(M_full[,1] - aux_c_lm$fitted.values - aux_c_lm$residuals))






# these col-linear variables are different
aux_c$collin.var
colnames(M_full)[is.na(aux_c_lm$coefficients)]
fit_c$collin.var
colnames(M_full)[is.na(fit_c_lm$coefficients)]

# check some col-linear variables in the data 
plot(dat_target$emp_share_transport, dat_target$emp_share_transport_intp)


# if fe <- " | statefip[year + I((year^2)/1000)] + year^south", it is not easy to recover it with lm()














#------ Ljung-Box test for each state ---------#:
h.seq_Q = c(10:20)

# --- with controls --- #

# check residuals
Q_res_c = matrix(0, nrow = ( ncol(res_c_wide )-1 ), ncol = length(h.seq_Q)) # store Ljung-Box test statistics, 47 rows for states, each column represent the Qstat at lag h
Q_p_res_c = matrix(0, nrow = ( ncol(res_c_wide )-1 ), ncol = length(h.seq_Q)) # store Ljung-Box p values, 47 rows for states, each column represent the Qstat at lag h
rownames(Q_res_c) = colnames(res_c_wide[,-1])
for (loop.s in 1:( ncol(res_c_wide )-1 ) ) {
  v = res_c_wide[,(loop.s+1)]
  v <- v[!is.na(v)]
  for (loop.h in 1:length(h.seq_Q)) {
    h = h.seq_Q[loop.h]
    LBt = Box.test(v, lag = h, type = "Ljung-Box")  
    Q_res_c[loop.s, loop.h] = LBt$statistic
    Q_p_res_c[loop.s, loop.h] = LBt$p.value
  }
}
rm(loop.s,loop.h, h,v, LBt)


# check Xperp
Q_Xperp_c = matrix(0, nrow = ( ncol(Xperp_c_wide )-1 ), ncol = length(h.seq_Q)) # store Ljung-Box test statistics, 47 rows for states, each column represent the Qstat at lag h
Q_p_Xperp_c = matrix(0, nrow = ( ncol(Xperp_c_wide )-1 ), ncol = length(h.seq_Q)) # store Ljung-Box p values, 47 rows for states, each column represent the Qstat at lag h
rownames(Q_Xperp_c) = colnames(Xperp_c_wide[,-1])
rownames(Q_p_Xperp_c) = colnames(Xperp_c_wide[,-1])
for (loop.s in 1:( ncol(Xperp_c_wide )-1 ) ) {
  v = Xperp_c_wide[,(loop.s+1)]
  v <- v[!is.na(v)]
  for (loop.h in 1:length(h.seq_Q)) {
    h = h.seq_Q[loop.h]
    LBt = Box.test(v, lag = h, type = "Ljung-Box")  
    Q_Xperp_c[loop.s, loop.h] = LBt$statistic
    Q_p_Xperp_c[loop.s, loop.h] = LBt$p.value
  }
}
rm(loop.s,loop.h, h,v, LBt)

Q_p_Xperp_c

boxplot(cbind(Q_Xperp_c, Q_res_c))

for (loop.h in 1:length(h.seq_Q)) {
  h = h.seq_Q[loop.h]
  plot(x=Q_Xperp_c[,loop.h], y = Q_res_c[,loop.h], 
       xlab = expression(X_(perp)), ylab = expression(hat(epsilon)),
       main = paste0("Ljung-Box Q statsitic with lag h=",h)
  )
  abline(a=0,b=1,col=2)
}



































year = dat_target$year
statefip = dat_target$statefip

##################################################################
########### Some preparation for the simulation ##################
##################################################################
library(zoo) # impute missing values with function na.approx()

rm(list = setdiff(ls(), c(
  "M_full", "y", "fit_c_lm", "fit_c_lm_coef", 
  "aux_c_lm", "aux_c_lm_coef", "Q_res_c", "Q_Xperp_c", "Q_p_res_c", "Q_p_Xperp_c",
  "year", "statefip", "res_c_wide", "Xperp_c_wide",
  "dat_target", "res_c_acf_mat", "res_c_pacf_mat", "Xperp_c_acf_mat", "Xperp_c_pacf_mat"
)))


Z = M_full[,-1]
X_perp_ori = aux_c_lm$residuals
eps_perp_ori = fit_c_lm$residuals

eta = aux_c_lm_coef
bet = fit_c_lm_coef[1]
gamm = fit_c_lm_coef[-1]


rm(list = setdiff(ls(), c(
  "Z", "X_perp_ori", "eps_perp_ori", "eta", 
  "bet", "gamm", "Q_res_c", "Q_Xperp_c", "Q_p_res_c", "Q_p_Xperp_c",
  "year", "statefip", "res_c_wide", "Xperp_c_wide",
  "dat_target", "res_c_acf_mat", "res_c_pacf_mat", "Xperp_c_acf_mat", "Xperp_c_pacf_mat"
)))


res_c_wide


length(unique(year))
range(year)


















# Although I have res_c_wide, which organize the residuals into columns of states. But I still want to do it by myself
eps_perp_ori_wide = matrix(NA, nrow = nrow(res_c_wide), ncol=length(unique(statefip)))

rownames(eps_perp_ori_wide) = unique(year)
colnames(eps_perp_ori_wide) = unique(statefip)
for (loop.s in unique(statefip) ) {
  for (loop.y in unique(year) ) {
    tmp = eps_perp_ori[(year==loop.y)&(statefip==loop.s)]
    if (length(tmp)>0) {
      eps_perp_ori_wide[as.character(loop.y), as.character(loop.s)] = tmp
    }
  }
}
rm(loop.s, loop.y, tmp)

eps_perp_ori_wide

eps_perp_ori_wide[,"27"] - res_c_wide[,"27"] # these two matrices are the same.

rm(res_c_wide)

# impute the missing values
tmp_imp <- na.approx(eps_perp_ori_wide, na.rm = FALSE, rule=2)          # interpolate internal gaps, rule=2 will fill the end NA entries

rownames(tmp_imp) <- rownames(eps_perp_ori_wide)

colnames(tmp_imp) == colnames(eps_perp_ori_wide)

tmp_imp == eps_perp_ori_wide

tmp_imp[,"32"]
eps_perp_ori_wide[,"32"]
plot(tmp_imp[,"32"])

eps_perp_ori_wide = tmp_imp
rm(tmp_imp)




























# Although I have Xperp_c_wide, which organize the residuals into columns of states. But I still want to do it by myself
Xperp_ori_wide = matrix(NA, nrow = nrow(Xperp_c_wide), ncol=length(unique(statefip)))

rownames(Xperp_ori_wide) = unique(year)
colnames(Xperp_ori_wide) = unique(statefip)
for (loop.s in unique(statefip) ) {
  for (loop.y in unique(year) ) {
    tmp = X_perp_ori[(year==loop.y)&(statefip==loop.s)]
    if (length(tmp)>0) {
      Xperp_ori_wide[as.character(loop.y), as.character(loop.s)] = tmp
    }
    
  }
}
rm(loop.s, loop.y, tmp)

Xperp_ori_wide - Xperp_c_wide[,-1]

Xperp_ori_wide[,"1"] - Xperp_c_wide[,"1"] #
Xperp_ori_wide[,"4"] - Xperp_c_wide[,"4"] #
Xperp_ori_wide[,"27"] - Xperp_c_wide[,"27"] # these two matrices are the same.

rm(Xperp_c_wide)

# impute the missing values
tmp_imp <- na.approx(Xperp_ori_wide, na.rm = FALSE, rule=2)          # interpolate internal gaps, rule=2 will fill the end NA entries

rownames(tmp_imp) <- rownames(Xperp_ori_wide)

colnames(tmp_imp) == colnames(Xperp_ori_wide)

tmp_imp == Xperp_ori_wide

tmp_imp[,"32"]
Xperp_ori_wide[,"32"]
plot(tmp_imp[,"32"])

Xperp_ori_wide = tmp_imp
rm(tmp_imp)


Xperp_ori_wide == eps_perp_ori_wide



















rm(list = setdiff(ls(), c("Z",
  "X_perp_ori", "eps_perp_ori", "eta", 
  "bet", "gamm", "Q_res_c", "Q_Xperp_c", "Q_p_res_c", "Q_p_Xperp_c",
  "year", "statefip", "eps_perp_ori_wide", "Xperp_ori_wide",
  "dat_target", "res_c_acf_mat", "res_c_pacf_mat", "Xperp_c_acf_mat", "Xperp_c_pacf_mat"
)))


dim(Z) # 3550  218







# ---- Impute Z ----- #
# --- Since I have already imputed Xperp_ori_wide and eps_perp_ori_wide, I also need to impute Z --- #
# Z columns ncol(Z) =  1(skill_share)+6+12(emp_xxx)+47(state)+76(years)+76(south*year)
# So columns of Z has two parts: the first 1+6+12 columns, which can be imputed by na.approx in zoo.
# the other fe part from 13:218, should be reproduce by the same way as M_fe above


unique(year) - sort(unique(year))
unique(statefip) - sort(unique(statefip))

# first impute year
year_ipt = rep(unique(year), times = length( unique(statefip) ) )
# then we impute the fe part of Z
statefip_ipt = rep(unique(statefip), each = length(unique(year)) )
south_states <- c(10,11,12,13,24,37,51,54, 1, 21,28,47,5, 22,40,48)
south_ipt = 1 * (statefip_ipt %in% south_states)

year_state_df = data.frame(year = year_ipt, state = statefip_ipt, south = south_ipt)


Z_fe <- model.matrix(
  ~ 0 +
    factor(state) +
    factor(year) + factor(south) + factor(year):factor(south),
  data = year_state_df
)


colnames(Z_fe)

# remove south due to collinearity

Z_fe = Z_fe[, !( colnames(Z_fe)%in%c("factor(south)1") ) ]

dim(Z_fe) # 3619  199, where nrow(Z_fe) = 77 years * 47 states; ncol(X_fe) = 47(states) + 76(number of years-1) +  76(for south:year, which is number of years-1)

rm(year_state_df, south_ipt)

Z_fe = cbind(Z_fe, year_ipt, statefip_ipt) 

Z_fe[1:80,c(1,2, 200,201)] # currently, the first 77 rows in Z_fe is 1938-2014 at statefip 1, the second 77 rows are those 77 years for statefip 2,...
# This satisfies Prof. Lei's idea for the weight matrix diag(W1,...,W47)
######################################


# remove the year and state columns in Z_fe
Z_fe = Z_fe[,!(colnames(Z_fe) %in% c("year_ipt", "statefip_ipt"))]
dim(Z_fe) # 3619  199








# impute the first part of Z
Z_1.old = Z[,1:(1+6+12)]
dim(Z_1.old) # 3550  19
Z_1.old = cbind(Z_1.old, year, statefip)



Z_1.new = matrix(NA, nrow = length(year_ipt), ncol = (1+6+12) )
dim(Z_1.new) # 3619   19
Z_1.new = cbind(Z_1.new, year_ipt, statefip_ipt)
colnames(Z_1.new)[1:(1+6+12)] = colnames(Z_1.old)[1:(1+6+12)]
colnames(Z_1.new)


k1.old = paste(Z_1.old[,"year"],  Z_1.old[, "statefip"])
k1.new = paste(Z_1.new[, "year_ipt"],  Z_1.new[, "statefip_ipt"])

# 2) For each Z2 row, find the matching row in Z1 (NA if no match)
idx <- match(k1.old, k1.new)

# 3) Copy Z1’s first 10 columns into Z2’s first 10 columns (keeps NA where no match)
Z_1.new[idx, 1:(1+6+12)] <- Z_1.old[, 1:(1+6+12)]

rm(k1.old, k1.new, idx)

Z_1.new[1:77,] - Z_1.old[1:77,]
Z_1.new[78:154,] - Z_1.old[78:154,]

rm(Z_1.old)

sum(is.na(Z_1.new))
dim(Z_1.new)

tmp_imp = na.approx(Z_1.new, na.rm = FALSE, rule=2)

tmp_imp - Z_1.new

tmp_imp[78:154, ] - Z_1.new[78:154,]

sum(is.na(tmp_imp))

Z_1.new = tmp_imp
rm(tmp_imp)

Z_1.new[1:80, c(1,2,20,21)] # currently, the first 77 rows in Z_1.new is 1938-2014 at statefip 1, the second 77 rows are those 77 years for statefip 4,...
# This satisfies Prof. Lei's idea for the weight matrix diag(W1,...,W47)

# remove the year and state columns in Z_1.new
Z_1.new = Z_1.new[,!(colnames(Z_1.new) %in% c("year_ipt", "statefip_ipt"))]
dim(Z_1.new) # 3619  19



rm(Z)

# This new Z is imputed Z, with blocks of the same year
Z = cbind(Z_1.new, Z_fe)

dim(Z) # 3619  218

rm(Z_fe, Z_1.new)






length(gamm)
Z.gamm = Z%*%gamm

length(eta)
Z.eta = Z%*%eta


















rm(X_perp_ori, eps_perp_ori)











# --- These two vectors are loadings for some aggregate shock --- #
northeast_states = c(9, 23, 25, 33, 44, 50, 34, 36, 42)
west_states = c(4, 8, 16, 30, 32, 35, 49, 56, 2, 6, 15, 41, 53)

# the north east states and year >=2009 after financial crisis
loading1 = ( 1 * (statefip_ipt %in% northeast_states) ) * ( 1 * (year_ipt >= 2009) )
# the west states and year >=1983 when the internet was invented
loading2 = ( 1 * (statefip_ipt %in% west_states) ) * ( 1 * (year_ipt >= 1983) )

# Check whether loadings can be expressed by Z?
(lm(loading1~Z-1))$residuals # the residuals are not zeros, so loading 1 is not in the column space of Z

tmp = cbind(Z, loading1)
(lm(loading2~tmp-1))$residuals[1:10] # the residuals are not zeros, so loading 2 is not in the column space of cbind(Z, loading1)

rm(tmp)

# we need to redefine loading1 and loading2, because we have imputed Z.eta and Z.gamm, so they have 77*47 entries
# Yet now Z, and year, statefip only contains 3550 entries with several missing values, which does not match the imputed Z.gamm, Z.eta, Xperp_ori_wide, eps_perp_ori_wide
rm(loading1, loading2)

loading1_wide = (1*(unique(year_ipt) >=2009)) %*% t(1 * unique(statefip_ipt) %in% northeast_states)
rownames(loading1_wide) = unique(year_ipt)
colnames(loading1_wide) = unique(statefip_ipt) 
# loading1_wide has 77 rows as years, and 47 columns as states.

loading2_wide = (1*(unique(year_ipt) >=1983)) %*% t(1 * unique(statefip_ipt) %in% west_states)
rownames(loading2_wide) = unique(year_ipt)
colnames(loading2_wide) = unique(statefip_ipt)
# loading2_wide has 77 rows as years, and 47 columns as states.


rm(year, statefip)


save(list = c("dat_target", "eps_perp_ori_wide", 
              "loading1_wide", "loading2_wide",
              "Q_res_c", "Q_Xperp_c", "Q_p_res_c", "Q_p_Xperp_c",
              "res_c_acf_mat", "res_c_pacf_mat", "Xperp_c_acf_mat", "Xperp_c_pacf_mat",
              "Xperp_ori_wide", "Z", 
              "Z.eta", "Z.gamm",
              "bet", "eta", "gamm",
              "northeast_states", "south_states", "statefip_ipt", "west_states", "year_ipt"),
     file = "1108_PermSimu_Setting_c.Rdata")











































###########################################
### In the stored setting above, Z is arranged in this order:
### the first 77 rows in Z are for statefip 1, the second 77 rows of Z are for statefiip 4,....
### This satisfies Prof. Lei's suggestion that weight matrix diag(W1,...,W47)
###########################################




























































































































































































































































































































































































#########################################################################
########### Permute Simulation change Lag ############################### 
########### Permut eps_perp and Xperp together ##########################
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


# ----- start simulation --- #

run = 300

loading.strength = 1
Lag.seq = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,30, 35, 40, 45, 50, 55, 60, 65, 70, 75)

t_HAC = matrix(0, nrow = run, ncol = length(Lag.seq))
se2_HAC = matrix(0, nrow = run, ncol = length(Lag.seq))

t_HAC_HC2 = matrix(0, nrow = run, ncol = length(Lag.seq))
se2_HAC_HC2 = matrix(0, nrow = run, ncol = length(Lag.seq))

t_HAC_HC3 = matrix(0, nrow = run, ncol = length(Lag.seq))
se2_HAC_HC3 = matrix(0, nrow = run, ncol = length(Lag.seq))


loading1 = c(loading1_wide)
loading2 = c(loading2_wide)

n_state = ncol(eps_perp_ori_wide)
n_year = nrow(eps_perp_ori_wide)

state.unq = unique(statefip_ipt)

set.seed(715)
for (loop.L in 1:length(Lag.seq)) {
  
  Lag = Lag.seq[loop.L]
  Wgt.state = matrix(0, n_year, n_year)
  # Newey-West with finite sample adjustment
  for (Wgtrow in 1:n_year) {
    for (Wgtcol in 1:n_year) {
      Wgt.state[Wgtrow, Wgtcol] = max( 1 - abs(Wgtrow-Wgtcol)/(Lag+1), 0  ) 
    }
  }
  rm(Wgtrow, Wgtcol)
  
  
  Wgt_HAC = kronecker(diag(n_state), Wgt.state)
  #--------------------------------------------------
  
  for (loop.run in 1:run) {
    
    # permute epsperp #---------------------------------------
    
    new_state1 = sample(state.unq, n_state, replace = F)
    
    eps_perm_wide = eps_perp_ori_wide[,as.character(new_state1)]
    
    # dim(eps_perm_wide)
    # eps_perm_wide[,"1"] - eps_perp_ori_wide[,"1"]
    # eps_perm_wide - eps_perp_ori_wide
    
    eps = c(eps_perm_wide)
    
    
    # add aggregate shock to generate eps
    f1 = rnorm(1)
    f2 = rnorm(1)
    
    eps = eps + loading.strength * ( f1*loading1 + f2*loading2 )
    
    # plot(eps)
    
    
    # permute Xperp
    new_state2 = sample(state.unq, n_state, replace = F) # Note that eps and Xperp can not use the same permutation, which will induce endogeneity.
    
    Xperp_perm_wide = Xperp_ori_wide[,as.character(new_state2)]
    
    # dim(Xperp_perm_wide)
    # max(abs(Xperp_perm_wide[,"1"] - Xperp_ori_wide[,"1"]))
    # Xperp_perm_wide - Xperp_ori_wide
    
    Xperp = c(Xperp_perm_wide)
    
    # plot(Xperp)
    
    X  = Z.eta + Xperp 
    
    # plot(X)
    
    # generate y #----------------------------------------
    y = bet * X +  Z.gamm + eps
    
    
    M = cbind(X, Z)
    colnames(M)[1] = "X"
    
    fit.ols=lm(y~M-1)
    
    # HAC using X groups (too slow) #-------------------------------------------------------
    # vcov.hac = solve(t(M)%*%M)%*%(t(M)%*%((fit.ols$residuals%*%t(fit.ols$residuals))*Wgt_HAC)%*%M)%*%solve(t(M)%*%M)
    # t_HAC[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(vcov.hac[1,1])
    # se2_HAC[loop.run, loop.L] = vcov.hac[1,1]
    #---------------------------------------------------------------------------
    
    
    # --- HAC --- #
    fit.X = lm(X~Z-1)
    
    est.Xperp = fit.X$residuals
    
    est.eps.hat = fit.ols$residuals
    
    Hardmard = t(est.Xperp)%*% ( ( est.eps.hat%*%t(est.eps.hat) ) * Wgt_HAC ) %*% est.Xperp
    
    se2_denom = ( sum(est.Xperp*est.Xperp) )^2
    
    est.se2 = Hardmard / se2_denom
    
    # est.se2 - vcov.hac[1,1]
    
    se2_HAC[loop.run, loop.L] = est.se2
    t_HAC[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2)
    # -------------------------------------------------------------# 
    
    
    
    
    
    # --- HC2 adjustment --- #
    h_M = leverages_chol(M) # h_M is the diagonal entries of the projection matrix of M: (M%*%solve(t(M)%*%M)%*%t(M))
    
    est.eps.hat.HC2 = est.eps.hat / sqrt(1 - h_M) 
    
    Hardmard.HC2 = t(est.Xperp)%*% ( ( est.eps.hat.HC2%*%t(est.eps.hat.HC2) ) * Wgt_HAC ) %*% est.Xperp
    
    est.se2.HC2 = Hardmard.HC2 / se2_denom
    
    se2_HAC_HC2[loop.run, loop.L] = est.se2.HC2
    t_HAC_HC2[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2.HC2)
    # ---------------------------- #
    
    
    
    # --- HC3 adjustment --- #
    est.eps.hat.HC3 = est.eps.hat / (1 - h_M) 
    
    # est.eps.hat.HC2 - est.eps.hat.HC3
    
    Hardmard.HC3 = t(est.Xperp)%*% ( ( est.eps.hat.HC3%*%t(est.eps.hat.HC3) ) * Wgt_HAC ) %*% est.Xperp
    
    # max(abs(Hardmard.HC3 - Hardmard.HC2))
    
    est.se2.HC3 = Hardmard.HC3 / se2_denom
    
    # est.se2;est.se2.HC2;est.se2.HC3
    
    se2_HAC_HC3[loop.run, loop.L] = est.se2.HC3
    t_HAC_HC3[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2.HC3)
    # ---------------------------- #
    
    
    print(loop.run)
    
  }
  
  print(loop.L)
  
}

rm(loop.L, loop.run, fit.ols, fit.X, est.Xperp, est.eps.hat, Hardmard, Hardmard.HC2, Hardmard.HC3,
   est.se2, est.se2.HC2, est.se2.HC3, se2_denom)

save(list = c(
  "dat_target", "eps_perp_ori_wide",
  "loading1_wide", "loading2_wide",
  "Q_p_res_c", "Q_p_Xperp_c", "Q_res_c", "Q_Xperp_c",
  "res_c_acf_mat", "res_c_pacf_mat",
  "se2_HAC", "se2_HAC_HC2", "se2_HAC_HC3",
  "t_HAC", "t_HAC_HC2", "t_HAC_HC3",
  "Xperp_c_acf_mat", "Xperp_c_pacf_mat",
  "Z", "Z.eta", "Z.gamm",
  "bet", "eta", "gamm", "Lag.seq", "loading.strength", 
  "loading1", "loading2", "n_state", "n_year",
  "northeast_states", "run", "south_states",
  "state.unq", "statefip_ipt",
  "west_states", "year_ipt"
), 
file = "1108_c_PermuteX&eps_epsShock_changeLag.Rdata"
)

load("1108_c_PermuteX&eps_epsShock_changeLag.Rdata")

# make plot about non-coverage prob #------------------
# --- HAC --- #
qtl = qnorm(c(0.15, 0.1,0.05, 0.025, 0.005), lower.tail = F) # non-coverage prob of 0.3, 0.2, 0.1, 0.05, 0.01

ncp = matrix(0, nrow = length(Lag.seq), ncol = length(qtl))
rownames(ncp) = Lag.seq
colnames(ncp) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.L in 1:length(Lag.seq)) {
    
    ncp[loop.L, loop.qtl] = sum( abs(t_HAC[,loop.L]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.L, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=Lag.seq, y=ncp[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp[,loop.qtl]),
       xlab = "Lag",
       ylab = "Non-coverage probability HAC",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}





# --- HAC HC2 --- #
ncp_HC2 = matrix(0, nrow = length(Lag.seq), ncol = length(qtl))
rownames(ncp_HC2) = Lag.seq
colnames(ncp_HC2) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.L in 1:length(Lag.seq)) {
    
    ncp_HC2[loop.L, loop.qtl] = sum( abs(t_HAC_HC2[,loop.L]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.L, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=Lag.seq, y=ncp_HC2[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp_HC2[,loop.qtl]),
       xlab = "Lag",
       ylab = "Non-coverage probability HAC HC2",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}

rm(loop.qtl)



# --- HAC HC3 --- #
ncp_HC3 = matrix(0, nrow = length(Lag.seq), ncol = length(qtl))
rownames(ncp_HC3) = Lag.seq
colnames(ncp_HC3) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.L in 1:length(Lag.seq)) {
    
    ncp_HC3[loop.L, loop.qtl] = sum( abs(t_HAC_HC3[,loop.L]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.L, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=Lag.seq, y=ncp_HC3[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp_HC3[,loop.qtl]),
       xlab = "Lag",
       ylab = "Non-coverage probability HAC HC3",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}

rm(loop.qtl, loop.L)


# --- qqplot --- #
prob = seq(0.01,0.99, by=0.01)

for (loop.L in 1:length(Lag.seq)) {
  qt.emp = quantile(t_HAC[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("Lag = ", Lag.seq[loop.L]))
  abline(a=0, b=1, col=2, lwd=3)
}

for (loop.L in 1:length(Lag.seq)) {
  qt.emp = quantile(t_HAC_HC2[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("Lag = ", Lag.seq[loop.L]))
  abline(a=0, b=1, col=2, lwd=3)
}

for (loop.L in 1:length(Lag.seq)) {
  qt.emp = quantile(t_HAC_HC3[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("Lag = ", Lag.seq[loop.L]))
  abline(a=0, b=1, col=2, lwd=3)
}





# --- Check the size of standard error --- #
plot(x=Lag.seq,
     y=apply(sqrt(se2_HAC), 2, mean),
     type = "b",
     pch=19,
     xlab = "Lag",
     ylab = " ",
     main = "HAC Standard Error")

plot(x=Lag.seq,
     y=apply(sqrt(se2_HAC_HC2), 2, mean),
     type = "b",
     pch=19,
     xlab = "Lag",
     ylab = " ",
     main = "HAC Standard Error HC2")

plot(x=Lag.seq,
     y=apply(sqrt(se2_HAC_HC3), 2, mean),
     type = "b",
     pch=19,
     xlab = "Lag",
     ylab = " ",
     main = "HAC Standard Error HC3")



# It seems that the best Lag is around 20




















































#########################################################################
########### Permute Simulation change Lag ############################### 
########### Permut eps_perp and Xperp together ##########################
########### Add aggregate shock by loading1 and loading2 to Xperp ####### 
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


# ----- start simulation --- #

run = 300

loading.strength = 1
Lag.seq = c(0,5,10,15,20,25,30, 35, 40, 45, 50, 55, 60, 65, 70, 75)


t_HAC = matrix(0, nrow = run, ncol = length(Lag.seq))
se2_HAC = matrix(0, nrow = run, ncol = length(Lag.seq))

t_HAC_HC2 = matrix(0, nrow = run, ncol = length(Lag.seq))
se2_HAC_HC2 = matrix(0, nrow = run, ncol = length(Lag.seq))

t_HAC_HC3 = matrix(0, nrow = run, ncol = length(Lag.seq))
se2_HAC_HC3 = matrix(0, nrow = run, ncol = length(Lag.seq))


loading1 = c(loading1_wide)
loading2 = c(loading2_wide)

n_state = ncol(eps_perp_ori_wide)
n_year = nrow(eps_perp_ori_wide)

state.unq = unique(statefip_ipt)

set.seed(715)
for (loop.L in 1:length(Lag.seq)) {
  
  Lag = Lag.seq[loop.L]
  Wgt.state = matrix(0, n_year, n_year)
  # Newey-West with finite sample adjustment
  for (Wgtrow in 1:n_year) {
    for (Wgtcol in 1:n_year) {
      Wgt.state[Wgtrow, Wgtcol] = max( 1 - abs(Wgtrow-Wgtcol)/(Lag+1), 0  ) 
    }
  }
  rm(Wgtrow, Wgtcol)
  
  
  Wgt_HAC = kronecker(diag(n_state), Wgt.state)
  #--------------------------------------------------
  
  for (loop.run in 1:run) {
    
    # permute epsperp #---------------------------------------
    
    new_state1 = sample(state.unq, n_state, replace = F)
    
    eps_perm_wide = eps_perp_ori_wide[,as.character(new_state1)]
    
    # dim(eps_perm_wide)
    # eps_perm_wide[,"1"] - eps_perp_ori_wide[,"1"]
    # eps_perm_wide - eps_perp_ori_wide
    
    eps = c(eps_perm_wide)
    
    # plot(eps)
    
    
    # permute Xperp
    new_state2 = sample(state.unq, n_state, replace = F) # Note that eps and Xperp can not use the same permutation, which will induce endogeneity.
    
    Xperp_perm_wide = Xperp_ori_wide[,as.character(new_state2)]
    
    # dim(Xperp_perm_wide)
    # max(abs(Xperp_perm_wide[,"1"] - Xperp_ori_wide[,"1"]))
    # Xperp_perm_wide - Xperp_ori_wide
    
    Xperp = c(Xperp_perm_wide)
    
    # add aggregate shock to generate eps
    f1 = rnorm(1)
    f2 = rnorm(1)
    
    Xperp = Xperp + loading.strength * ( f1*loading1 + f2*loading2 )
    
    # plot(Xperp)
    
    X  = Z.eta + Xperp 
    
    # plot(X)
    
    # generate y #----------------------------------------
    y = bet * X +  Z.gamm + eps
    
    
    M = cbind(X, Z)
    colnames(M)[1] = "X"
    
    fit.ols=lm(y~M-1)
    
    # HAC using X groups (too slow) #-------------------------------------------------------
    # vcov.hac = solve(t(M)%*%M)%*%(t(M)%*%((fit.ols$residuals%*%t(fit.ols$residuals))*Wgt_HAC)%*%M)%*%solve(t(M)%*%M)
    # t_HAC[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(vcov.hac[1,1])
    # se2_HAC[loop.run, loop.L] = vcov.hac[1,1]
    #---------------------------------------------------------------------------
    
    
    # --- HAC --- #
    fit.X = lm(X~Z-1)
    
    est.Xperp = fit.X$residuals
    
    est.eps.hat = fit.ols$residuals
    
    Hardmard = t(est.Xperp)%*% ( ( est.eps.hat%*%t(est.eps.hat) ) * Wgt_HAC ) %*% est.Xperp
    
    se2_denom = ( sum(est.Xperp*est.Xperp) )^2
    
    est.se2 = Hardmard / se2_denom
    
    # est.se2 - vcov.hac[1,1]
    
    se2_HAC[loop.run, loop.L] = est.se2
    t_HAC[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2)
    # -------------------------------------------------------------# 
    
    
    
    
    
    # --- HC2 adjustment --- #
    h_M = leverages_chol(M) # h_M is the diagonal entries of the projection matrix of M: (M%*%solve(t(M)%*%M)%*%t(M))
    
    est.eps.hat.HC2 = est.eps.hat / sqrt(1 - h_M) 
    
    Hardmard.HC2 = t(est.Xperp)%*% ( ( est.eps.hat.HC2%*%t(est.eps.hat.HC2) ) * Wgt_HAC ) %*% est.Xperp
    
    est.se2.HC2 = Hardmard.HC2 / se2_denom
    
    se2_HAC_HC2[loop.run, loop.L] = est.se2.HC2
    t_HAC_HC2[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2.HC2)
    # ---------------------------- #
    
    
    
    # --- HC3 adjustment --- #
    est.eps.hat.HC3 = est.eps.hat / (1 - h_M) 
    
    # est.eps.hat.HC2 - est.eps.hat.HC3
    
    Hardmard.HC3 = t(est.Xperp)%*% ( ( est.eps.hat.HC3%*%t(est.eps.hat.HC3) ) * Wgt_HAC ) %*% est.Xperp
    
    # max(abs(Hardmard.HC3 - Hardmard.HC2))
    
    est.se2.HC3 = Hardmard.HC3 / se2_denom
    
    # est.se2;est.se2.HC2;est.se2.HC3
    
    se2_HAC_HC3[loop.run, loop.L] = est.se2.HC3
    t_HAC_HC3[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2.HC3)
    # ---------------------------- #
    
    
    print(loop.run)
    
  }
  
  print(loop.L)
  
}

rm(loop.L, loop.run, fit.ols, fit.X, est.Xperp, est.eps.hat, Hardmard, Hardmard.HC2, Hardmard.HC3,
   est.se2, est.se2.HC2, est.se2.HC3, se2_denom)

save(list = c(
  "dat_target", "eps_perp_ori_wide",
  "loading1_wide", "loading2_wide",
  "Q_p_res_c", "Q_p_Xperp_c", "Q_res_c", "Q_Xperp_c",
  "res_c_acf_mat", "res_c_pacf_mat",
  "se2_HAC", "se2_HAC_HC2", "se2_HAC_HC3",
  "t_HAC", "t_HAC_HC2", "t_HAC_HC3",
  "Xperp_c_acf_mat", "Xperp_c_pacf_mat",
  "Z", "Z.eta", "Z.gamm",
  "bet", "eta", "gamm", "Lag.seq", "loading.strength", 
  "loading1", "loading2", "n_state", "n_year",
  "northeast_states", "run", "south_states",
  "state.unq", "statefip_ipt",
  "west_states", "year_ipt"
), 
file = "1108_c_PermuteX&eps_XShock_changeLag.Rdata"
)

load("1108_c_PermuteX&eps_XShock_changeLag.Rdata")

# make plot about non-coverage prob #------------------
# --- HAC --- #
qtl = qnorm(c(0.15, 0.1,0.05, 0.025, 0.005), lower.tail = F) # non-coverage prob of 0.3, 0.2, 0.1, 0.05, 0.01

ncp = matrix(0, nrow = length(Lag.seq), ncol = length(qtl))
rownames(ncp) = Lag.seq
colnames(ncp) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.L in 1:length(Lag.seq)) {
    
    ncp[loop.L, loop.qtl] = sum( abs(t_HAC[,loop.L]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.L, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=Lag.seq, y=ncp[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp[,loop.qtl]),
       xlab = "Lag",
       ylab = "Non-coverage probability HAC",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}





# --- HAC HC2 --- #
ncp_HC2 = matrix(0, nrow = length(Lag.seq), ncol = length(qtl))
rownames(ncp_HC2) = Lag.seq
colnames(ncp_HC2) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.L in 1:length(Lag.seq)) {
    
    ncp_HC2[loop.L, loop.qtl] = sum( abs(t_HAC_HC2[,loop.L]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.L, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=Lag.seq, y=ncp_HC2[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp_HC2[,loop.qtl]),
       xlab = "Lag",
       ylab = "Non-coverage probability HAC HC2",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}

rm(loop.qtl)



# --- HAC HC3 --- #
ncp_HC3 = matrix(0, nrow = length(Lag.seq), ncol = length(qtl))
rownames(ncp_HC3) = Lag.seq
colnames(ncp_HC3) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.L in 1:length(Lag.seq)) {
    
    ncp_HC3[loop.L, loop.qtl] = sum( abs(t_HAC_HC3[,loop.L]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.L, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=Lag.seq, y=ncp_HC3[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp_HC3[,loop.qtl]),
       xlab = "Lag",
       ylab = "Non-coverage probability HAC HC3",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}

rm(loop.qtl)


# --- qqplot --- #
prob = seq(0.01,0.99, by=0.01)

for (loop.L in 1:length(Lag.seq)) {
  qt.emp = quantile(t_HAC[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("Lag = ", Lag.seq[loop.L]))
  abline(a=0, b=1, col=2, lwd=3)
}

for (loop.L in 1:length(Lag.seq)) {
  qt.emp = quantile(t_HAC_HC2[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("Lag = ", Lag.seq[loop.L]))
  abline(a=0, b=1, col=2, lwd=3)
}

for (loop.L in 1:length(Lag.seq)) {
  qt.emp = quantile(t_HAC_HC3[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("Lag = ", Lag.seq[loop.L]))
  abline(a=0, b=1, col=2, lwd=3)
}





# --- Check the size of standard error --- #
plot(x=Lag.seq,
     y=apply(sqrt(se2_HAC), 2, mean),
     type = "b",
     pch=19,
     xlab = "Lag",
     ylab = " ",
     main = "HAC Standard Error")

plot(x=Lag.seq,
     y=apply(sqrt(se2_HAC_HC2), 2, mean),
     type = "b",
     pch=19,
     xlab = "Lag",
     ylab = " ",
     main = "HAC Standard Error HC2")

plot(x=Lag.seq,
     y=apply(sqrt(se2_HAC_HC3), 2, mean),
     type = "b",
     pch=19,
     xlab = "Lag",
     ylab = " ",
     main = "HAC Standard Error HC3")


# It seems that the best Lag is either 20, 40, or 45.














































































#########################################################################
########### Permute Simulation change Shock strength #################### 
########### Permute eps_perp and Xperp together ##########################
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


# ----- start simulation --- #

run = 300

Lag = 50

shock.seq = seq(0,1, by=0.1)

t_HAC = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC = matrix(0, nrow = run, ncol = length(shock.seq))

t_HAC_HC2 = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC_HC2 = matrix(0, nrow = run, ncol = length(shock.seq))

t_HAC_HC3 = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC_HC3 = matrix(0, nrow = run, ncol = length(shock.seq))


loading1 = c(loading1_wide)
loading2 = c(loading2_wide)

n_state = ncol(eps_perp_ori_wide)
n_year = nrow(eps_perp_ori_wide)

state.unq = unique(statefip_ipt)



Wgt.state = matrix(0, n_year, n_year)
# Newey-West with finite sample adjustment
for (Wgtrow in 1:n_year) {
  for (Wgtcol in 1:n_year) {
    Wgt.state[Wgtrow, Wgtcol] = max( 1 - abs(Wgtrow-Wgtcol)/(Lag+1), 0  ) 
  }
}
rm(Wgtrow, Wgtcol)
Wgt_HAC = kronecker(diag(n_state), Wgt.state)

set.seed(715)
for (loop.s in 1:length(shock.seq)) {
  
  loading.strength = shock.seq[loop.s]
  
  
  for (loop.run in 1:run) {
    
    # permute epsperp #---------------------------------------
    
    new_state1 = sample(state.unq, n_state, replace = F)
    
    eps_perm_wide = eps_perp_ori_wide[,as.character(new_state1)]
    
    # dim(eps_perm_wide)
    # eps_perm_wide[,"1"] - eps_perp_ori_wide[,"1"]
    # eps_perm_wide - eps_perp_ori_wide
    
    eps = c(eps_perm_wide)
    
    
    # add aggregate shock to generate eps
    f1 = rnorm(1)
    f2 = rnorm(1)
    
    eps = eps + loading.strength * ( f1*loading1 + f2*loading2 )
    
    # plot(eps)
    
    
    # permute Xperp
    new_state2 = sample(state.unq, n_state, replace = F) # Note that eps and Xperp can not use the same permutation, which will induce endogeneity.
    
    Xperp_perm_wide = Xperp_ori_wide[,as.character(new_state2)]
    
    # dim(Xperp_perm_wide)
    # max(abs(Xperp_perm_wide[,"1"] - Xperp_ori_wide[,"1"]))
    # Xperp_perm_wide - Xperp_ori_wide
    
    Xperp = c(Xperp_perm_wide)
    
    # plot(Xperp)
    
    X  = Z.eta + Xperp 
    
    # plot(X)
    
    # generate y #----------------------------------------
    y = bet * X +  Z.gamm + eps
    
    
    M = cbind(X, Z)
    colnames(M)[1] = "X"
    
    fit.ols=lm(y~M-1)
    
    # HAC using X groups (too slow) #-------------------------------------------------------
    # vcov.hac = solve(t(M)%*%M)%*%(t(M)%*%((fit.ols$residuals%*%t(fit.ols$residuals))*Wgt_HAC)%*%M)%*%solve(t(M)%*%M)
    # t_HAC[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(vcov.hac[1,1])
    # se2_HAC[loop.run, loop.L] = vcov.hac[1,1]
    #---------------------------------------------------------------------------
    
    
    # --- HAC --- #
    fit.X = lm(X~Z-1)
    
    est.Xperp = fit.X$residuals
    
    est.eps.hat = fit.ols$residuals
    
    Hardmard = t(est.Xperp)%*% ( ( est.eps.hat%*%t(est.eps.hat) ) * Wgt_HAC ) %*% est.Xperp
    
    se2_denom = ( sum(est.Xperp*est.Xperp) )^2
    
    est.se2 = Hardmard / se2_denom
    
    # est.se2 - vcov.hac[1,1]
    
    se2_HAC[loop.run, loop.s] = est.se2
    t_HAC[loop.run,loop.s] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2)
    # -------------------------------------------------------------# 
    
    
    
    
    
    # --- HC2 adjustment --- #
    h_M = leverages_chol(M) # h_M is the diagonal entries of the projection matrix of M: (M%*%solve(t(M)%*%M)%*%t(M))
    
    est.eps.hat.HC2 = est.eps.hat / sqrt(1 - h_M) 
    
    Hardmard.HC2 = t(est.Xperp)%*% ( ( est.eps.hat.HC2%*%t(est.eps.hat.HC2) ) * Wgt_HAC ) %*% est.Xperp
    
    est.se2.HC2 = Hardmard.HC2 / se2_denom
    
    se2_HAC_HC2[loop.run, loop.s] = est.se2.HC2
    t_HAC_HC2[loop.run,loop.s] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2.HC2)
    # ---------------------------- #
    
    
    
    # --- HC3 adjustment --- #
    est.eps.hat.HC3 = est.eps.hat / (1 - h_M) 
    
    # est.eps.hat.HC2 - est.eps.hat.HC3
    
    Hardmard.HC3 = t(est.Xperp)%*% ( ( est.eps.hat.HC3%*%t(est.eps.hat.HC3) ) * Wgt_HAC ) %*% est.Xperp
    
    # max(abs(Hardmard.HC3 - Hardmard.HC2))
    
    est.se2.HC3 = Hardmard.HC3 / se2_denom
    
    # est.se2;est.se2.HC2;est.se2.HC3
    
    se2_HAC_HC3[loop.run, loop.s] = est.se2.HC3
    t_HAC_HC3[loop.run,loop.s] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2.HC3)
    # ---------------------------- #
    
    
    print(loop.run)
    
  }
  
  print(loop.s)
  
}

rm(loop.s, loop.run, fit.ols, fit.X, est.Xperp, est.eps.hat, Hardmard, Hardmard.HC2, Hardmard.HC3,
   est.se2, est.se2.HC2, est.se2.HC3, se2_denom)

save(list = c(
  "dat_target", "eps_perp_ori_wide", 
  "loading1_wide", "loading2_wide",
  "Q_p_res_c", "Q_p_Xperp_c",
  "Q_res_c", "Q_Xperp_c",
  "res_c_acf_mat", "res_c_pacf_mat",
  "se2_HAC", "se2_HAC_HC2", "se2_HAC_HC3",
  "t_HAC", "t_HAC_HC2", "t_HAC_HC3",
  "Xperp_c_acf_mat", "Xperp_c_pacf_mat",
  "Xperp_ori_wide",
  "Z", "Z.eta", "Z.gamm", 
  "bet", "eta", "gamm", "Lag",
  "loading1", "loading2",
  "n_state", "n_year",
  "northeast_states", "run", "shock.seq", "south_states",
  "state.unq", "statefip_ipt", "west_states", "year_ipt"
), 
file = "1108_c_PermuteX&eps_epsShock_changeShock_L50.Rdata"
)

load("1108_c_PermuteX&eps_epsShock_changeShock_L20.Rdata")

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





# --- HAC HC2 --- #
ncp_HC2 = matrix(0, nrow = length(shock.seq), ncol = length(qtl))
rownames(ncp_HC2) = shock.seq
colnames(ncp_HC2) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.s in 1:length(shock.seq)) {
    
    ncp_HC2[loop.s, loop.qtl] = sum( abs(t_HAC_HC2[,loop.s]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.s, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=shock.seq, y=ncp_HC2[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp_HC2[,loop.qtl]),
       xlab = "shock strength",
       ylab = "Non-coverage probability HAC HC2",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}

rm(loop.s, loop.qtl)



# --- HAC HC3 --- #
ncp_HC3 = matrix(0, nrow = length(shock.seq), ncol = length(qtl))
rownames(ncp_HC3) = shock.seq
colnames(ncp_HC3) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.s in 1:length(shock.seq)) {
    
    ncp_HC3[loop.s, loop.qtl] = sum( abs(t_HAC_HC3[,loop.s]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.s, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=shock.seq, y=ncp_HC3[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp_HC3[,loop.qtl]),
       xlab = "shock.seq",
       ylab = "Non-coverage probability HAC HC3",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}

rm(loop.qtl)


# --- qqplot --- #
prob = seq(0.01,0.99, by=0.01)

for (loop.s in 1:length(shock.seq)) {
  qt.emp = quantile(t_HAC[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("shock = ", shock.seq[loop.s]))
  abline(a=0, b=1, col=2, lwd=3)
}

for (loop.s in 1:length(shock.seq)) {
  qt.emp = quantile(t_HAC_HC2[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("shock = ", shock.seq[loop.s]))
  abline(a=0, b=1, col=2, lwd=3)
}

for (loop.s in 1:length(shock.seq)) {
  qt.emp = quantile(t_HAC_HC3[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("shock = ", shock.seq[loop.L]))
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

plot(x=shock.seq,
     y=apply(sqrt(se2_HAC_HC2), 2, mean),
     type = "b",
     pch=19,
     xlab = "shock",
     ylab = " ",
     main = "HAC Standard Error HC2")

plot(x=shock.seq,
     y=apply(sqrt(se2_HAC_HC3), 2, mean),
     type = "b",
     pch=19,
     xlab = "Lag",
     ylab = " ",
     main = "HAC Standard Error HC3")













































#########################################################################
########### Permute Simulation change Shock strength #################### 
########### Permut eps_perp and Xperp together ##########################
########### Add aggregate shock by loading1 and loading2 to Xperp ####### 
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


# ----- start simulation --- #

run = 300

Lag = 50  

shock.seq = seq(0,1, by=0.1)

t_HAC = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC = matrix(0, nrow = run, ncol = length(shock.seq))

t_HAC_HC2 = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC_HC2 = matrix(0, nrow = run, ncol = length(shock.seq))

t_HAC_HC3 = matrix(0, nrow = run, ncol = length(shock.seq))
se2_HAC_HC3 = matrix(0, nrow = run, ncol = length(shock.seq))


loading1 = c(loading1_wide)
loading2 = c(loading2_wide)

n_state = ncol(eps_perp_ori_wide)
n_year = nrow(eps_perp_ori_wide)

state.unq = unique(statefip_ipt)



Wgt.state = matrix(0, n_year, n_year)
# Newey-West with finite sample adjustment
for (Wgtrow in 1:n_year) {
  for (Wgtcol in 1:n_year) {
    Wgt.state[Wgtrow, Wgtcol] = max( 1 - abs(Wgtrow-Wgtcol)/(Lag+1), 0  ) 
  }
}
rm(Wgtrow, Wgtcol)
Wgt_HAC = kronecker(diag(n_state), Wgt.state)

set.seed(715)
for (loop.s in 1:length(shock.seq)) {
  
  loading.strength = shock.seq[loop.s]
  #--------------------------------------------------
  
  for (loop.run in 1:run) {
    
    # permute epsperp #---------------------------------------
    
    new_state1 = sample(state.unq, n_state, replace = F)
    
    eps_perm_wide = eps_perp_ori_wide[,as.character(new_state1)]
    
    # dim(eps_perm_wide)
    # eps_perm_wide[,"1"] - eps_perp_ori_wide[,"1"]
    # eps_perm_wide - eps_perp_ori_wide
    
    eps = c(eps_perm_wide)
    
    
    
    
    # plot(eps)
    
    
    # permute Xperp
    new_state2 = sample(state.unq, n_state, replace = F) # Note that eps and Xperp can not use the same permutation, which will induce endogeneity.
    
    Xperp_perm_wide = Xperp_ori_wide[,as.character(new_state2)]
    
    # dim(Xperp_perm_wide)
    # max(abs(Xperp_perm_wide[,"1"] - Xperp_ori_wide[,"1"]))
    # Xperp_perm_wide - Xperp_ori_wide
    
    Xperp = c(Xperp_perm_wide)
    
    # add aggregate shock to generate eps
    f1 = rnorm(1)
    f2 = rnorm(1)
    
    Xperp = Xperp + loading.strength * ( f1*loading1 + f2*loading2 )
    
    # plot(Xperp)
    
    X  = Z.eta + Xperp 
    
    # plot(X)
    
    # generate y #----------------------------------------
    y = bet * X +  Z.gamm + eps
    
    
    M = cbind(X, Z)
    colnames(M)[1] = "X"
    
    fit.ols=lm(y~M-1)
    
    # HAC using X groups (too slow) #-------------------------------------------------------
    # vcov.hac = solve(t(M)%*%M)%*%(t(M)%*%((fit.ols$residuals%*%t(fit.ols$residuals))*Wgt_HAC)%*%M)%*%solve(t(M)%*%M)
    # t_HAC[loop.run,loop.L] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(vcov.hac[1,1])
    # se2_HAC[loop.run, loop.L] = vcov.hac[1,1]
    #---------------------------------------------------------------------------
    
    
    # --- HAC --- #
    fit.X = lm(X~Z-1)
    
    est.Xperp = fit.X$residuals
    
    est.eps.hat = fit.ols$residuals
    
    Hardmard = t(est.Xperp)%*% ( ( est.eps.hat%*%t(est.eps.hat) ) * Wgt_HAC ) %*% est.Xperp
    
    se2_denom = ( sum(est.Xperp*est.Xperp) )^2
    
    est.se2 = Hardmard / se2_denom
    
    # est.se2 - vcov.hac[1,1]
    
    se2_HAC[loop.run, loop.s] = est.se2
    t_HAC[loop.run,loop.s] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2)
    # -------------------------------------------------------------# 
    
    
    
    
    
    # --- HC2 adjustment --- #
    h_M = leverages_chol(M) # h_M is the diagonal entries of the projection matrix of M: (M%*%solve(t(M)%*%M)%*%t(M))
    
    est.eps.hat.HC2 = est.eps.hat / sqrt(1 - h_M) 
    
    Hardmard.HC2 = t(est.Xperp)%*% ( ( est.eps.hat.HC2%*%t(est.eps.hat.HC2) ) * Wgt_HAC ) %*% est.Xperp
    
    est.se2.HC2 = Hardmard.HC2 / se2_denom
    
    se2_HAC_HC2[loop.run, loop.s] = est.se2.HC2
    t_HAC_HC2[loop.run,loop.s] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2.HC2)
    # ---------------------------- #
    
    
    
    # --- HC3 adjustment --- #
    est.eps.hat.HC3 = est.eps.hat / (1 - h_M) 
    
    # est.eps.hat.HC2 - est.eps.hat.HC3
    
    Hardmard.HC3 = t(est.Xperp)%*% ( ( est.eps.hat.HC3%*%t(est.eps.hat.HC3) ) * Wgt_HAC ) %*% est.Xperp
    
    # max(abs(Hardmard.HC3 - Hardmard.HC2))
    
    est.se2.HC3 = Hardmard.HC3 / se2_denom
    
    # est.se2;est.se2.HC2;est.se2.HC3
    
    se2_HAC_HC3[loop.run, loop.s] = est.se2.HC3
    t_HAC_HC3[loop.run,loop.s] = ( summary(fit.ols)$coefficients[1, "Estimate"] - bet ) / sqrt(est.se2.HC3)
    # ---------------------------- #
    
    
    print(loop.run)
    
  }
  
  print(loop.s)
  
}

rm(loop.s, loop.run, fit.ols, fit.X, est.Xperp, est.eps.hat, Hardmard, Hardmard.HC2, Hardmard.HC3,
   est.se2, est.se2.HC2, est.se2.HC3, se2_denom)

save(list = c(
  "dat_target", "eps_perp_ori_wide",
  "loading1_wide", "loading2_wide",
  "Q_p_res_c", "Q_p_Xperp_c",
  "Q_res_c", "Q_Xperp_c",
  "res_c_acf_mat", "res_c_pacf_mat",
  "se2_HAC", "se2_HAC_HC2", "se2_HAC_HC3",
  "t_HAC", "t_HAC_HC2", "t_HAC_HC3",
  "Xperp_c_acf_mat", "Xperp_c_pacf_mat",
  "Xperp_ori_wide", 
  "Z", "Z.eta", "Z.gamm",
  "bet", "eta", "gamm", "Lag",
  "loading1", "loading2",
  "n_state", "n_year", "northeast_states",
  "run", "shock.seq", "south_states", "state.unq",
  "statefip_ipt", "west_states", "year_ipt"
), 
file = "1108_c_PermuteX&eps_XShock_changeShock_L50.Rdata"
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





# --- HAC HC2 --- #
ncp_HC2 = matrix(0, nrow = length(shock.seq), ncol = length(qtl))
rownames(ncp_HC2) = shock.seq
colnames(ncp_HC2) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.s in 1:length(shock.seq)) {
    
    ncp_HC2[loop.s, loop.qtl] = sum( abs(t_HAC_HC2[,loop.s]) > qtl[loop.qtl] ) / run
    
  }
  
}

rm(loop.s, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=shock.seq, y=ncp_HC2[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp_HC2[,loop.qtl]),
       xlab = "shock strength",
       ylab = "Non-coverage probability HAC HC2",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}

rm(loop.s, loop.qtl)



# --- HAC HC3 --- #
ncp_HC3 = matrix(0, nrow = length(shock.seq), ncol = length(qtl))
rownames(ncp_HC3) = shock.seq
colnames(ncp_HC3) = c("nominal=0.3", "nominal=0.2", "nominal=0.1", "nominal=0.05", "nominal=0.01")

for (loop.qtl in 1:length(qtl)){
  
  for (loop.s in 1:length(shock.seq)) {
    
    ncp_HC3[loop.s, loop.qtl] = sum( abs(t_HAC_HC3[,loop.s]) > qtl[loop.qtl] ) / run
    
  }
  
}


rm(loop.s, loop.qtl)

for (loop.qtl in 1:length(qtl)) {
  
  plot(x=shock.seq, y=ncp_HC3[,loop.qtl], 
       type = "b", lwd=2, pch=19,
       col=2,
       ylim = range(ncp_HC3[,loop.qtl]),
       xlab = "shock.seq",
       ylab = "Non-coverage probability HAC HC3",
       main = paste0("nominal level =  ", round( (2*pnorm(qtl[loop.qtl], lower.tail = F)), digits = 2 ) )
  )
  abline(h=2*pnorm(qtl[loop.qtl], lower.tail = F), lty=2, col=1)
  
}

rm(loop.qtl)


# --- qqplot --- #
prob = seq(0.01,0.99, by=0.01)

for (loop.s in 1:length(shock.seq)) {
  qt.emp = quantile(t_HAC[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("shock = ", shock.seq[loop.s]))
  abline(a=0, b=1, col=2, lwd=3)
}

for (loop.s in 1:length(shock.seq)) {
  qt.emp = quantile(t_HAC_HC2[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("shock = ", shock.seq[loop.s]))
  abline(a=0, b=1, col=2, lwd=3)
}

for (loop.s in 1:length(shock.seq)) {
  qt.emp = quantile(t_HAC_HC3[,loop.L], probs = prob)
  plot(x = qt.emp, y = qnorm(prob), main = paste0("shock = ", shock.seq[loop.L]))
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

plot(x=shock.seq,
     y=apply(sqrt(se2_HAC_HC2), 2, mean),
     type = "b",
     pch=19,
     xlab = "shock",
     ylab = " ",
     main = "HAC Standard Error HC2")

plot(x=shock.seq,
     y=apply(sqrt(se2_HAC_HC3), 2, mean),
     type = "b",
     pch=19,
     xlab = "shock",
     ylab = " ",
     main = "HAC Standard Error HC3")




boxplot(t_HAC)
boxplot(t_HAC_HC2)





















#######################################################
## There is a big problem: Xperp and eps from the real data have very different variances
## When I add the shock, the shock strengths are the same to be added on both Xperp or eps, which make them not comparable.
## Use the standard deviation of Xperp or eps from the real data, times the sqrt(1-loading strength) * permute + sqrt(loading strength) * shock
## This should be fair
#############################################################

###########################################################
# But I still do not know how why the SE is decreasing for the shock on Xperp as the shock strength increases.
############################################################
