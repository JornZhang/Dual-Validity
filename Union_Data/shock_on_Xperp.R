#########################################################################
## HAC simulation for one shock strength, designed for SLURM array jobs
##
## Usage on cluster:
##   Rscript shock_on_Xperp.R <task_id> <run> <Lag>
##
## Example:
##   Rscript shock_on_Xperp.R 3 10000 20
##
## task_id selects shock.seq[task_id], i.e., kappa. By default shock.seq has length 10.
## The script assumes that Data_Setting.Rdata is in the working
## directory unless INPUT_RDATA is set as an environment variable.
#########################################################################

rm(list = ls())
options(warn = 1)

## ----------------------------- arguments ----------------------------- ##
args <- commandArgs(trailingOnly = TRUE)

get_int_arg <- function(i, default) {
  if (length(args) >= i && nzchar(args[i])) {
    return(as.integer(args[i]))
  }
  as.integer(default)
}

## First argument: shock index. Prefer command-line arg; otherwise SLURM array id.
if (length(args) >= 1 && nzchar(args[1])) {
  loop.s <- as.integer(args[1])
} else {
  loop.s <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
}

run <- get_int_arg(2, Sys.getenv("RUN", unset = "10000"))
Lag <- get_int_arg(3, Sys.getenv("LAG", unset = "20"))

input_rdata <- Sys.getenv("INPUT_RDATA", unset = "1108_PermSimu_Setting_c.Rdata")
output_dir  <- Sys.getenv("OUTPUT_DIR", unset = "results_0707")
base_seed   <- as.integer(Sys.getenv("BASE_SEED", unset = "1130"))

## Set this to TRUE only if you want the generated AR(1) processes to have
## marginal standard deviations equal to sd(eps_perp_ori_wide) and sd(Xperp_ori_wide).
## The default FALSE reproduces your 0621.R convention, where those empirical SDs
## are used as innovation SDs.
use_marginal_sd <- identical(tolower(Sys.getenv("USE_MARGINAL_SD", unset = "false")), "true")

shock.seq <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

if (is.na(loop.s) || loop.s < 1 || loop.s > length(shock.seq)) {
  stop("loop.s/task_id must be an integer from 1 to ", length(shock.seq), ".")
}
if (is.na(run) || run < 1) stop("run must be a positive integer.")
if (is.na(Lag) || Lag < 1) stop("Lag must be a positive integer.")
if (!file.exists(input_rdata)) stop("Cannot find input RData file: ", input_rdata)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

loading.strength <- shock.seq[loop.s]
nominal <- 0.05

cat("============================================================\n")
cat("Shock_on_Xperp\n")
cat("task_id / loop.s      :", loop.s, "\n")
cat("loading.strength      :", loading.strength, "\n")
cat("run                   :", run, "\n")
cat("Lag                   :", Lag, "\n")
cat("input_rdata           :", input_rdata, "\n")
cat("output_dir            :", output_dir, "\n")
cat("base_seed             :", base_seed, "\n")
cat("use_marginal_sd       :", use_marginal_sd, "\n")
cat("============================================================\n")

## ------------------------------- load data --------------------------- ##
load(input_rdata)

required_objects <- c(
  "Z", "eta", "gamm", "bet",
  "eps_perp_ori_wide", "Xperp_ori_wide",
  "loading1_wide", "loading2_wide",
  "res_c_acf_mat", "Xperp_c_acf_mat"
)
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop("The input RData file is missing: ", paste(missing_objects, collapse = ", "))
}

Z <- as.matrix(Z)
eta <- as.numeric(eta)
gamm <- as.numeric(gamm)
bet <- as.numeric(bet)[1]

n_state <- ncol(eps_perp_ori_wide)
n_year  <- nrow(eps_perp_ori_wide)
N <- n_state * n_year

Z.eta  <- as.numeric(Z %*% eta)
Z.gamm <- as.numeric(Z %*% gamm)

if (length(Z.eta) != N || length(Z.gamm) != N) {
  stop("Length mismatch: length(Z %*% eta) or length(Z %*% gamm) is not n_year * n_state.")
}

loading1 <- as.numeric(c(loading1_wide))
loading2 <- as.numeric(c(loading2_wide))
if (length(loading1) != N || length(loading2) != N) {
  stop("Length mismatch: loading1/loading2 must have length n_year * n_state.")
}

sd.eps.raw    <- sd(as.numeric(c(eps_perp_ori_wide)), na.rm = TRUE)
sd.Xperp.raw <- sd(as.numeric(c(Xperp_ori_wide)), na.rm = TRUE)

clip_phi <- function(x, lower = -0.95, upper = 0.95) {
  x <- as.numeric(x)[1]
  if (!is.finite(x)) stop("Non-finite AR(1) coefficient.")
  max(min(x, upper), lower)
}

phi_eps    <- clip_phi(mean(res_c_acf_mat$ACF_mat[2, ], na.rm = TRUE))
phi_Xperp <- clip_phi(mean(Xperp_c_acf_mat$ACF_mat[2, ], na.rm = TRUE))

if (use_marginal_sd) {
  sd.eps.innov    <- sd.eps.raw    * sqrt(1 - phi_eps^2)
  sd.Xperp.innov <- sd.Xperp.raw * sqrt(1 - phi_Xperp^2)
} else {
  sd.eps.innov    <- sd.eps.raw
  sd.Xperp.innov <- sd.Xperp.raw
}

## Rank used for HC1 correction. This is fixed because Z is fixed and X is generated
## outside the column space of Z in the simulation.
k_full <- ncol(Z) + 1
if (k_full >= N) stop("k_full >= N, so HC1 correction is not well-defined.")
hc1_t_factor <- sqrt((N - k_full) / N)

cat("n_year               :", n_year, "\n")
cat("n_state              :", n_state, "\n")
cat("N                    :", N, "\n")
cat("k_full               :", k_full, "\n")
cat("phi_eps              :", phi_eps, "\n")
cat("phi_Xperp            :", phi_Xperp, "\n")
cat("sd.eps.raw           :", sd.eps.raw, "\n")
cat("sd.Xperp.raw         :", sd.Xperp.raw, "\n")
cat("sd.eps.innov         :", sd.eps.innov, "\n")
cat("sd.Xperp.innov       :", sd.Xperp.innov, "\n")
cat("HC1 t multiplier     :", hc1_t_factor, "\n")
cat("============================================================\n")

## ------------------------------ functions ---------------------------- ##
simulate_ar1 <- function(n_year, n_state, phi, sigma) {
  innovations <- matrix(rnorm(n_year * n_state, sd = sigma),
                        nrow = n_year, ncol = n_state)
  x <- matrix(NA_real_, nrow = n_year, ncol = n_state)
  x[1, ] <- innovations[1, ] / sqrt(1 - phi^2)
  for (tt in 2:n_year) {
    x[tt, ] <- phi * x[tt - 1, ] + innovations[tt, ]
  }
  x
}

make_bartlett_weight <- function(Tn, Lag) {
  idx <- seq_len(Tn)
  outer(idx, idx, function(i, j) pmax(1 - abs(i - j) / (Lag + 1), 0))
}

hac_quad_bartlett <- function(g, Lag) {
  g <- as.numeric(g)
  Tn <- length(g)
  if (Tn < 1) return(0)
  Lag <- min(as.integer(Lag), Tn - 1)

  out <- sum(g^2)
  if (Lag >= 1) {
    for (ell in seq_len(Lag)) {
      w <- 1 - ell / (Lag + 1)
      out <- out + 2 * w * sum(g[(ell + 1):Tn] * g[1:(Tn - ell)])
    }
  }
  as.numeric(out)
}

hac_quad_prewhite_ar1 <- function(g, Lag) {
  g <- as.numeric(g)
  Tn <- length(g)
  if (Tn <= 2) return(hac_quad_bartlett(g, Lag = min(Lag, Tn - 1)))

  g_lag <- g[-Tn]
  g_now <- g[-1]

  denom <- sum(g_lag^2)
  if (!is.finite(denom) || denom <= 1e-12) {
    rho_hat <- 0
  } else {
    rho_hat <- sum(g_now * g_lag) / denom
  }
  if (!is.finite(rho_hat)) rho_hat <- 0
  rho_hat <- max(min(rho_hat, 0.95), -0.95)

  v <- g_now - rho_hat * g_lag

  ## v has length Tn - 1. The factor below is a small length correction.
  as.numeric((Tn / length(v)) *
               hac_quad_bartlett(v, Lag = min(Lag, length(v) - 1)) /
               (1 - rho_hat)^2)
}

satterthwaite_df <- function(q_state, max_df) {
  q_state <- as.numeric(q_state)
  q_state <- q_state[is.finite(q_state)]
  if (length(q_state) == 0) return(2)

  ## The state-level variance contributions should be nonnegative. This guards
  ## against tiny numerical negative values.
  q_pos <- pmax(q_state, 0)
  s1 <- sum(q_pos)
  s2 <- sum(q_pos^2)
  if (!is.finite(s1) || !is.finite(s2) || s1 <= 0 || s2 <= 0) return(2)

  df <- 2 * s1^2 / s2
  df <- max(min(df, max_df), 2)
  as.numeric(df)
}

safe_t <- function(beta_hat, beta_true, se2) {
  if (!is.finite(se2) || se2 <= 0) return(NA_real_)
  as.numeric((beta_hat - beta_true) / sqrt(se2))
}

calc_coverage <- function(t_stat, df = NULL) {
  ok <- is.finite(t_stat)
  if (!is.null(df)) ok <- ok & is.finite(df)
  if (!any(ok)) return(NA_real_)

  if (is.null(df)) {
    crit <- qnorm(1 - nominal / 2)
  } else {
    crit <- qt(1 - nominal / 2, df = df[ok])
  }
  mean(abs(t_stat[ok]) <= crit)
}

## -------------------------- pre-computed weight ---------------------- ##
Wgt.state <- make_bartlett_weight(n_year, Lag)

## ---------------------------- simulation ----------------------------- ##
## Same seed for every shock level. This makes the Monte Carlo draws matched
## across SLURM array tasks, differing only through loading.strength.
set.seed(base_seed)

## Usual HAC vectors
se2_HAC <- numeric(run)
t_HAC_raw <- numeric(run)
t_HAC_HC1 <- numeric(run)
df_HAC <- numeric(run)

## Prewhitened HAC vectors
se2_HAC_PW <- numeric(run)
t_HAC_PW_raw <- numeric(run)
t_HAC_PW_HC1 <- numeric(run)
df_HAC_PW <- numeric(run)

for (loop.run in seq_len(run)) {
  ## ----- permute eps_perp via AR(1) ----- ##
  eps_ar1_wide <- simulate_ar1(
    n_year  = n_year,
    n_state = n_state,
    sigma   = sd.eps.innov,
    phi     = phi_eps
  )
  eps <- c(eps_ar1_wide)
  

  ## ----- permute Xperp via AR(1) ----- ##
  Xperp_ar1_wide <- simulate_ar1(
    n_year  = n_year,
    n_state = n_state,
    sigma   = sd.Xperp.innov,
    phi     = phi_Xperp
  )
  Xperp_ar1 <- c(Xperp_ar1_wide)
  
  ## add aggregate shock
  f1 <- rnorm(1)
  f2 <- rnorm(1)
  
  Xperp <- sqrt(1 - loading.strength) * Xperp_ar1 +
    sqrt(loading.strength) * sd.Xperp.raw * (f1 * loading1 + f2 * loading2)
  
  

  ## Construct X and y.
  X <- Z.eta + Xperp
  y <- bet * X + Z.gamm + eps

  ## OLS fit using lm.fit for speed.
  M <- cbind(X = X, Z)
  fit.ols <- lm.fit(x = M, y = y)
  beta_hat <- as.numeric(fit.ols$coefficients[1])
  est.eps.hat <- as.numeric(fit.ols$residuals)

  ## Residualize X on Z.
  fit.X <- lm.fit(x = Z, y = X)
  est.Xperp <- as.numeric(fit.X$residuals)

  ## HAC score representation, avoiding the huge N x N Hadamard matrix.
  score <- est.Xperp * est.eps.hat
  score_wide <- matrix(score, nrow = n_year, ncol = n_state)
  se2_denom <- sum(est.Xperp^2)^2

  ## Usual Bartlett HAC.
  q_state <- colSums(score_wide * (Wgt.state %*% score_wide))
  hardmard <- sum(q_state)
  se2 <- as.numeric(hardmard / se2_denom)

  se2_HAC[loop.run] <- se2
  t_HAC_raw[loop.run] <- safe_t(beta_hat, bet, se2)
  t_HAC_HC1[loop.run] <- t_HAC_raw[loop.run] * hc1_t_factor
  df_HAC[loop.run] <- satterthwaite_df(q_state, max_df = n_state - 1)

  ## Prewhitened HAC on state-level score series.
  q_state_PW <- apply(score_wide, 2, hac_quad_prewhite_ar1, Lag = Lag)
  hardmard_PW <- sum(q_state_PW)
  se2_PW <- as.numeric(hardmard_PW / se2_denom)

  se2_HAC_PW[loop.run] <- se2_PW
  t_HAC_PW_raw[loop.run] <- safe_t(beta_hat, bet, se2_PW)
  t_HAC_PW_HC1[loop.run] <- t_HAC_PW_raw[loop.run] * hc1_t_factor
  df_HAC_PW[loop.run] <- satterthwaite_df(q_state_PW, max_df = n_state - 1)

  if (loop.run %% 100 == 0 || loop.run == run) {
    cat("task", loop.s,
        "kappa", loading.strength,
        "rep", loop.run, "of", run,
        "mean_df", round(mean(df_HAC[seq_len(loop.run)], na.rm = TRUE), 2),
        "mean_df_PW", round(mean(df_HAC_PW[seq_len(loop.run)], na.rm = TRUE), 2),
        "\n")
  }
}

## ------------------------------ summary ------------------------------ ##
summary_df <- data.frame(
  task_id = loop.s,
  loading_strength = loading.strength,
  Lag = Lag,
  run = run,
  n_year = n_year,
  n_state = n_state,
  N = N,
  k_full = k_full,
  phi_eps = phi_eps,
  phi_Xperp = phi_Xperp,
  sd_eps_raw = sd.eps.raw,
  sd_Xperp_raw = sd.Xperp.raw,
  sd_eps_innov = sd.eps.innov,
  sd_Xperp_innov = sd.Xperp.innov,
  method = c(
    "HAC_raw_norm",
    "HAC_HC1_norm",
    "HAC_HC1_Satterthwaite",
    "HAC_PW_raw_norm",
    "HAC_PW_HC1_norm",
    "HAC_PW_HC1_Satterthwaite"
  ),
  coverage = c(
    calc_coverage(t_HAC_raw),
    calc_coverage(t_HAC_HC1),
    calc_coverage(t_HAC_HC1, df = df_HAC),
    calc_coverage(t_HAC_PW_raw),
    calc_coverage(t_HAC_PW_HC1),
    calc_coverage(t_HAC_PW_HC1, df = df_HAC_PW)
  ),
  mean_df = c(NA_real_, NA_real_, mean(df_HAC, na.rm = TRUE),
              NA_real_, NA_real_, mean(df_HAC_PW, na.rm = TRUE)),
  median_df = c(NA_real_, NA_real_, median(df_HAC, na.rm = TRUE),
                NA_real_, NA_real_, median(df_HAC_PW, na.rm = TRUE)),
  mean_se2 = c(mean(se2_HAC, na.rm = TRUE), mean(se2_HAC, na.rm = TRUE), mean(se2_HAC, na.rm = TRUE),
               mean(se2_HAC_PW, na.rm = TRUE), mean(se2_HAC_PW, na.rm = TRUE), mean(se2_HAC_PW, na.rm = TRUE))
)
summary_df$noncoverage <- 1 - summary_df$coverage

settings <- list(
  task_id = loop.s,
  loading_strength = loading.strength,
  shock.seq = shock.seq,
  run = run,
  Lag = Lag,
  nominal = nominal,
  input_rdata = input_rdata,
  output_dir = output_dir,
  base_seed = base_seed,
  same_seed_across_shocks = TRUE,
  use_marginal_sd = use_marginal_sd,
  n_year = n_year,
  n_state = n_state,
  N = N,
  k_full = k_full,
  hc1_t_factor = hc1_t_factor,
  phi_eps = phi_eps,
  phi_Xperp = phi_Xperp,
  sd.eps.raw = sd.eps.raw,
  sd.Xperp.raw = sd.Xperp.raw,
  sd.eps.innov = sd.eps.innov,
  sd.Xperp.innov = sd.Xperp.innov
)

reslt <- list(
  settings = settings,
  summary = summary_df,
  t = t_HAC_raw,
  t_HC1 = t_HAC_HC1,
  se2 = se2_HAC,
  df = df_HAC,
  t_PW = t_HAC_PW_raw,
  t_PW_HC1 = t_HAC_PW_HC1,
  se2_PW = se2_HAC_PW,
  df_PW = df_HAC_PW
)

kappa_label <- formatC(loading.strength, format = "f", digits = 1)
rds_file <- file.path(output_dir,
                      sprintf("Shock_on_Xperp_L%d_task%02d_kappa%s.rds", Lag, loop.s, kappa_label))
csv_file <- file.path(output_dir,
                      sprintf("Shock_on_Xperp_L%d_task%02d_kappa%s_summary.csv", Lag, loop.s, kappa_label))

saveRDS(reslt, file = rds_file)
write.csv(summary_df, file = csv_file, row.names = FALSE)

cat("============================================================\n")
cat("Finished task", loop.s, "with kappa", loading.strength, "\n")
cat("Saved RDS    :", rds_file, "\n")
cat("Saved summary:", csv_file, "\n")
print(summary_df)
cat("============================================================\n")
