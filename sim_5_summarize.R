# =============================================================================
# sim_5_summarize.R
# Compute Bias, SD, ESE, CP95 and generate LaTeX tables
# for ECM, Gibbs, and Fang methods
#
# Table structure: theta1 in {-0.5, 0, 0.5} x theta2 in {0, 0.5}
#                  alpha/beta in {1, 0.95, 0.90, 0.85}
#                  500 replicates per setting
# =============================================================================

library(xtable)

# =============================================================================
# 1. Load results
# =============================================================================
load_results <- function(path, n_batch = 24) {
  theta1 <- do.call(rbind, lapply(1:n_batch, function(i)
    read.table(paste0(path, "theta1.", i, ".txt"), header = TRUE)))
  theta2 <- do.call(rbind, lapply(1:n_batch, function(i)
    read.table(paste0(path, "theta2.", i, ".txt"), header = TRUE)))
  list(theta1 = as.matrix(theta1), theta2 = as.matrix(theta2))
}

ECM   <- load_results("simulation/results/ECM/")
Gibbs <- load_results("simulation/results/Gibbs/")
Fang  <- load_results("simulation/results/Fang/")

# =============================================================================
# 2. Parameter settings (matches sim_data_generation.R)
# =============================================================================
# 12000 datasets: 24 settings x 500 replicates
# Setting order:
#   theta1 in {-0.5, 0, 0.5} x theta2 in {0, 0.5} = 6 combinations
#   x alpha/beta in {1, 0.95, 0.9, 0.85} = 24 settings
#   each repeated 500 times
#
# IT index:
#   alpha=1.00: IT 1    ~ 3000  (settings 1~6,   each 500 reps)
#   alpha=0.95: IT 3001 ~ 6000  (settings 7~12,  each 500 reps)
#   alpha=0.90: IT 6001 ~ 9000  (settings 13~18, each 500 reps)
#   alpha=0.85: IT 9001 ~ 12000 (settings 19~24, each 500 reps)
#
# Within each alpha block, 6 settings x 500 reps:
#   theta1=-0.5, theta2=0   : 1   ~ 500
#   theta1=-0.5, theta2=0.5 : 501 ~ 1000
#   theta1=0,    theta2=0   : 1001~ 1500
#   theta1=0,    theta2=0.5 : 1501~ 2000
#   theta1=0.5,  theta2=0   : 2001~ 2500
#   theta1=0.5,  theta2=0.5 : 2501~ 3000

theta1_vals <- c(-0.5, -0.5, 0, 0, 0.5, 0.5)
theta2_vals <- c(0, 0.5, 0, 0.5, 0, 0.5)
alpha_vals  <- c(1, 0.95, 0.9, 0.85)

# =============================================================================
# 3. Compute metrics for one method
# =============================================================================
# theta_mat: 12000 x 4 matrix (EST, SE, lb95, ub95) for one parameter
# true_val: true parameter value
# setting_idx: which 500 rows to use (1~500, 501~1000, etc.)

compute_metrics <- function(theta_mat, true_val, rows) {
  est  <- theta_mat[rows, 1]
  se   <- theta_mat[rows, 2]
  lb   <- theta_mat[rows, 3]
  ub   <- theta_mat[rows, 4]

  # Remove failed cases (est == 0 and se == 0, or NA in se)
  ok   <- !(est == 0 & se == 0) & !is.na(se) & !is.na(est)
  n_ok <- sum(ok)
  est  <- est[ok]; se <- se[ok]; lb <- lb[ok]; ub <- ub[ok]

  bias <- mean(est) - true_val
  sd_  <- sd(est)
  ese  <- mean(se, na.rm = TRUE)
  cp95 <- mean(lb <= true_val & true_val <= ub, na.rm = TRUE)

  c(Bias = round(bias, 3),
    SD   = round(sd_,  3),
    ESE  = round(ese,  3),
    CP95 = round(cp95, 3))
}

# =============================================================================
# 4. Build summary table for one method and one alpha/beta pair
# =============================================================================
build_table <- function(res, alpha_idx) {
  # alpha_idx: 1=1.00, 2=0.95, 3=0.90, 4=0.85
  base <- (alpha_idx - 1) * 3000

  rows_list <- lapply(1:6, function(s) base + (500*(s-1)+1):(500*s))

  out <- data.frame()
  for (s in 1:6) {
    t1_true <- theta1_vals[s]
    t2_true <- theta2_vals[s]
    rows    <- rows_list[[s]]

    m1 <- compute_metrics(res$theta1, t1_true, rows)
    m2 <- compute_metrics(res$theta2, t2_true, rows)

    out <- rbind(out,
      data.frame(theta1=t1_true, theta2=t2_true, param="theta1",
                 Bias=m1["Bias"], SD=m1["SD"], ESE=m1["ESE"], CP95=m1["CP95"]),
      data.frame(theta1=t1_true, theta2=t2_true, param="theta2",
                 Bias=m2["Bias"], SD=m2["SD"], ESE=m2["ESE"], CP95=m2["CP95"])
    )
  }
  rownames(out) <- NULL
  out
}

# =============================================================================
# 5. Generate LaTeX table — new format (3 methods x theta1/theta2 panels)
# =============================================================================
make_latex_table <- function(alpha_idx, alpha_label) {

  t_ecm   <- build_table(ECM,   alpha_idx)
  t_gibbs <- build_table(Gibbs, alpha_idx)
  t_fang  <- build_table(Fang,  alpha_idx)

  cat(sprintf("\n%% ============================================================\n"))
  cat(sprintf("%% alpha/beta = %s\n", alpha_label))
  cat(sprintf("%% ============================================================\n"))

  cat("\\begin{table}[]\n")
  cat("\\begin{adjustbox}{center}\n")
  cat("\\begin{tabular}{ccclrrrrcrrrr}\\toprule\\toprule\n")
  cat(sprintf(" &  &  &  & \\multicolumn{4}{c}{$\\hat{\\theta}_1$} & & \\multicolumn{4}{c}{$\\hat{\\theta}_2$} \\\\\n"))
  cat("\\cline{5-8} \\cline{10-13}\n")
  cat("$\\theta_1$ & $\\theta_2$ & & Method & Bias & SD & ESE & CP95 & & Bias & SD & ESE & CP95 \\\\ \\hline\\rule{0pt}{12pt}\n")

  prev_t1 <- NA
  for (s in 1:6) {
    t1_val <- c(-0.5, -0.5, 0, 0, 0.5, 0.5)[s]
    t2_val <- c(0, 0.5, 0, 0.5, 0, 0.5)[s]

    r_ecm   <- t_ecm[t_ecm$theta1 == t1_val & t_ecm$theta2 == t2_val, ]
    r_gibbs <- t_gibbs[t_gibbs$theta1 == t1_val & t_gibbs$theta2 == t2_val, ]
    r_fang  <- t_fang[t_fang$theta1 == t1_val & t_fang$theta2 == t2_val, ]

    ecm1   <- r_ecm[r_ecm$param == "theta1", ]
    ecm2   <- r_ecm[r_ecm$param == "theta2", ]
    gibbs1 <- r_gibbs[r_gibbs$param == "theta1", ]
    gibbs2 <- r_gibbs[r_gibbs$param == "theta2", ]
    fang1  <- r_fang[r_fang$param == "theta1", ]
    fang2  <- r_fang[r_fang$param == "theta2", ]

    t1_str <- if (!is.na(prev_t1) && t1_val == prev_t1) "" else as.character(t1_val)
    prev_t1 <- t1_val

    cat(sprintf("%s & %s & & ECM   & %.3f & %.3f & %.3f & %.3f & & %.3f & %.3f & %.3f & %.3f \\\\\n",
                t1_str, t2_val,
                ecm1$Bias, ecm1$SD, ecm1$ESE, ecm1$CP95,
                ecm2$Bias, ecm2$SD, ecm2$ESE, ecm2$CP95))
    cat(sprintf(" &  & & Gibbs & %.3f & %.3f & %.3f & %.3f & & %.3f & %.3f & %.3f & %.3f \\\\\n",
                gibbs1$Bias, gibbs1$SD, gibbs1$ESE, gibbs1$CP95,
                gibbs2$Bias, gibbs2$SD, gibbs2$ESE, gibbs2$CP95))
    cat(sprintf(" &  & & Li et al. & %.3f & %.3f & %.3f & %.3f & & %.3f & %.3f & %.3f & %.3f \\\\\n",
                fang1$Bias, fang1$SD, fang1$ESE, fang1$CP95,
                fang2$Bias, fang2$SD, fang2$ESE, fang2$CP95))
  }

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\end{adjustbox}\n")
  cat(sprintf("\\caption{Simulation results with $g(t)=\\log(t/20)+t$, based on 500 data sets with $(\\alpha, \\beta)=(%s, %s)$; Bias is the average bias, SD is the sample standard deviation, ESE is the average of the estimated standard errors, and CP95 the 95 percent coverage probability.}\n",
              alpha_label, alpha_label))
  cat("\\end{table}\n\n")
}

# =============================================================================
# 6. Generate all tables
# =============================================================================
if (sink.number() > 0) sink(type = "output")
sink("simulation/results/tables.tex")

make_latex_table(1, "1")
make_latex_table(2, "0.95")
make_latex_table(3, "0.90")
make_latex_table(4, "0.85")

sink()
if (sink.number() > 0) sink(type = "output")

cat("Done! Tables saved to simulation/results/tables.tex\n")

# =============================================================================
# 7. NA rate summary
# =============================================================================
alpha_labels <- c("1.00", "0.95", "0.90", "0.85")

na_summary <- function(path, method) {
  cat(sprintf("\n=== %s: NA rate by alpha ===\n", method))
  for (i in 1:4) {
    batches <- (3*(i-1)+1):(3*i)
    theta1  <- do.call(rbind, lapply(batches, function(b)
      read.table(paste0(path, "theta1.", b, ".txt"), header = TRUE)))
    n_total <- nrow(theta1)
    failed  <- (theta1[,1] == 0 & theta1[,2] == 0) | is.na(theta1[,2])
    n_fail  <- sum(failed)
    cat(sprintf("  alpha=beta=%s: %d / %d failed (%.1f%%)\n",
                alpha_labels[i], n_fail, n_total, 100*n_fail/n_total))
  }
}

na_summary("simulation/results/ECM/",   "ECM")
na_summary("simulation/results/Gibbs/", "Gibbs")
na_summary("simulation/results/Fang/",  "Fang")

# =============================================================================
# 8. Computing time summary + LaTeX table
# =============================================================================
get_times <- function(path) {
  times <- unlist(lapply(1:24, function(i)
    read.table(paste0(path, "comp.time.", i, ".txt"), header = TRUE)[, 1]))
  times[times > 0]
}

print_time_summary <- function(times, method) {
  cat(sprintf("\n=== %s: Computing time (seconds per dataset) ===\n", method))
  cat(sprintf("  Mean : %.2f sec\n", mean(times)))
  cat(sprintf("  SD   : %.2f sec\n", sd(times)))
  cat(sprintf("  Min  : %.2f sec\n", min(times)))
  cat(sprintf("  Max  : %.2f sec\n", max(times)))
  cat(sprintf("  Total: %.1f hours\n", sum(times) / 3600))
}

ecm_t   <- get_times("simulation/results/ECM/")
gibbs_t <- get_times("simulation/results/Gibbs/")
fang_t  <- get_times("simulation/results/Fang/")

print_time_summary(ecm_t,   "ECM")
print_time_summary(gibbs_t, "Gibbs")
print_time_summary(fang_t,  "Fang")

cat("\n=== LaTeX time comparison table ===\n")
cat("\\begin{table}[]\n")
cat("\\centering\n")
cat("\\begin{tabular}{lcc}\\toprule\n")
cat("Method & Mean time (sec) & SD \\\\ \\midrule\n")
cat(sprintf("ECM          & %.2f & %.2f \\\\\n", mean(ecm_t),   sd(ecm_t)))
cat(sprintf("Gibbs        & %.2f & %.2f \\\\\n", mean(gibbs_t), sd(gibbs_t)))
cat(sprintf("Fang et al.\\ & %.2f & %.2f \\\\\n", mean(fang_t),  sd(fang_t)))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\caption{Average computing time per dataset (seconds).}\n")
cat("\\end{table}\n")

# =============================================================================
# 9. Append NA rate and time tables to tables.tex
# =============================================================================
sink("simulation/results/tables.tex", append = TRUE)

# NA rate table
cat("\n% ============================================================\n")
cat("% Convergence failure rates\n")
cat("% ============================================================\n")
cat("\\begin{table}[]\n")
cat("\\centering\n")
cat("\\begin{tabular}{lcccc}\\toprule\n")
cat("Method & $\\alpha,\\beta=1$ & $\\alpha,\\beta=0.95$ & $\\alpha,\\beta=0.90$ & $\\alpha,\\beta=0.85$ \\\\ \\midrule\n")

get_na_rates <- function(path) {
  sapply(1:4, function(i) {
    batches <- (3*(i-1)+1):(3*i)
    theta1  <- do.call(rbind, lapply(batches, function(b)
      read.table(paste0(path, "theta1.", b, ".txt"), header = TRUE)))
    failed  <- (theta1[,1] == 0 & theta1[,2] == 0) | is.na(theta1[,2])
    round(100 * mean(failed), 1)
  })
}

ecm_na   <- get_na_rates("simulation/results/ECM/")
gibbs_na <- get_na_rates("simulation/results/Gibbs/")
fang_na  <- get_na_rates("simulation/results/Fang/")

cat(sprintf("ECM          & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\\n",
            ecm_na[1], ecm_na[2], ecm_na[3], ecm_na[4]))
cat(sprintf("Gibbs        & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\\n",
            gibbs_na[1], gibbs_na[2], gibbs_na[3], gibbs_na[4]))
cat(sprintf("Fang et al.\\ & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\\n",
            fang_na[1], fang_na[2], fang_na[3], fang_na[4]))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\caption{Convergence failure rates (\\%) by method and sensitivity/specificity setting.}\n")
cat("\\end{table}\n\n")

# Time table
cat("% ============================================================\n")
cat("% Computing time comparison\n")
cat("% ============================================================\n")
cat("\\begin{table}[]\n")
cat("\\centering\n")
cat("\\begin{tabular}{lrr}\\toprule\n")
cat("Method & Mean (sec) & SD (sec) \\\\ \\midrule\n")
cat(sprintf("ECM          & %.2f & %.2f \\\\\n", mean(ecm_t),   sd(ecm_t)))
cat(sprintf("Gibbs        & %.2f & %.2f \\\\\n", mean(gibbs_t), sd(gibbs_t)))
cat(sprintf("Fang et al.\\ & %.2f & %.2f \\\\\n", mean(fang_t),  sd(fang_t)))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\caption{Average computing time per dataset (seconds) based on 12,000 simulated datasets.}\n")
cat("\\end{table}\n")

sink()
if (sink.number() > 0) sink(type = "output")
cat("Done! NA and time tables appended to simulation/results/tables.tex\n")