# =============================================================================
# sim_6_baseplots.R
# Figures 1-2: averaged baseline survival estimates vs. the truth.
#
# Reads the per-dataset baseline survival curves written by sim_2_ECM.R and
# sim_3_Gibbs_array.R (simulation/results/{ECM,Gibbs}/S0.*.txt), averages them
# within each (alpha,beta) setting, and overlays the true S0(t).
#
# Run with gtcs_probit/ as the working directory.
#   Output: baseplot1.png (ECM), baseplot2.png (Gibbs)
# =============================================================================

grid       <- seq(0, 2, length.out = 100)          # same grid as the estimators
true_S0    <- 1 - pnorm(log(grid / 20) + grid)      # S0(t) = 1 - Phi(g(t)), g(t)=log(t/20)+t
alpha_vals <- c(1, 0.95, 0.90, 0.85)                # row blocks: 1-3000, 3001-6000, ...
titles     <- lapply(c("1", "0.95", "0.90", "0.85"),
                     function(v) bquote(alpha * "=" * beta * "=" * .(v)))

# Load all S0 curves for a method into one (12000 x 100) matrix.
load_S0 <- function(path, n_batch = 24)
  do.call(rbind, lapply(1:n_batch,
          function(i) as.matrix(read.table(paste0(path, "S0.", i, ".txt"), header = TRUE))))

# Mean estimated S0 within each alpha block (3000 datasets each).
mean_curves <- function(S0)
  lapply(seq_along(alpha_vals),
         function(k) colMeans(S0[((k - 1) * 3000 + 1):(k * 3000), ], na.rm = TRUE))

make_plot <- function(curves, file) {
  png(file, width = 1200, height = 800, res = 130)
  par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
  for (k in seq_along(alpha_vals)) {
    plot(grid, true_S0, type = "l", lwd = 2, col = "black",
         ylim = c(0.65, 1), xlab = "t", ylab = expression(S[0](t)), main = titles[k])
    lines(grid, curves[[k]], lwd = 2, lty = 2, col = "red")
    legend("bottomleft", legend = c("True", "Estimated"),
           col = c("black", "red"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.9)
  }
  dev.off()
}

make_plot(mean_curves(load_S0("simulation/results/ECM/")),   "baseplot1.png")
make_plot(mean_curves(load_S0("simulation/results/Gibbs/")), "baseplot2.png")
cat("Wrote baseplot1.png (ECM) and baseplot2.png (Gibbs)\n")
