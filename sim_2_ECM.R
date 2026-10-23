# =============================================================================
# sim_2_ECM.R
# Simulation study using ECM algorithm
# Chapter 2: Regression Analysis of Group-Tested Current Status Data
#             Under Semiparametric Probit Model
#
# Requires:
#   - DATA list from sim_data_generation.R
#   - ECM() from R/ECM.R
#
# Output (saved every 500 iterations):
#   simulation/results/ECM/theta1.{ind}.txt
#   simulation/results/ECM/theta2.{ind}.txt
#   simulation/results/ECM/S0.{ind}.txt
# =============================================================================

source("ECM.R")

# Load simulation datasets
if (!exists("DATA")) {
  if (file.exists("DATA.rds")) {
    cat("Loading DATA from DATA.rds...\n")
    DATA <- readRDS("DATA.rds")
  } else {
    cat("DATA not found. Generating...\n")
    source("sim_1_data_generation.R")
    saveRDS(DATA, "DATA.rds")
    cat("DATA saved to DATA.rds\n")
  }
}

# Create output directory
dir.create("simulation/results/ECM", recursive = TRUE, showWarnings = FALSE)

# Alpha/Beta settings (matches sim_data_generation.R)
Alpha <- rep(c(1, 0.95, 0.9, 0.85), each = 3000)

# Result containers
theta1    <- matrix(0, 12000, 4)   # EST, SE, lb95, ub95
theta2    <- matrix(0, 12000, 4)
S0        <- matrix(0, 12000, 100)
comp.time <- numeric(12000)        # per-dataset computing time (seconds)
na.flag   <- NULL

total.start <- proc.time()

for (IT in 1:12000) {

  Data  <- DATA[[IT]]
  X     <- cbind(Data$X1, Data$X2)
  c     <- Data$C
  alpha <- beta <- Alpha[IT]

  t.start <- proc.time()

  result <- ECM(
    Y       = Data$Y,
    X       = X,
    c       = c,
    groupID = Data$groupID,
    alpha   = alpha,
    beta    = beta,
    grid    = seq(0, 2, length.out = 100),
    m       = 10
  )

  comp.time[IT] <- (proc.time() - t.start)["elapsed"]

  theta1[IT, ] <- as.numeric(result$summary.theta[1, ])
  theta2[IT, ] <- as.numeric(result$summary.theta[2, ])
  S0[IT, ]     <- result$St

  if (is.na(theta1[IT, 2]) | is.na(theta2[IT, 2])) {
    na.flag <- c(na.flag, IT)
    cat(sprintf("NA detected at IT = %d\n", IT))
  }

  cat(sprintf("Data %d complete | theta1=%.3f  theta2=%.3f | time=%.2fs\n",
              IT, theta1[IT, 1], theta2[IT, 1], comp.time[IT]))

  # Save every 500 iterations
  if (IT %% 500 == 0) {
    ind  <- ceiling(IT / 500)
    flag <- (500 * (ind - 1) + 1):(500 * ind)
    write.table(theta1[flag, ],     file = paste0("simulation/results/ECM/theta1.",     ind, ".txt"), row.names = FALSE)
    write.table(theta2[flag, ],     file = paste0("simulation/results/ECM/theta2.",     ind, ".txt"), row.names = FALSE)
    write.table(S0[flag, ],         file = paste0("simulation/results/ECM/S0.",         ind, ".txt"), row.names = FALSE)
    write.table(comp.time[flag],    file = paste0("simulation/results/ECM/comp.time.",  ind, ".txt"), row.names = FALSE)
    cat(sprintf(">>> Saved batch %d | mean time=%.2fs\n", ind, mean(comp.time[flag])))
  }
}

# Total elapsed time
total.time <- (proc.time() - total.start)["elapsed"]
cat(sprintf("\n=== ECM Simulation Complete ===\n"))
cat(sprintf("Total time     : %.1f min\n", total.time / 60))
cat(sprintf("Mean time/data : %.2f sec\n", mean(comp.time)))
cat(sprintf("SD   time/data : %.2f sec\n", sd(comp.time)))

# Save summary
time.summary <- data.frame(
  mean = mean(comp.time),
  sd   = sd(comp.time),
  min  = min(comp.time),
  max  = max(comp.time),
  total_min = total.time / 60
)
write.table(time.summary, file = "simulation/results/ECM/time_summary.txt", row.names = FALSE)

# Save NA flags
if (!is.null(na.flag)) {
  writeLines(as.character(na.flag), "simulation/results/ECM/na_flag.txt")
  cat(sprintf("Total NA cases: %d\n", length(na.flag)))
}