# =============================================================================
# sim_3_Gibbs_array.R
# SLURM job-array version of the Gibbs simulation study.
#
# One array task = one BATCH of consecutive datasets (default 500).
# 12,000 datasets / 500 = 24 batches  ->  --array=1-24
#
# Batch b handles datasets ((b-1)*BATCH_SIZE + 1) : (b*BATCH_SIZE).
# Within a task, datasets run in parallel via a FORK backend (shared memory:
# DATA is NOT copied per worker, unlike the PSOCK version in sim_Gibbs.R).
#
# Output (identical names/format to sim_Gibbs.R, so summarize_results.R works):
#   simulation/results/Gibbs/theta1.{b}.txt     (BATCH_SIZE x 4: mean,sd,q025,q975)
#   simulation/results/Gibbs/theta2.{b}.txt     (BATCH_SIZE x 4)
#   simulation/results/Gibbs/S0.{b}.txt         (BATCH_SIZE x 100)
#   simulation/results/Gibbs/comp.time.{b}.txt  (BATCH_SIZE x 1, per-dataset secs)
#
# SLURM sets:   SLURM_ARRAY_TASK_ID (batch), SLURM_CPUS_PER_TASK (cores).
# Optional env overrides (defaults = paper config; used for local smoke tests):
#   DATA_PATH (DATA.rds)  OUTDIR  BATCH_SIZE (500)  MAXITER (5000)  BURNIN (1000)
# =============================================================================

suppressMessages({ library(doParallel); library(foreach) })

geti <- function(name, default) {
  v <- Sys.getenv(name, unset = NA); if (is.na(v) || v == "") default else v
}

batch     <- as.integer(geti("SLURM_ARRAY_TASK_ID", "1"))
ncores    <- as.integer(geti("SLURM_CPUS_PER_TASK", "1"))
data.path <- geti("DATA_PATH", "DATA.rds")
outdir    <- geti("OUTDIR",    "simulation/results/Gibbs")
bsize     <- as.integer(geti("BATCH_SIZE", "500"))
maxiter   <- as.integer(geti("MAXITER", "5000"))
burnin    <- as.integer(geti("BURNIN",  "1000"))

source("Gibbs.R")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DATA  <- readRDS(data.path)
# alpha=beta settings, must match sim_data_generation.R
Alpha <- rep(c(1, 0.95, 0.9, 0.85), each = 3000)
if (length(Alpha) != length(DATA))
  Alpha <- rep(Alpha, length.out = length(DATA))   # keep local sub-samples working

flag <- ((batch - 1) * bsize + 1):(batch * bsize)
flag <- flag[flag <= length(DATA)]

# --- Memory hardening -------------------------------------------------------
# Keep ONLY this batch's datasets, then drop the rest and free memory. Without
# this, every task holds all ~12,000 datasets (~1.5 GB) in the parent; R's GC
# touches those pages, defeating fork copy-on-write and multiplying resident
# memory by the number of workers -> OOM at high core counts. After subsetting,
# the shared data is ~batch-sized (~60 MB), so worker memory stays flat.
alpha.batch <- Alpha[flag]           # per-dataset alpha, aligned to `flag`
DATA        <- DATA[flag]            # ~1.5 GB  ->  ~60 MB
invisible(gc())
# ---------------------------------------------------------------------------

cat(sprintf("batch %d | datasets %d..%d (n=%d) | %d cores | maxiter=%d burnin=%d\n",
            batch, min(flag), max(flag), length(flag), ncores, maxiter, burnin))

# FORK on Linux/macOS (shared memory); PSOCK fallback elsewhere.
use.fork <- .Platform$OS.type == "unix"
if (use.fork) {
  registerDoParallel(cores = ncores)
} else {
  cl <- makeCluster(ncores); registerDoParallel(cl)
}

t0 <- proc.time()
results <- foreach(
  k = seq_along(flag),
  .packages = c("truncnorm", "mvtnorm"),
  .export   = if (use.fork) NULL else c("DATA", "alpha.batch", "flag")
) %dopar% {
  if (!use.fork) source("Gibbs.R")
  Data  <- DATA[[k]]                  # DATA is already subset to this batch
  ts    <- proc.time()
  res   <- Gibbs(
    Y = Data$Y, X = cbind(Data$X1, Data$X2), c = Data$C, groupID = Data$groupID,
    alpha = alpha.batch[k], beta = alpha.batch[k],
    grid = seq(0, 2, length.out = 100), m = 10,
    maxiter = maxiter, burn.in = burnin
  )
  list(IT = flag[k],
       theta1 = as.numeric(res$summary.theta[1, ]),
       theta2 = as.numeric(res$summary.theta[2, ]),
       St     = res$St,
       ct     = as.numeric((proc.time() - ts)["elapsed"]))
}
if (!use.fork) stopCluster(cl)
elapsed <- (proc.time() - t0)["elapsed"]

# reassemble in dataset order
ord <- order(vapply(results, `[[`, numeric(1), "IT"))
results <- results[ord]
theta1 <- t(vapply(results, `[[`, numeric(4),   "theta1"))
theta2 <- t(vapply(results, `[[`, numeric(4),   "theta2"))
S0     <- t(vapply(results, `[[`, numeric(100), "St"))
ct     <- vapply(results, `[[`, numeric(1),     "ct")

write.table(theta1, file.path(outdir, paste0("theta1.",    batch, ".txt")), row.names = FALSE)
write.table(theta2, file.path(outdir, paste0("theta2.",    batch, ".txt")), row.names = FALSE)
write.table(S0,     file.path(outdir, paste0("S0.",        batch, ".txt")), row.names = FALSE)
write.table(ct,     file.path(outdir, paste0("comp.time.", batch, ".txt")), row.names = FALSE)

na <- sum(!is.finite(theta1[, 1]) | !is.finite(theta2[, 1]))
cat(sprintf("batch %d done | %.1f min wall | %.1fs/dataset (serial-equiv) | NA=%d\n",
            batch, elapsed / 60, mean(ct), na))
