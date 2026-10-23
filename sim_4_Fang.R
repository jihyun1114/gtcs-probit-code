# =============================================================================
# sim_4_Fang.R  -- Fang et al. (2027) method over the simulation datasets.
#
# >>> NOTE <<<
# The original Fang simulation driver was not kept. This script re-runs the Fang
# et al. (2027) estimator on DATA.rds so that the reported Table 1-4 Fang columns
# are exactly reproducible. For a FAIR comparison, Fang uses the SAME fixed spline
# specification as the proposed ECM/Gibbs -- degree 3, m=10 interior knots at
# equal-probability quantiles of the testing times (fang_code/Support_functions_matched.R)
# -- rather than its own AIC/BIC selection. One fit + profile-likelihood SE per
# dataset (~21 s), 24-batch SLURM array like sim_3_Gibbs_array.R.
#
# Output (matches the existing simulation/results/Fang/ files):
#   theta1.{b}.txt, theta2.{b}.txt  (cols: Est, SE, L95, U95; 500 rows)
#   comp.time.{b}.txt               (per-dataset seconds)
# =============================================================================
suppressMessages({ library(doParallel); library(foreach) })
source("fang_code/Support_functions_matched.R")   # Fang estimator with quantile knots (matched to ECM/Gibbs)

geti <- function(name, default) { v <- Sys.getenv(name, unset = NA); if (is.na(v) || v == "") default else v }
batch  <- as.integer(geti("SLURM_ARRAY_TASK_ID", "1"))
ncores <- as.integer(geti("SLURM_CPUS_PER_TASK", "1"))
bsize  <- as.integer(geti("BATCH_SIZE", "500"))
outdir <- geti("OUTDIR", "simulation/results/Fang")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DATA  <- readRDS(geti("DATA_PATH", "DATA.rds"))
Alpha <- rep(c(1, 0.95, 0.9, 0.85), each = 3000)          # alpha = beta, matches sim_1_data_generation.R
flag  <- ((batch - 1) * bsize + 1):(batch * bsize); flag <- flag[flag <= length(DATA)]
alpha.batch <- Alpha[flag]; DATA <- DATA[flag]; invisible(gc())

# Fixed spline settings matched to the proposed ECM/Gibbs -- degree 3, m=10 quantile knots
# (no AIC/BIC selection), for a fair comparison. Override via env if needed.
ORD  <- as.integer(geti("FANG_ORDER", "3"))
NINT <- as.integer(geti("FANG_NINT",  "10"))

registerDoParallel(cores = ncores)
cat(sprintf("batch %d | ncores=%d | %d datasets | getDoParWorkers=%d\n",
            batch, ncores, length(flag), getDoParWorkers()))
one_fang <- function(Data, a) {
  z <- cbind(Data$X1, Data$X2); xt <- Data$C; groups <- Data$groupID
  upool.id <- unique(groups); g.n <- as.integer(table(groups))
  Y <- as.integer(tapply(Data$Y, groups, function(v) v[1]))   # pool outcome (one per pool)
  nu <- rep(a, length(Y)); omega <- rep(a, length(Y)); p <- ncol(z)
  fit <- probit_GT(xt=xt, groups=groups, g.n=g.n, n.int=NINT, order=ORD,
                   z=z, upool.id=upool.id, Ys=Y, nu=nu, omega=omega)
  h <- 5*length(upool.id)^(-0.5); E <- h*diag(p); Me <- matrix(0, p, length(upool.id))
  for (i in 1:p) {
    pl1 <- prolik(xt=xt,groups=groups,g.n=g.n,n.int=NINT,order=ORD,z=z,upool.id=upool.id,Ys=Y,be=fit$beta+E[,i],g.est=fit$g,nu=nu,omega=omega)
    pl2 <- prolik(xt=xt,groups=groups,g.n=g.n,n.int=NINT,order=ORD,z=z,upool.id=upool.id,Ys=Y,be=fit$beta-E[,i],g.est=fit$g,nu=nu,omega=omega)
    Me[i,] <- (pl1-pl2)/(2*h)
  }
  se <- sqrt(diag(solve(Me %*% t(Me)))); b <- as.vector(fit$beta)
  list(theta1 = c(b[1], se[1], b[1]-1.96*se[1], b[1]+1.96*se[1]),
       theta2 = c(b[2], se[2], b[2]-1.96*se[2], b[2]+1.96*se[2]))
}

t0 <- proc.time()
res <- foreach(k = seq_along(flag), .packages = character(0)) %dopar% {
  ts <- proc.time()
  out <- tryCatch(one_fang(DATA[[k]], alpha.batch[k]), error = function(e) NULL)
  list(theta1 = if (is.null(out)) rep(NA,4) else out$theta1,
       theta2 = if (is.null(out)) rep(NA,4) else out$theta2,
       ct = as.numeric((proc.time()-ts)["elapsed"]))
}
theta1 <- t(sapply(res, `[[`, "theta1")); theta2 <- t(sapply(res, `[[`, "theta2")); ct <- sapply(res, `[[`, "ct")
write.table(theta1, file.path(outdir, paste0("theta1.",    batch, ".txt")), row.names = FALSE)
write.table(theta2, file.path(outdir, paste0("theta2.",    batch, ".txt")), row.names = FALSE)
write.table(ct,     file.path(outdir, paste0("comp.time.", batch, ".txt")), row.names = FALSE)
cat(sprintf("batch %d done | %.1f min | NA=%d\n", batch, (proc.time()-t0)["elapsed"]/60, sum(is.na(theta1[,1]))))
