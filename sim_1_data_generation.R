# =============================================================================
# sim_1_data_generation.R
# Data generation for simulation study
# Chapter 2: Regression Analysis of Group-Tested Current Status Data
#             Under Semiparametric Probit Model
#
# Setting:
#   N = 4000 individuals, n = 1000 groups, pool size J = 4
#   T ~ F(t | x1, x2) = Phi(g(t) + theta1*x1 + theta2*x2)
#   g(t) = log(t/20) + t  (true baseline)
#   x1 ~ N(0, 0.5),  x2 ~ Bernoulli(0.5)
#   C  ~ Uniform(0.1, 2)
#   (theta1, theta2) in {-0.5, 0, 0.5} x {0, 0.5}
#   (alpha, beta) in {(1,1), (0.95,0.95), (0.90,0.90), (0.85,0.85)}
#   => 12 parameter combinations x 4 (alpha,beta) settings = 12,000 datasets
# =============================================================================

# -----------------------------------------------------------------------------
# True baseline function
# -----------------------------------------------------------------------------
g <- function(t) log(t / 20) + t

# -----------------------------------------------------------------------------
# Simulation parameters
# -----------------------------------------------------------------------------
N   <- 4000          # total individuals
n   <- 1000          # number of groups
J   <- N / n         # pool size (= 4)
B   <- 500           # number of simulation replicates per setting (paper uses 500)

theta1_vec <- rep(c(-0.5,  0, 0.5), each = 1000, times = 4)
theta2_vec <- rep(c(0, 0.5),        each =  500, times = 12)
Alpha      <- rep(c(1, 0.95, 0.9, 0.85), each = 3000)

groupID <- rep(1:n, each = J)   # fixed group membership

# -----------------------------------------------------------------------------
# Data generation loop
# -----------------------------------------------------------------------------
set.seed(0)
DATA <- vector("list", 12000)

for (l in 1:12000) {

  alpha <- Alpha[l]
  beta  <- alpha                          # alpha = beta in all settings
  theta <- c(theta1_vec[l], theta2_vec[l])

  # Step 1: Generate covariates and testing times
  X1 <- rnorm(N, mean = 0,   sd = 0.5)
  X2 <- rbinom(N, size = 1,  prob = 0.5)
  C  <- runif(N, min = 0.1,  max = 2)

  # Step 2: Generate individual true disease status
  #   delta_ij = I(T_ij <= C_ij) ~ Bernoulli( Phi(g(C) + theta'x) )
  eta <- g(C) + theta[1] * X1 + theta[2] * X2
  dij <- rbinom(N, size = 1, prob = pnorm(eta))

  # Step 3: True group status  Delta_i = I(sum_j delta_ij > 0)
  Delta <- as.integer(tapply(dij, groupID, sum) > 0)   # n-vector

  # Step 4: Observed (possibly imperfect) group test result Y_i
  Yi <- Delta
  if (alpha != 1 || beta != 1) {
    pos <- Delta == 1
    neg <- Delta == 0
    Yi[pos] <- rbinom(sum(pos), size = 1, prob = alpha)
    Yi[neg] <- rbinom(sum(neg), size = 1, prob = 1 - beta)
  }

  # Store individual-level data frame
  DATA[[l]] <- data.frame(
    X1      = X1,
    X2      = X2,
    C       = C,
    dij     = dij,                    # individual true status
    Delta   = rep(Delta, each = J),   # true group status (individual-level)
    Y       = rep(Yi,    each = J),   # observed group test result (individual-level)
    groupID = groupID
  )

  if (l %% 500 == 0)
    cat(sprintf("theta=(%4.1f, %3.1f)  alpha=beta=%.2f  iter %5d done\n",
                theta[1], theta[2], alpha, l))
}

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
# Average left-censoring rate across all datasets (should be ~10.7%)
mean(sapply(DATA, function(d) mean(d$dij)))

# True baseline survival function
t_seq <- seq(0.1, 2, by = 0.01)
plot(t_seq, 1 - pnorm(g(t_seq)),
     type = "l", lwd = 2,
     xlab = "t", ylab = expression(S[0](t)),
     main = "True Baseline Survival Function")