# =============================================================================
# ECM.R
# ECM algorithm for semiparametric probit model
# under group-tested current status data
#
# Args:
#   Y       : individual-level group test result (length N); or NULL if Yi given
#   Yi      : group-level test result (length n); or NULL if Y given
#   X       : covariate matrix (N x p), WITHOUT intercept
#   c       : testing times (length N)
#   groupID : group membership vector (length N)
#   grid    : grid points for baseline survival function (NULL = auto)
#   n.grid  : number of grid points (used if grid = NULL)
#   m       : number of interior knots
#   init.theta : initial values for theta (length p); default zeros
#   gam0    : initial value for gamma_0; default -3
#   gam     : initial values for gamma_l's; default 0.1
#   alpha   : sensitivity (default 1 = perfect test)
#   beta    : specificity (default 1 = perfect test)
#   order   : spline order (default 3)
#   knots   : user-specified knot sequence (NULL = quantile-based)
#   quantile: logical; use quantile-based knots (default TRUE)
#   maxiter : maximum number of ECM iterations (default 5000)
#   tol     : convergence tolerance (default 1e-5)
# =============================================================================

ECM <- function(Y = NULL, Yi = NULL, X, c, groupID, grid = NULL, n.grid = 100,
                m = 10, init.theta = NULL, gam0 = -3, gam = NULL,
                alpha = 1, beta = 1, order = 3, knots = NULL, quantile = TRUE,
                maxiter = 5000, tol = 1e-5) {

  if (!is.null(knots)) message("Note: m will be redefined based on supplied knots.")
  if (is.null(Y) && is.null(Yi)) stop("Both Y and Yi are NULL. Please provide one.")

  # ---------------------------------------------------------------------------
  # I-spline basis function
  # ---------------------------------------------------------------------------
  Ispline <- function(x, order, knots) {
    m <- length(knots)
    k <- order + 1
    n <- m - 2 + k
    t <- c(rep(knots[1], k), knots[2:(m-1)], rep(knots[m], k))

    yy1 <- matrix(0, n + k - 1, length(x))
    for (l in k:n) {
      yy1[l, ] <- (x >= t[l] & x < t[l+1]) / (t[l+1] - t[l])
    }
    yytem1 <- yy1
    for (ii in 1:order) {
      yytem2 <- matrix(0, n + k - 1 - ii, length(x))
      for (i in (k-ii):n) {
        yytem2[i, ] <- (ii+1) * ((x - t[i]) * yytem1[i, ] +
                                   (t[i+ii+1] - x) * yytem1[i+1, ]) /
          (t[i+ii+1] - t[i]) / ii
      }
      yytem1 <- yytem2
    }
    index   <- sapply(x, function(xi) sum(t <= xi))
    ibases  <- matrix(0, n - 1, length(x))
    for (j in seq_along(x)) {
      for (i in 2:n) {
        if (i < index[j] - order + 1) {
          ibases[i-1, j] <- 1
        } else if (i <= index[j] && i >= index[j] - order + 1) {
          ibases[i-1, j] <- (t[(i+order+1):(index[j]+order+1)] - t[i:index[j]]) %*%
            yytem2[i:index[j], j] / (order + 1)
        }
      }
    }
    return(ibases)
  }

  # ---------------------------------------------------------------------------
  # Knots and grid setup
  # ---------------------------------------------------------------------------
  if (is.null(knots)) {
    if (quantile) {
      knots <- c(0, quantile(c, seq(1/(m+1), 1-1/(m+1), by = 1/(m+1))), max(c) + 0.01)
    } else {
      knots <- seq(0, max(c) + 0.01, length.out = m + 2)
    }
  } else {
    m <- length(knots) - 2
  }
  if (is.null(grid)) grid <- seq(0, max(c) + 0.01, length.out = n.grid)

  # ---------------------------------------------------------------------------
  # Setup
  # ---------------------------------------------------------------------------
  count     <- table(groupID)
  if (is.null(Y)) Y <- rep(Yi, count)
  groupsize <- as.integer(table(groupID))

  X     <- cbind(1, X)
  N     <- nrow(X); p <- ncol(X)
  gamma <- alpha + beta - 1
  k     <- m + order

  I  <- Ispline(c,    order, knots)
  I1 <- Ispline(grid, order, knots)

  # ---------------------------------------------------------------------------
  # Initialise parameters
  # ---------------------------------------------------------------------------
  if (is.null(init.theta)) init.theta <- rep(0, p - 1)
  theta <- c(gam0, init.theta)          # (p)-vector: (gam0, theta_1, ..., theta_{p-1})
  if (is.null(gam)) gam <- rep(0.1, k)

  # ---------------------------------------------------------------------------
  # ECM iterations
  # ---------------------------------------------------------------------------
  iter <- 0; err <- 1

  while (iter < maxiter && err > tol) {
    theta0 <- theta
    u      <- drop(t(gam %*% I)) + drop(X %*% theta)

    m1    <- as.numeric(aggregate(1 - pnorm(u), list(groupID), FUN = prod)[, 2])
    m1rep <- rep(m1, times = groupsize)

    # E(delta_ij | observed data)
    dstar <- ifelse(Y == 1,
                    alpha * pnorm(u) / (alpha - gamma * m1rep),
                    (1 - alpha) * pnorm(u) / (1 - alpha + gamma * m1rep))

    # E(Z_ij | observed data) — numerically stable, matches paper equation
    phi_u <- dnorm(u)
    Phi_u <- pnorm(u)
    eps   <- 1e-300
    zstar <- u + dstar * phi_u / (Phi_u + eps) -
                 (1 - dstar) * phi_u / (1 - Phi_u + eps)

    if (any(is.nan(zstar))) stop("NaN detected in zstar — check model specification.")

    # CM-step 1: update theta (closed form)
    theta <- solve(t(X) %*% X) %*% t(X) %*% (zstar - drop(t(gam %*% I)))

    # CM-step 2: update gamma_l sequentially (non-negativity enforced)
    for (l in 1:k) {
      resid  <- zstar - t(I[-l,]) %*% gam[-l] - X %*% theta
      gam[l] <- max(sum(resid * I[l, ]) / sum(I[l, ]^2), 0)
    }

    err  <- max(abs(theta[-1] - theta0[-1]))
    iter <- iter + 1
    if (iter %% 100 == 0) cat(sprintf("iteration %d complete\n", iter))
  }

  if (iter == maxiter) warning("ECM reached maxiter without converging.")

  # ---------------------------------------------------------------------------
  # Variance estimation via Louis's method
  # ---------------------------------------------------------------------------
  u     <- drop(t(gam %*% I)) + drop(X %*% theta)
  m1    <- as.numeric(aggregate(1 - pnorm(u), list(groupID), FUN = prod)[, 2])
  m1rep <- rep(m1, times = groupsize)

  Edij  <- ifelse(Y == 1,
                  alpha * pnorm(u) / (alpha - gamma * m1rep),
                  (1 - alpha) * pnorm(u) / (1 - alpha + gamma * m1rep))

  phi_u <- dnorm(u); Phi_u <- pnorm(u)
  eps   <- 1e-300
  ratio <- phi_u / (Phi_u * (1 - Phi_u) + eps)

  Vdij <- Edij * (1 - Edij)
  Vzij <- Vdij * ratio^2 + 1 +
          u * phi_u / (1 - Phi_u + eps) -
          (phi_u / (1 - Phi_u + eps))^2 -
          Edij * (u * ratio + (1 - 2 * Phi_u) * ratio^2)

  Q1 <- -t(X) %*% X
  Q2 <- -I %*% t(I)
  Q3 <- -I %*% X
  Q  <- rbind(cbind(Q1, t(Q3)), cbind(Q3, Q2))

  Vlc1 <- t(X) %*% (drop(Vzij) * X)
  Vlc2 <- I %*% (drop(Vzij) * t(I))
  Vlc3 <- I %*% (drop(Vzij) * X)
  Vlc  <- rbind(cbind(Vlc1, t(Vlc3)), cbind(Vlc3, Vlc2))

  info <- (-Q - Vlc)[1:p, 1:p]
  se   <- sqrt(diag(solve(info))[2:p])

  # ---------------------------------------------------------------------------
  # Baseline survival function
  # ---------------------------------------------------------------------------
  gam0_hat <- theta[1]
  St       <- 1 - pnorm(drop(c(gam0_hat, gam) %*% rbind(1, I1)))

  plot(grid, St, xlab = "t", ylab = expression(S[0](t)),
       lwd = 2, col = 2, type = "l")

  # ---------------------------------------------------------------------------
  # Summary table
  # ---------------------------------------------------------------------------
  est    <- theta[-1]
  rnames <- paste0("theta", seq_len(p - 1))
  summary.theta <- data.frame(
    EST  = est,
    SE   = se,
    lb95 = est - 1.96 * se,
    ub95 = est + 1.96 * se,
    row.names = rnames
  )

  return(list(
    iter          = iter,
    gam           = gam,
    gam0          = gam0_hat,
    theta.est     = est,
    se            = se,
    alpha         = alpha,
    beta          = beta,
    grid          = grid,
    St            = St,
    summary.theta = summary.theta
  ))
}