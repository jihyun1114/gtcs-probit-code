# =============================================================================
# Gibbs.R
# Gibbs sampler for semiparametric probit model
# under group-tested current status data
#
# Args:
#   Y          : individual-level group test result (length N); or NULL if Yi given
#   Yi         : group-level test result (length n); or NULL if Y given
#   X          : covariate matrix (N x p), WITHOUT intercept
#   c          : testing times (length N)
#   groupID    : group membership vector (length N)
#   grid       : grid points for baseline survival function (NULL = auto)
#   n.grid     : number of grid points (used if grid = NULL)
#   m          : number of interior knots
#   init.theta : initial values for theta (length p); default zeros
#   eta        : initial value for hyperparameter eta; default 1
#   gam0       : initial value for gamma_0; default -3
#   gam        : initial values for gamma_l's; default 0.1
#   theta0     : prior mean for theta; default zero vector
#   Sigma0     : prior covariance for theta; default unit information matrix
#   m0, v0     : prior mean and precision for gam0
#   a0, b0     : hyperprior shape and rate for eta ~ Ga(a0, b0)
#   alpha      : sensitivity (default 1 = perfect test)
#   beta       : specificity (default 1 = perfect test)
#   order      : spline order (default 3)
#   knots      : user-specified knot sequence (NULL = quantile-based)
#   quantile   : logical; use quantile-based knots (default TRUE)
#   maxiter    : total number of Gibbs iterations (default 5000)
#   burn.in    : number of burn-in iterations (default 2000)
# =============================================================================

Gibbs = function(
Y = NULL, Yi = NULL, X, c, groupID, grid = NULL, n.grid = 100, m = 10,
                 init.theta = numeric(ncol(X)), eta = 1, gam0 = -3, gam = NULL,
                 theta0 = NULL, Sigma0 = NULL, m0 = -3, v0 = 1, a0 = 1, b0 = 1,
                 alpha = 1, beta = 1, order = 3, knots = NULL, quantile = TRUE,
                 maxiter = 5000, burn.in = 2000) {

  if(is.null(knots)==F) {print("m will be redefined when knots are specified")}
  if(is.null(Y) && is.null(Yi)) stop("Warning: Both Y and Yi are not provided.")

  Ispline<-function(x,order,knots){
    m=length(knots)
    k=order+1
    n=m-2+k
    t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m])
    yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
    for (l in k:n){
      yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
    }
    yytem1=yy1
    for (ii in 1:order){
      yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
      for (i in (k-ii):n){
        yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
      }
      yytem1=yytem2
    }
    index=rep(0,length(x))
    for (i in 1:length(x)){
      index[i]=sum(t<=x[i])
    }
    ibases=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))
    if (order==1){
      for (i in 2:n){
        ibases[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
      }
    }else{
      for (j in 1:length(x)){
        for (i in 2:n){
          if (i<(index[j]-order+1)){
            ibases[i-1,j]=1
          }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
            ibases[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
          }else{
            ibases[i-1,j]=0
          }
        }
      }
    }
    return(ibases)
  }

  require(truncnorm); require(mvtnorm)

  if (is.null(knots)) {
    if (quantile == TRUE) {knots = c(0, quantile(c, seq(0.1, 0.9, length.out = m)), max(c) + 0.01)
    } else{knots = seq(0, max(c)+ 0.01, length.out = m + 2)}
  } else{m = length(knots) - 2}
  if (is.null(grid)) {grid = seq(0, max(c)+0.01, length.out = n.grid)}
  n.knots = length(knots); n.grid = length(grid)

  groupsize = table(groupID); if (is.null(Y)) {Y = rep(Yi, groupsize)}
  N = nrow(X); p = ncol(X); gamma = alpha + beta - 1
  I = Ispline(c, order, knots = knots); I1 = Ispline(grid, order, knots = knots)
  dij = numeric(N); zij = dij

  if(is.null(theta0)) {theta0 = numeric(ncol(X))}
  if(is.null(Sigma0)) {Sigma0 = nrow(X)*solve(t(X)%*%X)}
  k = m + order

  theta = init.theta; gam = rep(0.1,k)

  theta.mat = matrix(0, maxiter, k + p + 1)
  g.mat = matrix(0, maxiter, n.grid)

  idx.by.pool = split(seq_len(N), groupID)   # individuals grouped by pool (fixed)

  for (iter in 1:maxiter) {
    u = t(gam0 + gam %*% I) + X %*% theta
    Yi = aggregate(Y, list(groupID), FUN = mean)[, 2]
    Di = Yi

    m1 = aggregate(1 - pnorm(u), list(groupID), FUN = prod)[, 2]
    prY1 = alpha - gamma * m1; prD1 = 1 - m1
    prY0 = 1 - prY1; prD0 = 1 - prD1
    prob1 = alpha * prD1 / prY1; prob2 = 1 - beta * prD0 / prY0
    prob1[prob1 > 1] = 1; prob2[prob2 > 1] = 1

    Di[Yi == 1] = rbinom(sum(Yi==1), 1, prob1[Yi == 1])
    Di[Yi == 0] = rbinom(sum(Yi==0), 1, prob2[Yi == 0])
    D = rep(Di, times = groupsize)

    ## Sample delta_ij | Delta_i : independent Bernoulli(pi_ij), pi_ij = Phi(u_ij),
    ## restricted to { sum_j delta_ij >= 1 }.  Drawn sequentially (rejection-free).
    pi.ind = pnorm(u)                        # pi_ij = Phi(u_ij)   (length N)
    dij = numeric(N)
    for (i in which(Di == 1)) {              # only pools with Delta_i = 1
      idx = idx.by.pool[[i]]
      pp  = pi.ind[idx]                       # pool member probs (local; NOT ncol(X) p)
      Ji  = length(pp)
      S   = rev(cumprod(rev(1 - pp)))         # S[j] = prod_{j'=j}^{Ji} (1 - pp[j'])
      d   = integer(Ji); hasPos = FALSE
      for (j in seq_len(Ji)) {
        prob = if (hasPos) pp[j] else {
          denom = 1 - S[j]                   # 1 - prod_{j'>=j}(1 - pp)
          if (denom > 0) pp[j] / denom else 1
        }
        d[j] = rbinom(1L, 1L, min(prob, 1))
        if (d[j] == 1L) hasPos = TRUE
      }
      dij[idx] = d
    }

    zij[dij == 0] = rtruncnorm(sum(dij==0),-Inf,0,u[dij==0],1)
    if (any(dij == 1)) {zij[dij == 1] = rtruncnorm(sum(dij==1),0, Inf,u[dij==1],1)}

    W0 = v0 + N; E0 = (v0 * m0 + sum(unlist(zij) - t(I) %*% gam - X %*% theta)) / W0
    gam0 = rnorm(1, E0, sqrt(1 / W0))

    W = rowSums(I ^ 2)
    for (l in 1:length(W)) {
      if (W[l] == 0) {
        gam[l] = rexp(1, eta)
      } else{
        El = (sum(I[l,] * (zij - gam0 - t(I[-l,]) %*% gam[-l] - X %*% theta)) - eta) / W[l]
        gam[l] = rtruncnorm(1, 0, Inf, El, sqrt(1 / W[l]))
      }
    }

    Sigmatilde = solve(solve(Sigma0) + t(X) %*% X)
    thetatilde = c(Sigmatilde %*% (solve(Sigma0) %*% theta0 + t(X) %*% (zij - gam0 - t(gam %*% I))))
    theta = c(rmvnorm(1, thetatilde, Sigmatilde))

    eta = rgamma(1, a0 + k, b0 + sum(gam))

    theta.mat[iter, 1:p] = theta
    theta.mat[iter, p + 1] = gam0
    theta.mat[iter,-(1:(p + 1))] = gam
    g.mat[iter,] = t(gam0 + gam %*% I1)

    if (iter %% 100 == 0) {cat(paste("iteration", iter, "complete\n"))}
  }

  St = colMeans(1 - pnorm(g.mat[-(1:burn.in),]))
  plot(grid, St, xlab = "t", ylab = expression(S[0](t)), lwd = 2, col = 2, type = "l")

  mean  = colMeans(theta.mat[-(1:burn.in), 1:p])
  asd   = apply(theta.mat[-(1:burn.in), 1:p], 2, sd)
  q.025 = apply(theta.mat[-(1:burn.in), 1:p], 2, function(x) {quantile(x, 0.025)})
  q.975 = apply(theta.mat[-(1:burn.in), 1:p], 2, function(x) {quantile(x, 0.975)})
  rname = NULL
  for (ii in 1:p) {rname = c(rname, paste0("theta", ii))}
  summary.theta = as.data.frame(cbind(mean, asd, q.025, q.975), row.names = rname)

  return(list(theta.mat = theta.mat, theta.mean = mean, theta.q = c(q.025, q.975), alpha = alpha, beta = beta,
      g.mat = g.mat, grid = grid, St = St, summary.theta = summary.theta))
}