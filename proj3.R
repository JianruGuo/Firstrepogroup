# Contribution of each group member:
# Zixuan Qiu s2777279: create a function 'build_designs()' to evaluate X, X tilde and S
# Jianru Guo s2806788: define likelihood and gradient functions, fit model using BFGS optimization and select optimal smoothness
# Fengyu Shen s2798951: refit the model at the best gamma, ran 200 bootstraps for 95% CIs, and plotted fitted deaths and infection curve

# Code Structure Overview
# 1. build_designs() builds spline, convolution matrices and penalty matrix.
# 2. pen_neg_loglik() and pen_neg_loglik_grad() define the penalized likelihood and gradient.
# 3. fit_one_lambda() fits the Poisson model using BFGS optimization.
# 4. A grid search selects the best λ by minimizing BIC.
# 5. The model is refitted at the best λ, and bootstrap resampling (B = 200) computes 95% confidence bands for fitted f(t).

# GitHub Link: https://github.com/JianruGuo/Firstrepogroup.git

# getwd()
# setwd("C:/Users/cheryl/Desktop/PG-Study/Extended Statistical Programming/Assignment/project3")
# Function to build design matrices for COVID-19 death deconvolution model
# This implements a Poisson GLM where deaths are modeled as a convolution
# of infection trajectory f(t) with infection to death duration distribution pi(j)

library(splines)

build_designs <- function(data_path = "engcov.txt",
                          K = 80,         # number of spline basis functions
                          t_lag = 30,     # number of past days
                          edur = 3.151,   # lognormal meanlog for infection to death duration
                          sdur = 0.469) { # lognormal sdlog for infection to death duration
  
  engcov <- read.table(data_path, header = TRUE)
  t_vec <- engcov$julian # Julian day numbers
  y <- engcov$deaths # observed daily deaths
  n <- length(t_vec) # number of observation days
  
  # the model probability of each fatal disease duration 
  # from 1 day up to the maximum of 80 days
  # i.e. pi(j)
  d <- 1:80
  pd <- dlnorm(d, meanlog = edur, sdlog = sdur)
  pd <- pd/sum(pd) # normalization
  
  # prepare evaluation points for spline basis functions
  # from t1-30 to tn
  eval_t <- seq(from = t_vec[1] - t_lag, to = t_vec[n], by = 1)
  
  # create K+4 evenly spaced knots for B-spline basis
  # middle K-2 knots cover the evaluation interval
  # This means 3 knots extend beyond each boundary
  t_min <- min(eval_t)
  t_max <- max(eval_t)
  t_range <- t_max - t_min
  # K-2 interior knots divide the range into K-3 intervals
  spacing <- t_range / (K - 3)
  # create K+4 knots
  knots <- seq(from = t_min - 3 * spacing, 
               to = t_max + 3 * spacing, 
               length.out = K + 4)
  
  # construct B-spline basis matrix X tilde: (n + t_lag) × K
  # each column is a cubic B-spline basis function b_k(t)
  Xtilde <- splineDesign(knots = knots, 
                         x = eval_t, 
                         ord = 4)
  
  # convolution to get X
  # initialize a (n x K) matrix
  # row: days, column: spline basis function
  X <- matrix(0, nrow = n, ncol = K)
  
  for (i in 1:n) {
    jmax <- min(29 + i, 80)
    
    for (j in 1:jmax) {
      row_idx <- (t_lag + i) - j
      X[i, ] <- X[i, ] + pd[j] * Xtilde[row_idx, ]
    }
  }
  
  # second-order difference penalty matrix S
  S <- crossprod(diff(diag(K), differences = 2))
  
  # return all components needed for model fitting
  list(
    Xtilde = Xtilde,
    X = X,
    S = S,
    y = y,
    t = t_vec,
    pd = pd,
    eval_t = eval_t
  )
}

# construct penalized negative log likelihood and its gradient, then test the derivative function by finite differencing
pen_neg_loglik <- function(gamma, X, S, y, lambda) {
  beta <- exp(gamma)                # β = exp(γ) 
  mu   <- as.vector(X %*% beta)     # μ = X β
  
  # avoid log(0)
  mu[mu < 1e-8] <- 1e-8
  
  # construct Poisson negative log-likelihood for minimization(don't need the constant term)
  nll <- sum(mu - y * log(mu))
  
  # construct penalty
  pen <- 0.5 * lambda * as.numeric(t(beta) %*% S %*% beta)
  
  nll + pen
}

# construct gradient of the pen_neg_loglik function
pen_neg_loglik_grad <- function(gamma, X, S, y, lambda) {
  beta <- exp(gamma)                
  mu   <- as.vector(X %*% beta)     
  mu[mu < 1e-8] <- 1e-8           
  
  # gradient of the term about negative log likelihood
  w <- 1 - y / mu
  X_beta <- X * rep(beta, each = nrow(X))  # dimension(n×K):using j-th column in X to multiply each corresponding β
  g_ll   <- as.vector(t(X_beta) %*% w)     # dimension(K×1)
  
  # gradient of the term about penalty
  Sbeta   <- as.vector(S %*% beta)
  g_pen <- lambda * beta * Sbeta
  
  g_ll + g_pen
}

# test the derivative function by finite differencing
des <- build_designs("engcov.txt", K = 80)

X <- des$X
S <- des$S
y <- des$y
t_vec  <- des$t        # length n, Julian dates for deaths
Xtilde <- des$Xtilde   # model matrix for f(t)
eval_t <- des$eval_t

n <- length(y)
K <- ncol(X)
lambda <- 5e-5          # smooth parameter λ specified by the question
th0<- rep(0, ncol(X))  # initialize γ ( β = exp(0) = 1)
fd  <- th0     #initialize              

nll0 <- pen_neg_loglik(th0, X = X, S = S, y = y, lambda = lambda)
eps  <- 1e-6                 # finite difference step

for (i in 1:length(th0)) {   # loop over parameters
  th1 <- th0
  th1[i] <- th1[i] + eps    # increase th0[i] by eps
  nll1 <- pen_neg_loglik(th1, X = X, S = S, y = y, lambda = lambda)
  fd[i] <- (nll1 - nll0) / eps   
}

g_analytic <- pen_neg_loglik_grad(th0, X = X, S = S, y = y, lambda = lambda)

cat("Max abs difference:", max(abs(fd - g_analytic)), "\n")

cbind(
  fd        = fd[1:10],
  analytic  = g_analytic[1:10],
  diff      = fd[1:10] - g_analytic[1:10]
)

fit_one_lambda <- function(lambda, th0, X, S, y, Xtilde){

 fit <- optim(
   par     = th0,                
   fn      = pen_neg_loglik,        # penalized negative log-likelihood function
   gr      = pen_neg_loglik_grad,   # gradient of the penalized NLL
   method  = "BFGS",                # quasi-Newton optimization method
   X = X, S = S, y = y, lambda = lambda,  
   control = list(maxit = 1000, trace = 1) # allow up to 1000 iterations and show progress
 )

 if (fit$convergence != 0) {
   warning("optim did not fully converge for lambda =", lambda)
 }
 fit$convergence   
 fit$value         # final minimized penalized negative log likelihood

 # compute fitted values for deaths and infections 
 gamma_hat <- fit$par          # optimized γ values
 beta_hat  <- exp(gamma_hat)   # convert to β = exp(γ), ensures positivity

 # predicted deaths: μ̂ = X * β
 mu_hat <- as.vector(X %*% beta_hat)

 # fitted infection curve
 f_hat <- as.vector(Xtilde %*% beta_hat)
 
 loglik <- sum(y * log(mu_hat)-mu_hat)
 
 # compute Hessians H0 and H_lambda wrt beta
 w <- y / (mu_hat^2)       # w
 # diag(W) %*% X = X * w (row-wise multiply)
 WX  <- X * w             
 H0  <- crossprod(X, WX)   
 Hlam <- H0 + lambda * S   # H_lambda = H0 + lambda S
 
 # EDF = trace(H_lambda^{-1} H0)
 V   <- solve(Hlam, H0)    # Hlam %*% V = H0
 EDF <- sum(diag(V))
 
 # BIC(λ) = -2 l(β_hat) + log(n) * EDF
 BIC_val <- -2 * loglik + log(n) * EDF
 
 list(
   gamma_hat = gamma_hat,
   beta_hat  = beta_hat,
   mu_hat    = mu_hat,
   f_hat     = f_hat,
   loglik    = loglik,
   EDF       = EDF,
   BIC       = BIC_val,
   value     = fit$value   # penalized NLL
 )
}

# visualization
# plot actual and fitted deaths against time
fit_T3 <- fit_one_lambda(lambda, th0, X, S, y, Xtilde)

gamma_hat <- fit_T3$gamma_hat
beta_hat  <- fit_T3$beta_hat
mu_hat    <- fit_T3$mu_hat
f_hat     <- fit_T3$f_hat

# two panels, tidy margins
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1)) 

# top panel: actual vs fitted deaths
plot(t_vec, y, type = "p", pch = 16, col = "grey40",
     xlab = "Julian day", ylab = "Daily deaths",
     main = paste("Actual vs fitted deaths (lambda =", lambda, ")"))
# add fitted curve 
lines(t_vec, mu_hat, col = "red", lwd = 2)
legend("topleft", legend = c("Observed", "Fitted"), 
       col = c("grey40", "red"), pch = c(16, NA), lty = c(NA, 1), bty = "n")

# bottom panel: fitted infection curve 
plot(eval_t, f_hat, type = "l", lwd = 2, col = "blue",
     xlab = "Julian day", ylab = "Infection intensity f(t)",
     main = "Fitted infection curve f(t)")
# vertical lines mark the range where deaths are observed
abline(v = range(t_vec), lty = 3, col = "grey70")

# define lambda grid
log_lambda_grid <- seq(-13, -7, length.out = 50)
lambda_grid     <- exp(log_lambda_grid)

# storage
BIC_vals   <- numeric(length(lambda_grid))
EDF_vals   <- numeric(length(lambda_grid))
loglik_vals <- numeric(length(lambda_grid))
gamma_store <- matrix(NA, nrow = length(lambda_grid), ncol = K)

# initial gamma: reuse from previous λ
th0 <- rep(0, K)  

for (i in seq_along(lambda_grid)) {
  lambda <- lambda_grid[i]
  cat("Fitting for lambda =", lambda, " (i =", i, "of", length(lambda_grid), ")\n")
  
  fit_i <- fit_one_lambda(lambda, th0, X, S, y, Xtilde)
  
  BIC_vals[i]    <- fit_i$BIC
  EDF_vals[i]    <- fit_i$EDF
  loglik_vals[i] <- fit_i$loglik
  gamma_store[i, ] <- fit_i$gamma_hat
  
  th0 <- fit_i$gamma_hat 
}

# find the lambda that minimums BIC
best_idx <- which.min(BIC_vals)
lambda_best <- lambda_grid[best_idx]
log_lambda_best <- log_lambda_grid[best_idx]

cat("Best lambda =", lambda_best, " (log lambda =", log_lambda_best, ")\n")

# refit with the best lambda
best_fit <- fit_one_lambda(lambda_best,
                           th0 = gamma_store[best_idx, ],
                           X = X, S = S, y = y, Xtilde = Xtilde)
# extract fitted components
gamma_hat_best <- best_fit$gamma_hat
beta_hat_best  <- best_fit$beta_hat
mu_hat_best    <- best_fit$mu_hat
f_hat_best     <- best_fit$f_hat

cat("Refitted model at best lambda =", format(lambda_best, scientific = TRUE), "\n")

# define bootstrap function 
bootstrap_f <- function(B_times, y, X, S, lambda, Xtilde, gamma_init) {
  n <- length(y) # sample size
  K <- ncol(X) # number of basis functions
  m <- nrow(Xtilde) # number of evaluation points
  f_mat <- matrix(NA_real_, nrow = m, ncol = B_times) # storage matrix
  gamma0 <- gamma_init # initialize gamma
  
  for (b in 1:B_times) {
    # nonparametric bootstrap: sample indices with replacement
    w <- tabulate(sample.int(n, n, replace = TRUE), nbins = n)
    # fit penalized model using BFGS optimization
    optb <- optim(par = gamma0,
                  fn = pen_neg_loglik,
                  gr = pen_neg_loglik_grad,
                  method = "BFGS",
                  X = X, S = S, y = y, lambda = lambda,
                  control = list(maxit = 1000, reltol = 1e-8))
    # update initial gamma for next iteration to improve convergence
    gamma0 <- optb$par
    beta_b <- exp(gamma0)
    f_mat[, b] <- as.vector(Xtilde %*% beta_b)
    # print progress every 10 replicates
    if (b %% 10 == 0)
      cat("Bootstrap replicate", b, "of", B_times, "done\n")
  }
  f_mat
}

# set random seed for reproducibility
set.seed(123)
B_reps <- 200   
cat("Running bootstrap with", B_reps, "replicates...\n")

f_boot <- bootstrap_f(B_times = B_reps,
                      y = y, X = X, S = S,
                      lambda = lambda_best, # keep λ fixed
                      Xtilde = Xtilde,
                      gamma_init = gamma_hat_best)

# compute 95% confidence intervals
f_lo <- apply(f_boot, 1, quantile, probs = 0.025, na.rm = TRUE)
f_hi <- apply(f_boot, 1, quantile, probs = 0.975, na.rm = TRUE)

# final visualization
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# top panel: observed vs fitted deaths
plot(t_vec, y, pch = 16, col = "grey40",
     main = sprintf("Final fitted deaths (lambda = %.2e)", lambda_best),
     xlab = "Julian day (2020)", ylab = "Daily NHS deaths")
lines(t_vec, mu_hat_best, col = "red", lwd = 2)
legend("topleft", legend = c("Observed", "Fitted"),
       col = c("grey40", "red"), pch = c(16, NA), lty = c(NA, 1), bty = "n")

# bottom panel: infection curve f(t) with 95% bootstrap band
ok <- is.finite(f_lo) & is.finite(f_hi) & is.finite(f_hat_best)
yr <- range(c(f_lo[ok], f_hi[ok], f_hat_best[ok]), na.rm = TRUE)
# create empty plot
plot(eval_t[ok], f_hat_best[ok], type = "n",
     ylim = yr,
     main = "Estimated infection curve f(t) with 95% bootstrap interval",
     xlab = "Julian day (2020)", ylab = "Infection intensity f(t)")
# draw 95% CI band (light blue polygon)
polygon(c(eval_t[ok], rev(eval_t[ok])),
        c(f_lo[ok],   rev(f_hi[ok])),
        col = grDevices::adjustcolor("dodgerblue", alpha.f = 0.28),
        border = NA)
# draw fitted mean and boundary lines
lines(eval_t[ok], f_hat_best[ok], col = "blue", lwd = 2)
lines(eval_t[ok], f_hi[ok],      col = "blue", lty = 2, lwd = 1.5)
lines(eval_t[ok], f_lo[ok],      col = "blue", lty = 2, lwd = 1.5)
# add vertical grey lines marking observed data range
abline(v = range(t_vec), lty = 3, col = "grey70")

# check CI width
range(f_hi - f_lo, na.rm = TRUE)
#The 95% bootstrap interval is extremely narrow (<0.1 across all t), 
#which shows that the estimated infection trajectory is very stable and robust.
#The infection peak occurs around Julian day 65–70, about 3–4 weeks before the death peak,
#consistent with the expected infection-to-death delay.


















