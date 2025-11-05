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
  
  # Construct B-spline basis matrix X tilde: (n + t_lag) × K
  # Each column is a cubic B-spline basis function b_k(t)
  Xtilde <- splineDesign(knots = knots, 
                         x = eval_t, 
                         ord = 4)

  
  # Convolution to get X
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
  
  
  # Second-order difference penalty matrix S
  S <- crossprod(diff(diag(K), differences = 2))
  
  # Return all components needed for model fitting
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







