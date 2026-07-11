############## packages ###########
require(mixtools)
require(Matrix)
library(splines)
library(locfit)
library(cluster)
library(igraph)
library(genlasso)
library(abind)
############################################################
## Adaptive-Lasso Trajectory Merging (K > 2)
## + Merge-aware B-spline Gamma Fitting
## - Adaptive lasso on LEVELS via weighted soft-thresholding
## - Penalized spline pilot for weights
## - Produces sub_mat compatible with gamma_spline
############################################################

suppressPackageStartupMessages({
  library(splines)
  library(igraph)
})

############################################################
## 0. Utilities
############################################################
make_time_grid <- function(m, range = c(0, 1)) {
  seq(range[1], range[2], length.out = m)
}

############################################################
## 1. Penalized spline pilot (P-spline style)
############################################################
fit_pspline <- function(y, t, df_spline = 10, lambda_pilot = 1e-2,
                        degree = 3, penalty_order = 2) {
  
  B <- bs(t, df = df_spline, degree = degree, intercept = TRUE)
  q <- ncol(B)
  
  D_theta <- diff(diag(q), differences = penalty_order)
  Omega <- t(D_theta) %*% D_theta
  
  theta_hat <- solve(crossprod(B) + lambda_pilot * Omega,
                     crossprod(B, y))
  as.vector(B %*% theta_hat)
}

############################################################
## 2. Adaptive-lasso weights on LEVELS (offline)
############################################################
weights_offline_level <- function(t, d_tilde,
                                  gamma_adapt = 1, eps = 1e-6, h = NULL) {
  
  m <- length(d_tilde)
  if (is.null(h)) {
    dt <- median(diff(t))
    h <- 2 * dt
  }
  
  S <- numeric(m)
  for (j in seq_len(m)) {
    idx <- which(abs(t - t[j]) <= h)
    S[j] <- mean(abs(d_tilde[idx]))
  }
  
  1 / (S^gamma_adapt + eps)
}

############################################################
## 3. Adaptive lasso on LEVELS (closed form)
############################################################
soft_threshold_weighted <- function(d, w, lambda) {
  sign(d) * pmax(abs(d) - lambda * w, 0)
}

############################################################
## 4. Adaptive-lasso merge detection for K trajectories
############################################################
detect_merges_adalasso_level <- function(
    gamma, time,
    lambda_level = 0.2,
    tol_zero = 1e-3,
    df_spline = 8,
    lambda_pilot = 1e-2,
    gamma_adapt = 1,
    eps = 1e-6,
    degree = 3,
    penalty_order = 2,
    h = NULL
) {
  
  K <- dim(gamma)[1]
  m <- dim(gamma)[2]
  p <- dim(gamma)[3]
  
  sub_mat <- matrix(NA_integer_, nrow = K, ncol = m)
  
  ## precompute pairwise delta_hat(t)
  pair_delta <- list()
  
  for (k in 1:(K - 1)) {
    for (l in (k + 1):K) {
      
      # pairwise distance trajectory
      d <- sqrt(rowSums((gamma[k, , ] - gamma[l, , ])^2))
      
      # spline pilot for weights
      d_tilde <- fit_pspline(
        y = d, t = time,
        df_spline = df_spline,
        lambda_pilot = lambda_pilot,
        degree = degree,
        penalty_order = penalty_order
      )
      
      # adaptive weights
      w <- weights_offline_level(
        t = time, d_tilde = d_tilde,
        gamma_adapt = gamma_adapt, eps = eps, h = h
      )
      
      # adaptive lasso on levels
      delta_hat <- soft_threshold_weighted(d, w, lambda = lambda_level)
      
      pair_delta[[paste(k, l, sep = "_")]] <- delta_hat
    }
  }
  
  ## build merge graph at each time
  for (tt in 1:m) {
    adj <- matrix(0, K, K)
    for (k in 1:(K - 1)) {
      for (l in (k + 1):K) {
        key <- paste(k, l, sep = "_")
        if (abs(pair_delta[[key]][tt]) <= tol_zero) {
          adj[k, l] <- adj[l, k] <- 1
        }
      }
    }
    g <- graph_from_adjacency_matrix(adj, mode = "undirected")
    sub_mat[, tt] <- components(g)$membership
  }
  
  sub_mat
}

############################################################
## 5. Constant-merge segments & boundary knots
############################################################
get_constant_merge_segments <- function(sub_mat) {
  m <- ncol(sub_mat)
  breaks <- c(1)
  for (t in 2:m) {
    if (!all(sub_mat[, t] == sub_mat[, t - 1])) breaks <- c(breaks, t)
  }
  breaks <- c(breaks, m + 1)
  lapply(seq_len(length(breaks) - 1), function(i) {
    seq(breaks[i], breaks[i + 1] - 1)
  })
}

segment_knots_from_submat <- function(
    sub_mat,
    time,
    degree = 3,
    min_knots_per_segment = NULL,
    global_df = NULL
) {
  m <- length(time)
  
  # default: enough knots to identify degree-d spline
  if (is.null(min_knots_per_segment)) {
    min_knots_per_segment <- max(1, degree - 1)
  }
  
  segments <- get_constant_merge_segments(sub_mat)
  
  knots <- c()
  
  for (seg in segments) {
    
    t_seg <- time[seg]
    seg_len <- length(t_seg)
    
    # always add boundary knots (merge points)
    knots <- c(knots, t_seg[1], t_seg[length(t_seg)])
    
    # determine number of interior knots for this segment
    if (!is.null(global_df)) {
      # proportional allocation
      n_int <- max(
        min_knots_per_segment,
        round(global_df * seg_len / m)
      )
    } else {
      # minimal safe choice
      n_int <- min_knots_per_segment
    }
    
    # add interior knots if segment long enough
    if (seg_len >= (degree + 2) && n_int > 0) {
      interior_idx <- round(
        seq(2, seg_len - 1, length.out = n_int + 2)
      )[-c(1, n_int + 2)]
      
      knots <- c(knots, t_seg[interior_idx])
    }
  }
  
  # remove boundary knots (bs() adds them automatically)
  knots <- sort(unique(knots))
  knots <- knots[knots > min(time) & knots < max(time)]
  
  return(knots)
}

############################################################
## 6. Merge-aware B-spline fitting (shared coefficients)
############################################################
fit_piecewise_splines <- function(
    x, cl_mat, gamma_raw, sub_mat, time,
    degree = 8, lambda_ridge = 1e-6
) {
  
  K <- dim(gamma_raw)[1]
  m <- dim(gamma_raw)[2]
  p <- dim(gamma_raw)[3]
  
  gamma_new <- gamma_raw
  knots_all <- segment_knots_from_submat(sub_mat, time)
  segments <- get_constant_merge_segments(sub_mat)
  
  for (seg in segments) {
    
    t_idx <- seg
    t_vals <- time[t_idx]
    groups <- split(1:K, sub_mat[, t_idx[1]])
    
    for (grp in groups) {
      
      y_pool <- matrix(NA_real_, nrow = length(t_idx), ncol = p)
      
      for (ii in seq_along(t_idx)) {
        tt <- t_idx[ii]
        ids <- which(cl_mat[, tt] %in% grp)
        if (length(ids) == 0) next
        y_pool[ii, ] <- colMeans(x[ids, tt, ], na.rm = TRUE)
      }
      
      for (j in 1:p) {
        y <- y_pool[, j]
        ok <- is.finite(y)
        if (sum(ok) < (degree + 2)) next
        
        B_ok <- bs(t_vals[ok], knots = knots_all,
                   degree = degree, intercept = TRUE)
        B_all <- bs(t_vals, knots = knots_all,
                    degree = degree, intercept = TRUE)
        
        beta <- solve(
          crossprod(B_ok) + lambda_ridge * diag(ncol(B_ok)),
          crossprod(B_ok, y[ok])
        )
        
        fitted <- as.vector(B_all %*% beta)
        for (k in grp) gamma_new[k, t_idx, j] <- fitted
      }
    }
  }
  
  gamma_new
}

############################################################
## 7. Main gamma_spline with adaptive-lasso merging
############################################################
gamma_spline <- function(
    x, time,
    K = 3,
    m.iter = 50,
    eps = 1e-4,
    lambda_level = 0.5,
    tol_zero = 1e-3,
    df_spline_merge = 10,
    lambda_pilot = 1e-2,
    gamma_adapt = 1,
    degree_gamma = 8,
    lambda_ridge = 1e-6,
    verbose = TRUE
) {
  
  n <- dim(x)[1]
  m <- length(time)
  p <- dim(x)[3]
  
  ## init labels via kmeans
  Xstack <- matrix(0, n, m * p)
  for (i in 1:n) for (t in 1:m) {
    tmp <- x[i, t, ]; tmp[!is.finite(tmp)] <- 0
    Xstack[i, ((t - 1) * p + 1):(t * p)] <- tmp
  }
  cl0 <- kmeans(Xstack, centers = K, nstart = 10)$cluster
  cl_mat <- matrix(rep(cl0, m), n, m)
  
  ## init gamma
  gamma <- array(0, dim = c(K, m, p))
  for (t in 1:m) for (k in 1:K) {
    ids <- which(cl_mat[, t] == k)
    if (length(ids) > 0)
      gamma[k, t, ] <- colMeans(x[ids, t, ], na.rm = TRUE)
  }
  
  gamma_old <- gamma + 1
  
  for (iter in 1:m.iter) {
    
    gamma_old <- gamma
    
    ## adaptive-lasso merge detection
    sub_mat <- detect_merges_adalasso_level(
      gamma = gamma, time = time,
      lambda_level = lambda_level,
      tol_zero = tol_zero,
      df_spline = df_spline_merge,
      lambda_pilot = lambda_pilot,
      gamma_adapt = gamma_adapt
    )
    
    ## merge-aware spline fitting
    gamma <- fit_piecewise_splines(
      x = x, cl_mat = cl_mat,
      gamma_raw = gamma, sub_mat = sub_mat,
      time = time,
      degree = degree_gamma,
      lambda_ridge = lambda_ridge
    )
    
    ## relabel
    #for (t in 1:m) {
    #  Xt <- x[, t, ]
    #  dist_mat <- matrix(0, n, K)
    #  for (k in 1:K) {
    #    Gk <- matrix(gamma[k, t, ], n, p, byrow = TRUE)
    #    dist_mat[, k] <- rowMeans(abs(Xt - Gk), na.rm = TRUE)
    #  }
    #  cl_mat[, t] <- max.col(-dist_mat)
    #}
    
    # cl_vec: length-n vector (one label per individual)
    cl_vec <- integer(n)
    dist_mat <- matrix(0, nrow = n, ncol = K)
    
    for (k in 1:K) {
      
      # gamma_k over all time: m x p
      Gk_all <- gamma[k, , ]   # m x p
      
      for (i in 1:n) {
        # x_i over all time: m x p
        Xi <- x[i, , ]
        dist_mat[i, k] <- sum(abs(Xi - Gk_all), na.rm = TRUE)
      }
    }
    
    cl_vec <- max.col(-dist_mat)
    
    # replicate across time if downstream code expects n x m
    cl_mat <- matrix(rep(cl_vec, each = m), nrow = n, ncol = m, byrow = TRUE)
    
    diffv <- sqrt(mean((gamma - gamma_old)^2))
    #if (verbose) cat(sprintf("iter=%d | diff=%.6f\n", iter, diffv))
    if (diffv < eps) break
  }
  
  list(gamma = gamma, cl_mat = cl_mat, sub_mat = sub_mat)
}


#########################################
MDSP_nonpara = function(Y, X, Z = NULL, time=NULL, m, n, h_t = NULL, l1 = NULL, l2 = 1e-3,
                        l_vec = NULL, kp_vec = NULL, K=2, loc = 6, beta.int=NULL ,
                        degree_gamma = 10, gamma_adapt = 1){
  
  #############################################
  #############  ADMMM Algorithm ##############
  #############################################
  
  {
    p.x = dim(X)[3]
    
    if(is.null(Z)){
      
      p.z = 0
      
    }else{
      
      p.z = dim(Z)[3]
      
    }
    
    ### set tunning hyperparameters and initial values
    ## hyper para
    #kpp=kp # ADMM step
    # L1 regularization
    if(is.null(l_vec)){
      
      lm.s=rep(l1,m)
      
    }else{
      
      if(length(l_vec) != m){
        print("error with L1 regularization")
      }
      
      lm.s = l_vec
      
    }
    ## initials
    if (is.null(beta.int)){
      
      nu.int=array(0, dim = c(n, m, p.x))
      
      if(is.null(kp_vec)) print("beta.int and kp can not both be null")
      
    }else{
      
      nu.int=beta.int
      
      if(is.null(kp_vec)){
        
        kp_vec <- rep(NA,m)
        
        if(m > 1){
          
          for (i in 1:m) {
            
            
            beta_norm = mean(apply( beta.int[,i,,drop = FALSE],2, function(x) {
              sqrt(mean(x^2))}))
            
            rate = beta_norm/50
            
            kp_vec[i] = lm.s[i]/rate
            
          }
          
        }else{
          
          beta_norm = mean( sqrt(mean(beta.int^2)))
          
          rate = beta_norm/50
          
          kp_vec[1] = lm.s/rate
          
        }
      }
    }
    
    # once, outside the while-loop
    na_mask_beta <- is.na(beta.int)
    lm.int=array(0, dim = c(n, m, p.x))
    
    lm.h.new = lm.int
    nu.h.new = nu.int
    #a.h.new = rep(0,p.a)
    b.h.new = array(0, dim = c(n, m, p.x))#matrix(0,ncol= length(time), nrow = p.x)
    a.h.new = matrix(0, ncol= m, nrow = p.z)
    theta.h.new = array(0, dim = c(n, m, p.x + p.z))#matrix(0,ncol= length(time), nrow = p.x + p.z)
    theta.h = array(1, dim = c(n, m, p.x + p.z))#matrix(1,ncol= length(time), nrow = p.x + p.z)
    theta.h.new[ na_mask_beta] <- NA
    theta.h[ na_mask_beta] <- NA
    
    iter = 1
    m.iter = 500
    eps = 0.01
    error = rep(NA,m.iter-1)
    
  }
  
  
  ##### Main loop
  while(iter < m.iter & sqrt(mean((theta.h-theta.h.new)^2, na.rm= TRUE)) > eps){
    
    error[iter] <- sqrt(sum((theta.h - theta.h.new)^2/prod(dim(theta.h) ), na.rm = TRUE))
    print(error[iter])
    
    lm.h=lm.h.new
    nu.h=nu.h.new
    b.h=b.h.new
    a.h=a.h.new
    theta.h=theta.h.new
    
    
    
    ############## theta updates ################
    {
      
      
      if(p.z > 0){
        
        for (i in 1:m) {
          
          kpp <- kp_vec[i]
          temp_ind <- complete.cases(X[, i, ])
          X_m <- as.matrix(X[ temp_ind, i, ])
          Z_m <-  as.matrix(Z[ temp_ind, i, ])
          Y_m <-  as.matrix(Y[ temp_ind, i  ])
          n_m <- sum(temp_ind)
          
          X.td <- cbind(Z_m , t(KhatriRao(diag(n_m), t(X_m))))
          
          temp_theta = solve(t(X.td) %*% X.td + bdiag(list(matrix(0,p.z,p.z), kpp * diag(n_m*p.x)))) %*%
            ( t(X.td)%*%Y_m + c( rep(0,p.z),
                                 kpp* as.vector(t(nu.h[temp_ind,i,])) - as.vector(t(lm.h[temp_ind,i,])) ) )
          
          theta.h.new[temp_ind,i,1:p.z] <-  matrix(temp_theta[1:p.z,1], nrow = n_m, ncol = p.z, byrow = TRUE)
          theta.h.new[temp_ind,i,(p.z + 1) : (p.x + p.z)] <- t(matrix(temp_theta[(p.z + 1):nrow(temp_theta),], nrow = p.x))
          
          
        }
        
      }else{
        
        for (i in 1:m) {
          
          kpp <- kp_vec[i]
          temp_ind <- which(complete.cases(X[, i, ]))
          X_m <- as.matrix(X[ temp_ind, i, ])
          Y_m <-  as.matrix(Y[ temp_ind, i  ])
          n_m <- length((temp_ind))
          
          
          for (j in 1:n_m) {
            
            theta.h.new[temp_ind[j],i,] <- solve(   X_m[j,] %*%   t(X_m[j,]) + kpp * diag(p.x)) %*%
              ( X_m[j,] *  Y_m[j] + kpp * nu.h[temp_ind[j],i,] - lm.h[temp_ind[j],i,])
            
          }
          
        }
        
        
      }
    } #beta.true#
    
    
    
    b.h.new = theta.h.new[, ,(p.z+1):(p.z + p.x), drop = FALSE]
    
    
    ######### nu,gamma updates ###############
    
    B.h.new =  b.h.new
    Nu.h.new = nu.h
    Lm = lm.h
    g.h = array(0, dim=c(K,p.x,m))
    g.h.new = array(0, dim=c(K,p.x,m))
    
    
    sp.0 <- gamma_spline(b.h.new, time = time,
                         50, 0.001,  K =K, gamma_adapt = gamma_adapt)
    
    cl.h = sp.0$cl_mat
    g.h = sp.0$gamma - 0.05
    g.h.new = sp.0$gamma  # it cause the divergence
    
    iter.nu=1
    
    ### find gamma, nu
    
    while (iter.nu < 500 & sqrt(mean((g.h - g.h.new)^2, na.rm=TRUE)) > 0.01) {
      
      g.h <- g.h.new
      
      for (j in 1:m) {
        
        kpp <- kp_vec[j]
        
        b.t.new <- b.h.new[ , j, ] + kpp^{-1} * Lm[ , j, ]
        b.t.new <- as.matrix(b.t.new)
        
        for (i in 1:n) {
          
          temp_b <- as.matrix(b.t.new[i, ])
          
          temp_g <- as.matrix(g.h.new[cl.h[i, j], j, ])
          
          nu <- temp_g +
            sign(temp_b - temp_g) *
            as.numeric(abs(temp_b - temp_g) > lm.s[j] / kpp) *
            (abs(temp_b - temp_g) - lm.s[j] / kpp)
          
          Nu.h.new[i, j, ] <- nu
        }
      }
      
      ## Propagate NA pattern from beta to Nu.h.new
      #Nu.h.new[na_mask_beta] <- NA
      
      sp.0 <- gamma_spline(Nu.h.new, time = time,
                           50, 0.001,  K =K, gamma_adapt = gamma_adapt)
      
      g.h.new <- sp.0$gamma
      cl.h    <- sp.0$cl_mat
      
      iter.nu <- iter.nu + 1
      
      if (iter.nu == 500) {
        print("WARNING: maximum iteration (nu) achieved")
      }
      
    }
    
    
    group_alg <- cl.h
    
    ### lambda updates
    
    nu.h.new = Nu.h.new
    
    lm.h.new <- lm.h
    
    tau <- 0.9
    for (i in 1:m) {
      
      
      s <- sqrt(mean((nu.h[,i, ] - nu.h.new[,i, ])^2, na.rm= TRUE))
      r <- sqrt(mean((b.h.new[,i, ] - nu.h.new[,i, ])^2, na.rm= TRUE))
      
      if(r > s){
        
        kp_vec[i] <- tau * kp_vec[i]
        
      }else{
        
        if(r < s){
          
          kp_vec[i] <- (1/tau) * kp_vec[i]
          
        }
      }
      
      
      
      lm.h.new[, i, ] <- lm.h[, i, ] + kp_vec[i] * ( b.h.new[, i, ] - nu.h.new[, i, ])
      
    }
    
    
    
    iter = iter+1
    
    
    end_time <- Sys.time()
    
    if(iter==500) {print("WARNING: maxium iteration (nu) achieves")}
    
  }
  
  
  
  sp.h <- gamma_spline(Nu.h.new,time = time, 50, 0.001, K =K,
                       gamma_adapt = gamma_adapt)
  
  g.h.new = sp.h$gamma
  
  cl.h = sp.h$cl_mat
  sub_mat = sp.h$sub_mat
  
  return(list( beta.int = beta.int, alpha = a.h.new,  beta = b.h.new, Nu.h.new = Nu.h.new,
               gamma = g.h.new, group = cl.h, sub_mat = sub_mat, it.errors=error[1:(iter-1)]))
  
  
}



