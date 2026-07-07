library(splines)


homo_time <- function(X,Y,Z = NULL,time){

  p.x <- dim(X)[3]
  m <- dim(X)[2]
  beta_time <- matrix(NA, nrow = p.x, ncol = m)

  if(is.null(Z)){

    p.z <- 0

  }else{

    p.z <- dim(Z)[3]

  }

  for(i in 1:m){


    Y_m <- Y[,i , drop = FALSE]
    X_m <- cbind(Z[,i , ],X[,i ,])


    fit_lm <- lm(Y_m ~ X_m-1)

    beta_time[,i] <- t(fit_lm$coefficients[(p.z+1):(p.z+p.x)])
  }

  return(list(gamma = beta_time))
}



flexmix_fit <- function(X,Y,Z = NULL, m = m , K =2){


  p.x <- dim(X)[3]
  n <- dim(X)[1]

  if(is.null(Z)){

    p.z <- 0

  }else{

    p.z <- dim(Z)[3]

  }



  alpha_flex_time <- array(NA, dim = c(n, m,p.z))
  beta_flex_time <- array(NA, dim = c(n, m,p.x))
  gamma_flex_time <-  array(NA, dim = c(K, m,p.x))
  cl_flextime <- matrix(NA,nrow = n, ncol = m)

  ###########################################
  # Centers Estimation
  ###########################################

  for (i in 1:m) {


    y_m <- Y[,i]


    if(p.z > 0){

      X.td <- cbind(Z[,i,],  X[, i,])


    }else{

      X.td <- X[,i,]


    }

    temp_ind = complete.cases(X.td)

    flex_model <- flexmix( y_m[temp_ind] ~ X.td[temp_ind,]-1, k =K)
    prob_mat <- flex_model@posterior$scaled
    cl_cur <- rep(NA, n)

    if(ncol(prob_mat) > 1){

      b1<- parameters(flex_model, component = 1)
      b2<- parameters(flex_model, component = 2)

      cl_cur[temp_ind] <- as.factor(flex_model@cluster)


    }else{

      b1 <- b2 <- parameters(flex_model, component = 1)

      t_id <- sample(1:length(temp_ind),length(temp_ind)/2, replace = FALSE)
      cl_cur[temp_ind[ t_id ] ] <- 1
      cl_cur[temp_ind[-t_id] ]<- 2

    }


    ##############################
    # Resort the centers
    ##############################
    knorm <- rep(NA,K)
    for (j in 1:K) {

      var.name <- paste0("b",j)
      knorm[j] <- mean((get(var.name)[p.z+seq(1,p.x)]))

    }

    or <- order(knorm, decreasing = FALSE )
    center <- cbind(b1, b2)[,or, drop = FALSE]
    levels(cl_cur) <- or
    cl_cur <- as.numeric( as.vector( cl_cur ) )


    if(p.z > 0){

      alpha_flex_time[, i,] <- center[ 1:p.z,cl_cur]

    }

    #
    beta_flex_time[,i,] <- t(center[p.z + 1:p.x,cl_cur])
    #
    gamma_flex_time[,i,] <- t( center[p.z + 1:p.x, ])
    cl_flextime[,i] <-  cl_cur

  }


  return(list(gamma = gamma_flex_time, beta = beta_flex_time,
              alpha =  alpha_flex_time, group= cl_flextime))

}



Res_cluster <- function(Y,X,Z = NULL,m = m, k=2){


  p.x <- dim(X)[3]
  n <- dim(X)[1]

  if(is.null(Z)){

    p.z <- 0

  }else{

    p.z <- dim(Z)[3]

  }

  if(p.z > 0){

    X.td <- cbind(Z[,i,],  X[, i,])

  }else{

    X.td <- X[,i,]

  }


  fit_flex <- flexmix_fit(X = X, Y = Y, Z = Z, time = time)
  cl <- fit_flex$group[,1]
  g <- matrix(0, nrow = m, ncol = p.z + k*p.x)
  g_new <- g + 1
  b <- matrix(NA, nrow = p.x, ncol = length(time))
  a <- matrix(NA, nrow = p.z, ncol = length(time))
  #cl <- rep(NA,length(time))
  cl_new <- rep(NA,length(time))

  while(sqrt(mean((g - g_new)^2)) > 0.01) {


    g_new <- g

    # Update
    for(i in 1:m){

      ind <- which(time == t_ind[i])
      temp_diag <- matrix(0, nrow = length(ind), ncol = p.z + k * p.x)

      if(p.z > 0){

        temp_diag[,1:p.z] <- Z[ind,1:p.z]

      }

      xx <- X[ind,,drop=FALSE]
      yy = Y[ind,,drop=FALSE]
      for (j in 1:k) {

        ind_cl <- which(cl[ind] == j)
        temp_diag[ind_cl, p.z + seq(1,p.x) + (j-1) * p.x] <- xx[ind_cl,]

      }

      temp_fit <- lm(yy ~ temp_diag)
      g[i,] <- solve(t(temp_diag) %*% temp_diag) %*% t(temp_diag) %*% Y[ind,]

    }


    # Initial

    for (t in 1:length(time)) {

      #for (j in 1:n) {
      i <- which(t_ind == time[t])

      if(p.z > 0){
        r1 <- as.numeric( Y[t] - X.td[t,,drop=FALSE] %*% c( g[i, c(1:(p.z + p.x))]))
        r2 <- as.numeric( Y[t] - X.td[t,,drop=FALSE] %*% c( g[i, c(1:p.z, p.z+ p.x + seq(1:p.x))]))
      }else{

        r1 <- as.numeric( Y[t] - X.td[t,,drop=FALSE] %*% c( g[i, c(1:( p.x))]))
        r2 <- as.numeric( Y[t] - X.td[t,,drop=FALSE] %*% c( g[i, c( p.x + seq(1:p.x))]))
      }
      if(mean(r1^2) < mean(r2^2)){

        b[,t] <- g[i, p.z + seq(1:p.x)]
        cl[t] <- 1

      }else{

        b[,t] <- g[i, seq(1:p.x) + p.x + p.z]
        cl[t] <- 2

      }

      if(p.z > 0){

        a[,t] <- g[i,seq(1:p.z)]

      }


    }

  }


  g <- g[ ,c( seq(p.z+1, p.z + k*p.x) )]

  return(list(beta= b, alpha = a, group = cl, gamma = g))

}


DS_nonpara <- function(Y, X , Z = NULL, time ,m , n , l1, beta.int, gap = 3, type = 1, loc =6){


  fit_time <- MDSP_nonpara(Y = Y, X = X, Z = Z, time =time,
                           m = m, n =n,  l1 = l1, beta.int = beta.int)

  beta_s = fit_time$beta
  p.x = dim(X)[3]

  for (count in 1:2) {

    b <- fit_time$beta
    g <- fit_time$gamma
    temp_cl <- fit_time$group
    w <- 0.5

    for (i in 1:n) {

      id_ind <- which(complete.cases(as.matrix(b[i,,])))

      if(length(id_ind) > gap){


        cl_ad <- rep(NA, length(id_ind))

        for (j in 1:length(id_ind)) {

          if(j > gap){

            cl_ad[j]<- mean(temp_cl[i, (j-gap) : (j-1)])

          }else{

            cl_ad[j] <- temp_cl[i,j]

          }
        }

        cl_ad <- as.numeric( cl_ad >= 1.5 ) + 1
        temp_g <- matrix(NA, ncol = p.x, nrow = length(id_ind))

        for (tt in 1:length(id_ind)) {


          if(cl_ad[tt] == 1){

            temp_g[tt, ] <- g[1,id_ind[tt],]


          }else{

            temp_g[tt,] <- g[2,id_ind[tt],]


          }

        }


        bi <- w * temp_g + (1-w) * b[i,,]
        beta_s[i,,]  <- bi



      }else{

        beta_s[i,,] <- b[i,,]

      }

    }

    fit_time <- MDSP_nonpara(Y = Y, X = X, Z = Z, time =time,
                             m = m, n =n, h_t = h_t, l1 = l1, beta.int = beta_s)
  }

  for (i in 1:3) {

    beta_s <- fit_time$beta
    fit_time <- MDSP_nonpara(Y = Y, X = X, Z = Z, time =time,
                             m = m, n =n, h_t = h_t, l1 = l1, beta.int = beta_s)
  }

  return(fit_time)
}


DS_para <- function(Y, X , Z = NULL, h_t, time ,m , n , l1, beta.int, gap = 3, type = 1, loc =6){


  fit_time <- MDSP_para(
    Y = Y, X = X, Z=Z, m = m, n = n, h_t = h_t,
    l1 = l1, beta.int = beta.int
  )

  beta_s = fit_time$beta
  p.x = dim(X)[3]

  for (count in 1:2) {

    b <- fit_time$beta
    g <- fit_time$gamma
    temp_cl <- fit_time$group
    w <- 0.5

    for (i in 1:n) {

      id_ind <- which(complete.cases(as.matrix(b[i,,])))

      if(length(id_ind) > gap){


        cl_ad <- rep(NA, length(id_ind))

        for (j in 1:length(id_ind)) {

          if(j > gap){

            cl_ad[j]<- mean(temp_cl[i, (j-gap) : (j-1)])

          }else{

            cl_ad[j] <- temp_cl[i,j]

          }
        }

        cl_ad <- as.numeric( cl_ad >= 1.5 ) + 1
        temp_g <- matrix(NA, ncol = p.x, nrow = length(id_ind))

        for (tt in 1:length(id_ind)) {


          if(cl_ad[tt] == 1){

            temp_g[tt, ] <- g[1,id_ind[tt],]


          }else{

            temp_g[tt,] <- g[2,id_ind[tt],]


          }

        }


        bi <- w * temp_g + (1-w) * b[i,,]
        beta_s[i,,]  <- bi



      }else{

        beta_s[i,,] <- b[i,,]

      }

    }

    fit_time <- MDSP_para(
      Y = Y, X = X, Z=Z, m = m, n = n, h_t = h_t,
      l1 = l1, beta.int = beta_s
    )
  }

  for (i in 1:3) {

    beta_s <- fit_time$beta
    fit_time <- MDSP_para(
      Y = Y, X = X, Z=Z, m = m, n = n, h_t = h_t,
      l1 = l1, beta.int = beta_s
    )
  }

  return(fit_time)
}


individual_trajectory <- function(Y, X, Z = NULL, time, K = 2, n_knots = 5) {

  # --- Combine Z and X if needed ---
  if (is.null(Z)) {
    p.z <- 0
  } else {
    p.z <- dim(Z)[3]
    X   <- abind::abind(Z, X, along = 3)
  }

  p.x <- dim(X)[3]
  n   <- dim(X)[1]
  m   <- dim(X)[2]

  # --- Spline basis over time ---
  # Use n_knots as df so it actually matters
  s <- splines::ns(as.vector(time), intercept = TRUE, df = n_knots)
  d <- ncol(s)

  # --- Storage ---
  tra_b    <- array(NA, dim = c(n, m, p.x))
  tra_beta <- array(NA, dim = c(n, m, p.x - p.z))
  tra_g    <- array(NA, dim = c(K, m, p.x - p.z))
  tra_cl   <- matrix(NA, nrow = n, ncol = m)
  alpha_est <- matrix(NA, nrow = p.x * d, ncol = n)

  #################################################
  # 1. Individual Trajectory (subject–level alpha)
  #################################################
  for (i in 1:n) {

    t_i <- which(stats::complete.cases(X[i, , ]))
    m_i <- length(t_i)

    xi <- matrix(0, nrow = m_i, ncol = p.x * m_i)
    ai <- matrix(0, nrow = p.x * m_i, ncol = p.x * d)

    for (j in 1:m_i) {
      xind <- seq_len(p.x) + (j - 1L) * p.x
      xi[j, xind] <- X[i, t_i[j], ]

      for (jj in 1:p.x) {
        aind_x <- seq_len(d) + (jj - 1L) * d
        ai[xind[jj], aind_x] <- s[t_i[j], ]
      }
    }

    XB  <- xi %*% ai
    est <- solve(t(XB) %*% XB + diag(1e-4, ncol(XB))) %*% t(XB) %*% Y[i, t_i]
    alpha_est[, i] <- est
  }

  #################################################
  # 2. Temporal Centers (time-wise K-means)
  #################################################
  for (i in 1:m) {

    rind <- which(stats::complete.cases(X[, i, ]))
    base <- s[i, , drop = FALSE]

    # reconstruct beta_{it} from alpha_est
    for (j in 1:p.x) {
      indx <- seq_len(d) + (j - 1L) * d
      b    <- base %*% alpha_est[indx, rind, drop = FALSE]
      tra_b[rind, i, j] <- b
    }

    # K-means on heterogeneous covariates only
    c <- stats::kmeans(tra_b[rind, i, (p.z + 1):p.x],
                       centers = K)

    c_cl   <- as.factor(c$cluster)
    c_cent <- c$centers

    # reorder clusters by norm (optional, for identifiability)
    knorm <- apply(c_cent, 1, sum)
    or    <- order(knorm, decreasing = FALSE)
    c_cent <- c_cent[or, , drop = FALSE]
    levels(c_cl) <- or
    c_cl <- as.numeric(c_cl)

    tra_cl[rind, i]      <- c_cl
    tra_g[, i, ]         <- c_cent
    tra_beta[rind, i, ]  <- c_cent[c_cl, , drop = FALSE]
  }

  list(beta = tra_beta, gamma = tra_g, cl = tra_cl, a = alpha_est)
}



vec_trajectory_kmeans <- function(Y, X, Z = NULL, time, ID = NULL,
                                  m = NULL, K = 2, n_knots = 5) {

  # --- Combine Z and X if needed ---
  if (is.null(Z)) {
    p.z <- 0
  } else {
    p.z <- dim(Z)[3]
    X   <- abind::abind(Z, X, along = 3)
  }

  p.x <- dim(X)[3]
  n   <- dim(X)[1]
  m   <- dim(X)[2]   # override m from data

  # --- Spline basis ---
  s <- splines::ns(as.vector(time), intercept = TRUE, df = n_knots)
  d <- ncol(s)

  # --- Storage ---
  tra_b     <- array(NA, dim = c(n, m, p.x))
  tra_beta  <- array(NA, dim = c(n, m, p.x - p.z))
  tra_g     <- array(NA, dim = c(K, m, p.x - p.z))
  tra_cl    <- matrix(NA, nrow = n, ncol = m)
  alpha_est <- matrix(NA, nrow = p.x * d, ncol = n)

  #################################################
  # 1. Individual alpha per subject
  #################################################
  for (i in 1:n) {

    t_i <- which(stats::complete.cases(X[i, , ]))
    m_i <- length(t_i)

    xi <- matrix(0, nrow = m_i, ncol = p.x * m_i)
    ai <- matrix(0, nrow = p.x * m_i, ncol = p.x * d)

    for (j in 1:m_i) {
      xind <- seq_len(p.x) + (j - 1L) * p.x
      xi[j, xind] <- X[i, t_i[j], ]

      for (jj in 1:p.x) {
        aind_x <- seq_len(d) + (jj - 1L) * d
        ai[xind[jj], aind_x] <- s[t_i[j], ]
      }
    }

    XB  <- xi %*% ai
    est <- solve(t(XB) %*% XB + diag(1e-4, ncol(XB))) %*% t(XB) %*% Y[i, t_i]
    alpha_est[, i] <- est
  }

  #################################################
  # 2. Reconstruct beta_{it} for all (i,t)
  #################################################
  for (t in 1:m) {
    base <- s[t, , drop = FALSE]
    for (j in 1:p.x) {
      indx <- seq_len(d) + (j - 1L) * d
      b    <- base %*% alpha_est[indx, , drop = FALSE]
      tra_b[, t, j] <- b
    }
  }

  # vectorize full trajectories for heterogeneous covariates
  km_vec <- matrix(NA, nrow = n, ncol = m * (p.x - p.z))
  for (i in 1:n) {
    km_vec[i, ] <- as.vector(tra_b[i, , (p.z + 1):p.x, drop = FALSE])
  }

  #################################################
  # 3. K-means clustering on vectorized trajectories
  #################################################
  c <- stats::kmeans(km_vec, centers = K)

  c_cl   <- as.factor(c$cluster)
  c_cent <- c$centers  # K x [m * (p.x - p.z)]

  # reorder clusters by norm for label identifiability
  knorm <- apply(c_cent, 1, sum)
  or    <- order(knorm, decreasing = FALSE)
  c_cent <- c_cent[or, , drop = FALSE]
  c_cl_relab <- match(c_cl, or)

  # subject has same cluster across all time points
  for (t in 1:m) {
    tra_cl[, t] <- as.numeric(c_cl_relab)
  }

  # fill temporal centers tra_g and subject-specific beta
  for (t in 1:m) {
    idx_t <- seq_len(p.x - p.z) + (t - 1L) * (p.x - p.z)

    for (k_ind in 1:K) {
      tra_g[k_ind, t, ] <- c_cent[k_ind, idx_t]
    }

    tra_beta[, t, ] <- tra_g[as.numeric(c_cl_relab), t, ]
  }

  list(beta = tra_beta, gamma = tra_g, cl = tra_cl, a = alpha_est)
}

vec_trajectory_hierarchical <- function(Y, X, Z = NULL, time, ID = NULL,
                           m = NULL, K = 2, n_knots = 5) {

  # --- Combine Z and X if needed ---
  if (is.null(Z)) {
    p.z <- 0
  } else {
    p.z <- dim(Z)[3]
    X   <- abind::abind(Z, X, along = 3)
  }

  p.x <- dim(X)[3]
  n   <- dim(X)[1]
  m   <- dim(X)[2]   # override m from data

  # --- Spline basis ---
  s <- splines::ns(as.vector(time), intercept = TRUE, df = n_knots)
  d <- ncol(s)

  # --- Storage ---
  tra_b     <- array(NA, dim = c(n, m, p.x))
  tra_beta  <- array(NA, dim = c(n, m, p.x - p.z))
  tra_g     <- array(NA, dim = c(K, m, p.x - p.z))
  tra_cl    <- matrix(NA, nrow = n, ncol = m)
  alpha_est <- matrix(NA, nrow = p.x * d, ncol = n)

  #################################################
  # 1. Individual alpha per subject
  #################################################
  for (i in 1:n) {

    t_i <- which(stats::complete.cases(X[i, , ]))
    m_i <- length(t_i)

    xi <- matrix(0, nrow = m_i, ncol = p.x * m_i)
    ai <- matrix(0, nrow = p.x * m_i, ncol = p.x * d)

    for (j in 1:m_i) {
      xind <- seq_len(p.x) + (j - 1L) * p.x
      xi[j, xind] <- X[i, t_i[j], ]

      for (jj in 1:p.x) {
        aind_x <- seq_len(d) + (jj - 1L) * d
        ai[xind[jj], aind_x] <- s[t_i[j], ]
      }
    }

    XB  <- xi %*% ai
    est <- solve(t(XB) %*% XB + diag(1e-4, ncol(XB))) %*% t(XB) %*% Y[i, t_i]
    alpha_est[, i] <- est
  }

  #################################################
  # 2. Reconstruct beta_{it} for all (i,t)
  #################################################
  for (t in 1:m) {
    base <- s[t, , drop = FALSE]
    for (j in 1:p.x) {
      indx <- seq_len(d) + (j - 1L) * d
      b    <- base %*% alpha_est[indx, , drop = FALSE]
      tra_b[, t, j] <- b
    }
  }

  # vectorize full trajectories for heterogeneous covariates
  km_vec <- matrix(NA, nrow = n, ncol = m * (p.x - p.z))
  for (i in 1:n) {
    km_vec[i, ] <- as.vector(tra_b[i, , (p.z + 1):p.x, drop = FALSE])
  }

  #################################################
  # 3. Hierarchical clustering on vectorized trajectories
  #################################################
  D  <- dist(as.matrix(km_vec))
  hc <- hclust(D, method = "ward.D2")
  cl <- cutree(hc, k = K)   # length = n

  # Compute cluster centers in the vectorized space
  c_cent <- matrix(NA, nrow = K, ncol = ncol(km_vec))
  for (k_ind in 1:K) {
    if (sum(cl == k_ind) == 1) {
      c_cent[k_ind, ] <- km_vec[cl == k_ind, ]
    } else {
      c_cent[k_ind, ] <- colMeans(km_vec[cl == k_ind, , drop = FALSE])
    }
  }

  # reorder clusters by norm for label identifiability
  knorm <- apply(c_cent, 1, sum)
  or    <- order(knorm, decreasing = FALSE)
  c_cent <- c_cent[or, , drop = FALSE]
  cl_relab <- match(cl, or)

  # subject has same cluster across all time points
  for (t in 1:m) {
    tra_cl[, t] <- cl_relab
  }

  # fill temporal centers tra_g and subject-specific beta
  for (t in 1:m) {
    idx_t <- seq_len(p.x - p.z) + (t - 1L) * (p.x - p.z)

    for (k_ind in 1:K) {
      tra_g[k_ind, t, ] <- c_cent[k_ind, idx_t]
    }

    tra_beta[, t, ] <- tra_g[cl_relab, t, ]
  }

  list(beta = tra_beta, gamma = tra_g, cl = tra_cl, a = alpha_est)
}



alpha_trajectory_hierarchical <- function(Y, X, Z = NULL, time, ID = NULL,
                             m = NULL, K = 2, n_knots = 5) {

  # --- Combine Z and X if needed ---
  if (is.null(Z)) {
    p.z <- 0
  } else {
    p.z <- dim(Z)[3]
    X   <- abind::abind(Z, X, along = 3)
  }

  p.x <- dim(X)[3]
  n   <- dim(X)[1]
  m   <- dim(X)[2]

  # --- Spline basis over time ---
  s <- splines::ns(as.vector(time), intercept = TRUE, df = n_knots)
  d <- ncol(s)

  # --- Storage ---
  tra_b     <- array(NA, dim = c(n, m, p.x))
  tra_beta  <- array(NA, dim = c(n, m, p.x - p.z))
  tra_g     <- array(NA, dim = c(K, m, p.x - p.z))
  tra_cl    <- matrix(NA, nrow = n, ncol = m)
  alpha_est <- matrix(NA, nrow = p.x * d, ncol = n)

  #################################################
  # 1. Individual alpha per subject
  #################################################
  for (i in 1:n) {
    t_i <- which(stats::complete.cases(X[i, , ]))
    m_i <- length(t_i)

    xi <- matrix(0, nrow = m_i, ncol = p.x * m_i)
    ai <- matrix(0, nrow = p.x * m_i, ncol = p.x * d)

    for (j in 1:m_i) {
      xind <- seq_len(p.x) + (j - 1L) * p.x
      xi[j, xind] <- X[i, t_i[j], ]

      for (jj in 1:p.x) {
        aind_x <- seq_len(d) + (jj - 1L) * d
        ai[xind[jj], aind_x] <- s[t_i[j], ]
      }
    }

    XB  <- xi %*% ai
    est <- solve(t(XB) %*% XB + diag(1e-4, ncol(XB))) %*% t(XB) %*% Y[i, t_i]
    alpha_est[, i] <- est
  }

  #################################################
  # 2. Time-wise hierarchical clustering on beta
  #################################################
  for (i in 1:m) {

    rind <- which(stats::complete.cases(X[, i, ]))
    base <- s[i, , drop = FALSE]

    # reconstruct beta_{it} from alpha_est
    for (j in 1:p.x) {
      indx <- seq_len(d) + (j - 1L) * d
      b    <- base %*% alpha_est[indx, rind, drop = FALSE]
      tra_b[rind, i, j] <- b
    }

    # matrix of heterogeneous coefficients at time i
    mat_i <- tra_b[rind, i, (p.z + 1):p.x]

    # ---- hierarchical clustering here ----
    D  <- dist(as.matrix(mat_i))                # distance between subjects
    hc <- hclust(D, method = "ward.D2")
    cl <- cutree(hc, k = K)                    # cluster labels (length = #rind)

    # cluster centers as mean of rows in each cluster
    c_cent <- matrix(NA, nrow = K, ncol = ncol(mat_i))
    for (k_ind in 1:K) {
      if (sum(cl == k_ind) == 1) {
        c_cent[k_ind, ] <- mat_i[cl == k_ind, ]
      } else {
        c_cent[k_ind, ] <- colMeans(mat_i[cl == k_ind, ])
      }
    }

    # optional: reorder clusters by norm to stabilize labels
    knorm <- apply(c_cent, 1, sum)
    or    <- order(knorm, decreasing = FALSE)
    c_cent <- c_cent[or, , drop = FALSE]

    # relabel cl according to 'or'
    cl_relab <- match(cl, or)

    tra_cl[rind, i]     <- cl_relab
    tra_g[, i, ]        <- c_cent
    tra_beta[rind, i, ] <- c_cent[cl_relab, , drop = FALSE]
  }

  list(beta = tra_beta, gamma = tra_g, cl = tra_cl, a = alpha_est)
}

alpha_trajectory_kmeans <- function(Y, X, Z = NULL, time, ID = NULL,
                                    m = NULL, K = 2, n_knots = 5) {

  # --- Combine Z and X if needed ---
  if (is.null(Z)) {
    p.z <- 0
  } else {
    p.z <- dim(Z)[3]
    X   <- abind::abind(Z, X, along = 3)
  }

  p.x <- dim(X)[3]
  n   <- dim(X)[1]
  m   <- dim(X)[2]

  # --- Spline basis over time ---
  s <- splines::ns(as.vector(time), intercept = TRUE, df = n_knots)
  d <- ncol(s)

  # --- Storage ---
  tra_b     <- array(NA, dim = c(n, m, p.x))
  tra_beta  <- array(NA, dim = c(n, m, p.x - p.z))
  tra_g     <- array(NA, dim = c(K, m, p.x - p.z))
  tra_cl    <- matrix(NA, nrow = n, ncol = m)
  alpha_est <- matrix(NA, nrow = p.x * d, ncol = n)

  #################################################
  # 1. Individual alpha per subject
  #################################################
  for (i in 1:n) {
    t_i <- which(stats::complete.cases(X[i, , ]))
    m_i <- length(t_i)

    xi <- matrix(0, nrow = m_i, ncol = p.x * m_i)
    ai <- matrix(0, nrow = p.x * m_i, ncol = p.x * d)

    for (j in 1:m_i) {
      xind <- seq_len(p.x) + (j - 1L) * p.x
      xi[j, xind] <- X[i, t_i[j], ]

      for (jj in 1:p.x) {
        aind_x <- seq_len(d) + (jj - 1L) * d
        ai[xind[jj], aind_x] <- s[t_i[j], ]
      }
    }

    XB  <- xi %*% ai
    est <- solve(t(XB) %*% XB + diag(1e-4, ncol(XB))) %*% t(XB) %*% Y[i, t_i]
    alpha_est[, i] <- est
  }

  #################################################
  # 2. Time-wise K-means clustering on beta
  #################################################
  for (i in 1:m) {

    rind <- which(stats::complete.cases(X[, i, ]))
    base <- s[i, , drop = FALSE]

    # reconstruct beta_{it} from alpha_est
    for (j in 1:p.x) {
      indx <- seq_len(d) + (j - 1L) * d
      b    <- base %*% alpha_est[indx, rind, drop = FALSE]
      tra_b[rind, i, j] <- b
    }

    # heterogeneous coefficients at time i
    mat_i <- tra_b[rind, i, (p.z + 1):p.x]

    # ---- K-means here ----
    c <- stats::kmeans(mat_i, centers = K)

    c_cl   <- as.factor(c$cluster)
    c_cent <- c$centers  # K x (p.x - p.z)

    # reorder clusters by norm
    knorm <- apply(c_cent, 1, sum)
    or    <- order(knorm, decreasing = FALSE)
    c_cent <- c_cent[or, , drop = FALSE]

    # relabel cl according to 'or'
    cl_relab <- match(c_cl, or)

    tra_cl[rind, i]     <- as.numeric(cl_relab)
    tra_g[, i, ]        <- c_cent
    tra_beta[rind, i, ] <- c_cent[as.numeric(cl_relab), , drop = FALSE]
  }

  list(beta = tra_beta, gamma = tra_g, cl = tra_cl, a = alpha_est)
}
