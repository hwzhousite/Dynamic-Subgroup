############## packages ###########
require(mixtools)
require(Matrix)
library(splines)
library(locfit)

#load MDSP_gamma function
gamma_spline <- function (x, time,  m.iter = 500, eps, sp = "ns", K = 2){


  #################################
  # p: p_x
  # n: number of identity
  # g: matrix for gamma
  # sp_mat: matrix stores regression and classification information
  #################################
  #while(1)

  {
    m <- length(time)
    p <- dim(x)[3]
    g <- array(0, dim = c(K,m,p))
    g_sort <- g
    #sp_mat <- data.frame(t = rep(seq(1,m),each = n), cl = rep(0, n*m))

    sp_time <- matrix(NA, nrow =length(time), ncol = 1)
    model <- list()
    cl_mat = matrix(NA, nrow = n, ncol =m)

    # Initialization of cluster labels
    for (i in 1:m) {


      temp_id <- which(complete.cases(as.matrix(x[,i,]) ))
      temp_x <- as.data.frame(x[  temp_id,i ,])

      cl_temp <- kmeans(temp_x, center= K)

      cl <- as.factor( cl_temp$cluster )
      cent <- cl_temp$centers

      knorm <- rep(NA,K)
      for (j in 1:K) {

        knorm[j] <- sum(cent[j,])

      }

      or <- order(knorm, decreasing = FALSE )
      cent <- cent[or, ,drop = FALSE]
      levels(cl) <- or

      cl_mat[ ,i] <- as.numeric( as.vector(cl))
      g[ ,i, ]  <- cent

    }


    iter = 0
    g.a = g; g.b = g + 1

    start_time <- Sys.time()
    # Iteration for Convergence
    while (sqrt(mean((g.a - g.b)^2)) > eps & iter < m.iter) {

      g.b <- g.a
      iter <- iter + 1

      # Update the center by spline regression
      if(m > 1){
          for (i in 1:K){

            temp_model <- list()

            for (j in 1:p) {

              data <- data.frame(y = g.a[ i , , j  ], t = time)

              sp_fit <- lm(y ~  ns(t,df = 7), data = data)

              g.a[i, ,j ] <- sp_fit$fitted.values

              temp_model[[j]] <- sp_fit

            }

            model[[i]] <- temp_model

          }
      }

      # Update the class by minimizing the Euclidean distance
      for (i in 1:m) {


        temp_id <- which(complete.cases(as.matrix( x[,i,] ) ))
        temp_x <- as.matrix((x[  temp_id, i,]))
        temp_l <- length(temp_id)

        dis <- matrix(NA, ncol = K, nrow =temp_l)
        temp_g <- as.matrix( g.a[ , i, ])

        for (j in 1:K) {

          cent_mat <- matrix(temp_g[j, ],nrow = length(temp_id),  ncol = p)

          dis[,j] <- apply(abs(temp_x - cent_mat),1,mean)

        }


        cl <- apply(dis, 1, which.min)
        cl_mat[,i] <- as.numeric( as.vector(cl))


      }

    }
  }

  end_time <- Sys.time()
  #print(end_time - start_time )

  g <- g.a



  return(list(cl_mat = cl_mat, gamma=g))
}



#########################################
MDSP_nonpara = function(Y, X, Z = NULL, time=NULL, m, n, h_t = NULL, l1 = NULL, l2 = 1e-3,
                        l_vec = NULL, kp_vec = NULL, K=2, loc = 6, beta.int=NULL ){

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
              sqrt(mean(x^2,na.rm = TRUE))}), na.rm = TRUE)

            rate = beta_norm

            kp_vec[i] = lm.s[i]/rate

          }

        }else{

          beta_norm = mean( sqrt(mean(beta.int^2)))

          rate = beta_norm

          kp_vec[1] = lm.s/rate


        }

      }

    }


    lm.int=array(0, dim = c(n, m, p.x))

    lm.h.new = lm.int
    nu.h.new = nu.int
    #a.h.new = rep(0,p.a)
    b.h.new = array(0, dim = c(n, m, p.x))#matrix(0,ncol= length(time), nrow = p.x)
    a.h.new = matrix(0, ncol= m, nrow = p.z)
    theta.h.new = array(0, dim = c(n, m, p.x + p.z))#matrix(0,ncol= length(time), nrow = p.x + p.z)
    theta.h = array(1, dim = c(n, m, p.x + p.z))#matrix(1,ncol= length(time), nrow = p.x + p.z)

    iter = 1
    m.iter = 500
    eps = 0.01
    error = rep(NA,m.iter-1)

  }


  ##### Main loop
  while(iter < m.iter & sqrt(mean((theta.h-theta.h.new)^2)) > eps){

    error[iter] <- sqrt(sum((theta.h - theta.h.new)^2/prod(dim(theta.h) )))
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
          n_m <- length(temp_ind)

          X.td <- cbind(Z_m , t(KhatriRao(diag(n_m), t(X_m))))

          temp_theta = solve(t(X.td) %*% X.td + bdiag(list(matrix(0,p.z,p.z), kpp * diag(n_m*p.x)))) %*%
            ( t(X.td)%*%Y_m + c( rep(0,p.z),
                                 kpp* as.vector(t(nu.h[temp_ind,i,])) - as.vector(t(lm.h[temp_ind,i,])) ) )

          theta.h.new[,i,1:p.z] <- temp_theta[1:p.z,1]
          theta.h.new[,i,(p.z + 1) : (p.x + p.z)] <- t(matrix(temp_theta[(p.z + 1):nrow(temp_theta),], nrow = p.x))


        }

      }else{

        for (i in 1:m) {
          kpp <- kp_vec[i]
          temp_ind <- which(complete.cases(X[,i, ]))
          X_m <- X[temp_ind  ,i,]
          Y_m <- Y[ temp_ind ,i]

          for (j in temp_ind) {

            theta.h.new[j,i,] <- solve(   X_m[j,] %*%   t(X_m[j,]) + kpp * diag(p.x)) %*%
              matrix(( X_m[j, ] *  Y_m[j] + kpp* nu.h[j,i,] - lm.h[j,i,]), ncol = 1)
          }

        }
      }
    } #beta.true#


    for (i in 1:n) {

      temp_x <- as.matrix(X[i, , ])
      temp_ind <-which(complete.cases(temp_x))

      if(length(temp_ind) > 6){

        for (g in (p.z + 1):(p.z + p.x)) {

          temp_data <- data.frame(y = theta.h.new[i, temp_ind,g], t = time[ temp_ind] )

          band <- min(0.8, loc/nrow(temp_data))
          temp_fit <-  locfit(y~lp(t,nn=band, deg = 2),data=temp_data)

          #temp_fit <- loess(y ~ t, data = temp_data, degree = 2, span = 0.5)
          #temp_fit <- locpoly(temp_data$t, temp_data$y, bandwidth=0.15, degree=2)
          #table(temp_data$t %in% temp_fit$x)
          theta.h.new[i, temp_ind,g] <- predict(temp_fit,  newdata = data.frame(t = temp_data$t))

        }
      }

    }

    b.h.new = theta.h.new[, ,(p.z+1):(p.z + p.x), drop = FALSE]


    ######### nu,gamma updates ###############

    B.h.new =  b.h.new
    Nu.h.new = nu.h
    Lm = lm.h
    g.h = array(0, dim=c(K,p.x,m))
    g.h.new = array(0, dim=c(K,p.x,m))


    sp.0 <- gamma_spline(b.h.new, time = time,
                             50, 0.001,  K =K)

    cl.h = sp.0$cl_mat
    g.h = sp.0$gamma - 0.05
    g.h.new = sp.0$gamma  # it cause the divergence

    iter.nu=1

    ### find gamma, nu

    while (iter.nu <500 & sqrt(mean((g.h-g.h.new)^2)) >0.01){

      g.h =  g.h.new

      for (j in 1:m) {

        kpp <- kp_vec[j]

        b.t.new <- as.matrix(b.h.new[ , j, ] + kpp^{-1} * Lm[ ,j, ])


        for (i in 1:n) {

          temp_b <- as.matrix(b.t.new[i, ])

          temp_g <-  as.matrix(g.h.new[cl.h[i,j] , j, ])

          # What is lm.s
          nu = temp_g +
            sign(temp_b - temp_g) *
            as.numeric(abs(temp_b - temp_g)>lm.s[j]/kpp) *
            (abs(temp_b - temp_g)-lm.s[j]/kpp)

          Nu.h.new[i,j,] <- nu

        }

      }

      sp.0 <- gamma_spline( Nu.h.new, time = time,  50, 0.001,K =K)
      g.h.new = sp.0$gamma
      cl.h = sp.0$cl_mat

      iter.nu <- iter.nu + 1

      if(iter.nu==500) {print("WARNING: maxium iteration (nu) achieves")}

    }


    group_alg <- cl.h

    ### lambda updates

    nu.h.new = Nu.h.new

    lm.h.new <- lm.h

    tau <- 0.9
    for (i in 1:m) {


      s <- sqrt(mean((nu.h[,i, ] - nu.h.new[,i, ])^2))
      r <- sqrt(mean((b.h.new[,i, ] - nu.h.new[,i, ])^2))

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



  sp.h <- gamma_spline(Nu.h.new,time = time, 50, 0.001, K =K)

  g.h.new = sp.h$gamma

  cl.h = sp.h$cl_mat


  return(list( beta.int = beta.int, alpha = a.h.new,  beta = b.h.new, Nu.h.new = Nu.h.new,
               gamma = g.h.new, group = cl.h, it.errors=error[1:(iter-1)]))


}


#' Random-split initialization for MDSP
#'
#' Runs `max_iter` random binary splits of the samples, keeps the split whose
#' two per-time OLS coefficient matrices are farthest apart, then derives the
#' residual-based clustering and the treatment-coefficient initialization array.
#'
#' @param Y          n x m response matrix
#' @param Z          n x m x q_ctrl control-covariate array
#' @param X          n x m x p_treat treatment-covariate array
#' @param ZX         n x m x (q_ctrl + p_treat) combined design array
#' @param K          number of clusters (assumed 2 below, as in the original code)
#' @param n,m        number of samples / time points
#' @param p_treat    number of treatment covariates
#' @param q_ctrl     number of control covariates
#' @param true_groups,gamma_true  ground truth used only for evaluation
#' @param max_iter   number of random splits to try
#'
#' @return list with gamma_mat, cl_mat1, accuracy_split,
#'         gamma_hat_split_gamma, gamma_rmse_split_gamma, and int
init_random_split <- function(Y, Z, X, ZX,
                              K, n, m, p_treat, q_ctrl,
                              true_groups, gamma_true,
                              max_iter = 100) {
  
  dist      <- 0
  gamma_mat <- NULL
  
  for (iter in 1:max_iter) {
    temp_cl    <- rbinom(n = n, size = 1, prob = 0.5) + 1
    temp_gamma <- array(0, dim = c(K, m, p_treat + q_ctrl))
    for (tt in 1:m) {
      X_t <- cbind(Z[, tt, ], X[, tt, ])
      for (kk in 1:K) {
        temp_x   <- X_t[temp_cl == kk, , drop = FALSE]
        temp_y   <- Y[temp_cl == kk, tt, drop = FALSE]
        temp_vec <- solve(t(temp_x) %*% temp_x) %*% (t(temp_x) %*% temp_y)
        temp_gamma[kk, tt, ] <- temp_vec
      }
    }
    temp_dist <- sum((temp_gamma[1, , ] - temp_gamma[2, , ])^2)
    if (temp_dist > dist) {
      dist      <- temp_dist
      gamma_mat <- temp_gamma  # K x m x (q_ctrl + p_treat)
    }
  }
  
  ## residuals under each cluster's coefficients
  res_mat <- array(NA_real_, dim = c(n, m, K))
  for (i in 1:n) {
    for (tt in 1:m) {
      coef_vec1 <- c(gamma_mat[1, tt, ])
      coef_vec2 <- c(gamma_mat[2, tt, ])
      res_mat[i, tt, 1] <- abs(Y[i, tt] - sum(ZX[i, tt, ] * coef_vec1))
      res_mat[i, tt, 2] <- abs(Y[i, tt] - sum(ZX[i, tt, ] * coef_vec2))
    }
  }
  
  ## per-sample cluster assignment (constant across time)
  cl_mat1 <- matrix(NA_integer_, nrow = n, ncol = m)
  for (i in 1:n) {
    cl_mat1[i, ] <- which.min(c(
      sum(abs(res_mat[i, , 1])),
      sum(abs(res_mat[i, , 2]))
    ))
  }
  
  accuracy_split <- evaluate_group_accuracy(
    cl_mat1, true_groups
  )$overall_accuracy
  
  ## gamma for "Random Split Gamma"
  gamma_hat_split_gamma  <- gamma_mat[, , (q_ctrl + 1):(q_ctrl + p_treat), drop = FALSE]
  gamma_rmse_split_gamma <- gamma_rmse(gamma_hat_split_gamma, gamma_true)
  
  ## MDSP initialization array from the random-split clustering
  int <- array(NA_real_, dim = c(n, m, p_treat))
  for (tt in 1:m) {
    for (j in 1:p_treat) {
      int[, tt, j] <- gamma_mat[cl_mat1[, tt], tt, q_ctrl + j]
    }
  }
  
  list(
    gamma_mat              = gamma_mat,
    cl_mat1                = cl_mat1,
    accuracy_split         = accuracy_split,
    gamma_hat_split_gamma  = gamma_hat_split_gamma,
    gamma_rmse_split_gamma = gamma_rmse_split_gamma,
    int                    = int
  )
}