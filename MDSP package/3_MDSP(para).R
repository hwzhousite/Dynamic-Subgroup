############## packages
require(mixtools)
require(Matrix)
library(splines)
library(locfit)
library(clue)
#load MDSP_gamma function
# install.packages("clue")  # if not installed


gamma_para <- function (x, h_t, m,  m.iter = 50, eps,  K = 2){


  #################################
  # p: p_x
  # n: number of identity
  # g: matrix for gamma
  # sp_mat: matrix stores regression and classification information
  #################################

  {

    p <- dim(x)[3]
    d <- dim(h_t)[2]
    g <- array(0, dim = c(K,m,p))
    G <- array(0, dim = c(K,p,d+1))
    cl_mat = matrix(NA, nrow = n, ncol =m)

    cluster_mat = matrix(0, nrow=n, ncol = m * p )

    # Initialization of cluster labels
    for (i in 1:n) {

      for (j in 1:m) {

        temp_x = x[i, j,]
        temp_x[is.na(temp_x)] = 0
        cluster_mat[i, seq(1,p) + (j-1)*p] <- x[i, j,]
      }
    }

    cl_temp <- kmeans( cluster_mat, center= K)
    cl <- as.factor( cl_temp$cluster )
    cl_mat[ ,1:m] <-  cl

    g = array(0, dim = c(K, m,p))
    for (i in 1:m) {

      for (j in 1:K) {

        g[j, i,] = apply( x[cl_mat[,i] == j,i, , drop=FALSE], 2 , function(x) mean(x, na.rm = TRUE) )

      }

    }

  }


  iter = 0
  g.a = g; g.b = g + 1

  start_time <- Sys.time()
  # Iteration for Convergence
  while (sqrt(mean((g.a - g.b)^2)) > eps & iter < m.iter) {

    g.b <- g.a
    iter <- iter + 1

    # Update the class by global info
    for (i in 1:K) {

      for (j in 1:p) {

        temp_fit <- lm(g.a[i, , j] ~ h_t)
        G[i, j,] <- temp_fit$coefficients
        g.a[i, ,j] <- temp_fit$fitted.values

      }

    }


    # Update the class by minimizing the Euclidean distance
    for (i in 1:m) {

      temp_g <- as.matrix(g.a[,i, ])
      temp_x <- as.matrix(x[,i, ])
      dis <- matrix(NA, ncol = K , nrow = n)


      for (j in 1:K) {

        cent_mat <- matrix(temp_g[j,], ncol = p, nrow = n, byrow = TRUE)

        dis[,j] <- apply(abs(temp_x - cent_mat),1,mean)

      }


      cl <- apply(dis, 1, which.min)
      cl_mat[, i] <- as.numeric( as.vector(cl))

    }

    cl = apply(cl_mat, 1, function(x) names(which.max(table(x))))
    cl_mat[,1:m] = as.numeric(cl)


  }


  end_time <- Sys.time()
  #print(end_time - start_time )

  g <- g.a



  return(list( gamma=g, G = G, cl_mat = cl_mat))
}
#########################################


MDSP_para = function(Y, X,  m, n, time = NULL, Z = NULL, h_t = NULL, l1 = NULL,
                     l_vec = NULL, kp_vec = NULL, K=2, cl.t = NULL, beta.int=NULL, tru_heter = NULL ){

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

      nu.int = beta.int

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


    lm.int=array(0, dim = c(n, m, p.x))

    lm.h.new = lm.int
    nu.h.new = nu.int
    #a.h.new = rep(0,p.a)
    b.h.new = array(0, dim = c(n, m, p.x))

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

          theta.h.new[,i,1:p.z] <-  matrix(temp_theta[1:p.z,1], nrow = n_m, ncol = p.z, byrow = TRUE)
          theta.h.new[,i,(p.z + 1) : (p.x + p.z)] <- t(matrix(temp_theta[(p.z + 1):nrow(temp_theta),], nrow = p.x))


        }

      }else{

        for (i in 1:m) {
          kpp <- kp_vec[i]
          X_m <- X[,i,]
          Z_m <- Z[,i,]
          Y_m <- Y[,i]
          n_m <- nrow(X_m)


          for (j in 1:n) {

            theta.h.new[j,i,] <- solve(   X_m[j,] %*%   t(X_m[j,]) + kpp * diag(p.x)) %*%
              ( X_m[j,] *  Y_m[j] + kpp * nu.h[j,i,] - lm.h[j,i,])

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

    sp.0 <- gamma_para(b.h.new, m =m , h_t = h_t,
                       m.iter =50, eps = 0.001, K =K)
    
    #hist(b.h.new[,1,],breaks = 50)

    cl.h = sp.0$cl_mat
    g.h = sp.0$gamma - 0.05
    g.h.new = sp.0$gamma # it cause the divergence

    iter.nu=1

    ### find gamma, nu

    while (iter.nu <500 & sqrt(mean((g.h-g.h.new)^2)) >0.01){

      g.h =  g.h.new

      for (j in 1:m) {

        kpp <- kp_vec[j]

        b.t.new <- b.h.new[ , j, ] + kpp^{-1} * Lm[ ,j, ]
        b.t.new <- as.matrix( b.t.new )

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



      sp.0 <- gamma_para(Nu.h.new, m =m , h_t = h_t,
                         m.iter =50, eps = 0.001, K =K)
      g.h.new = sp.0$gamma
      cl.h = sp.0$cl_mat

      iter.nu <- iter.nu + 1

      if(iter.nu==500) {print("WARNING: maxium iteration (nu) achieves")}

    }


    group_alg <- cl.h

    ### lambda updates

    nu.h.new = Nu.h.new

    lm.h.new <- lm.h

    tau <- 1.1
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



  sp.h <- gamma_para(nu.h.new, m =m , h_t = h_t,
                     m.iter =50, eps = 0.001, K =K)

  g.h.new = sp.h$gamma

  cl.h = sp.h$cl_mat

  G.model = sp.h$G

  return(list( beta.int = beta.int, alpha = a.h.new,  beta = b.h.new, nu = nu.h.new,
               gamma = g.h.new, group = cl.h, model =G.model, it.errors=error[1:(iter-1)]))


}


