stacking_design <- function(cost, d, n.init, epsilon, alpha, tt=2, MM=1, Lmax, f, k=2.5, xnew, n.max=1000, log.fg=TRUE, norm="L2",model.save=FALSE, save.n=FALSE, n.sobol=1000, ...){
  ## parameter setting 
  # cost: costs of simulations (a vector of length L); if NULL, it'll be provided when simulating the initial samples
  # d: dimension of X (scalar)
  # n.init: initial sample size (scalar)
  # epsilon: error tolerance (scalar)
  # alpha: converge rate of simulation (scalar)
  # tt: mesh size rate (scalar)
  # MM: largest mesh size (scalar)
  # Lmax: maximum number of the fidelity level (scalar)
  # f: simulation function (a function with arguments X and l)
  # xnew: test data (a data frame)
  # k: the smoothness parameter for the Matern kernel function. If NULL, it'll be estimated by LOOCV  
  # xnew: test data (a data frame)
  # n.max: the maximum number of sample size. If one of the sample size is too large, it will be very slow (scalar)
  # norm: which norm to be achieved when epsilon is given. Choose either "uniform" or "L2".
  # log.fg: print some logs? (logic)
  # model.save: whether the fitted models are saved (logic)
  # save.n: save n in each step? (logic)
  # n.sobol: number of low-discrepancy points to approximate the L2 or uniform norm 
  # ...: argument of GP.fit
  
  if(!is.null(k)) if(all(c(1,3,5,7,9)/2!=k)) stop("For Matern kernel, k must be 0.5, 1.5, 2.5, 3.5, 4.5, or NULL.")
  if(is.null(dim(xnew))) xnew <- matrix(xnew, ncol = 1)
  
  # initialization
  X <- D <- y <- z <- fit <- vector("list", Lmax)
  r <- rkhs.norm <- theta.norm <- rep(0, Lmax)
  P <- V <- B <- matrix(0, ncol=Lmax, nrow=nrow(xnew))
  if(is.null(cost)) cost <- rep(NA, Lmax)
  b <- rep(0,Lmax)
  alpha.set <- alpha
  if(save.n) n.save <- list() else n.save <- NULL
  var.save <- rep(0,Lmax)
  alpha.save <- rep(0,Lmax)
  
  # step (3): start with L=1
  L <- 0
  n <- n.old <- n.init # record previous sample size
  tol.est <- Inf
  while(tol.est > epsilon/2 & L<Lmax){ # if bias < e/2 then stop
    
    L <- L+1
    # step (4): initial estimation
    X[[L]] <- sobol(n.init, d, scrambling = 2) # domain: [-1,1]
    if(is.na(cost[L])) { # if no cost
      sim.out <- f(X[[L]], L, MM, tt, return.time=TRUE)
      y[[L]] <- sim.out$value
      cost[L] <- median(sim.out$time)
    }else{
      y[[L]] <- f(X[[L]], L, MM, tt)
    }
    
    if(L==1){
      z[[L]] <- y[[L]]
    }else{
      z[[L]] <- y[[L]] - y[[L-1]][1:n.init]
    }
    
    # if the sample size is too small, k will be given
    if(is.null(k) & n.init < 10*d) {
      fit[[L]] <- GP.fit(X[[L]], z[[L]], k=3.5, ...)
    }else{
      fit[[L]] <- GP.fit(X[[L]], z[[L]], k=k, ...)
    }
    
    k.min <- min(sapply(fit[1:L], function(x) x$k))
    for(l in 1:L){
      theta.norm[l] <- sqrt(sum((fit[[l]]$theta)^2))
      rkhs.norm[l] <- sqrt(t(fit[[l]]$y) %*% fit[[l]]$Ki %*% fit[[l]]$y)
      r[l] <- (theta.norm[l]^fit[[l]]$k * rkhs.norm[l] / cost[l])^(d/(k.min+d))
      
      if(l>1) if(r[l]>r[l-1]) r[l] <- r[l-1]
    }
    
    # step (5): determine sample size nl given error tolerance
    search.N <- function(N){
      n.out <- pmax(n.init, floor(N*r[1:L]))
      if(L>1) if(any(n.out[1:(L-1)] < n.old)) return(1000)
      out <- rep(0,L)
      for(l in 1:L){
        Xl <- sobol(n.out[l], d, scrambling = 2) 
        R <- sqrt(distance(t(t(Xl)*fit[[l]]$theta)))
        K <- matern.kernel(R, k=fit[[l]]$k)
        xnew <- sobol(n.sobol, d, scrambling=3)
        RX <- sqrt(distance(t(t(xnew)*fit[[l]]$theta), t(t(Xl)*fit[[l]]$theta)))
        KX <- matern.kernel(RX, k=fit[[l]]$k)
        
        power.X <- sqrt(pmax(0,1 - diag(KX %*% solve(K+diag(eps, n.out[l])) %*% t(KX))))
        if(norm == "L2"){
          out[l] <- median(power.X)*rkhs.norm[l]
        }else{
          out[l] <- max(power.X)*rkhs.norm[l]
        }
      }
      return(sum(out)-epsilon/2)
    }
    
    nn <- n.init - 1
    obj.out <- 1000
    while(nn <= (n.max-1) & obj.out >= 0){
      nn <- nn + 1
      obj.out <- search.N(nn/max(r[1:L]))
    }
    n <- pmin(n.max, floor(nn/max(r[1:L])*r[1:L]))
    n <- pmax(n.init, n)
    if(save.n) n.save[[L]] <- n
    var.save[L] <- obj.out + epsilon/2
    
    ### if previous n is larger then replace it
    if(L > 1) n[1:(L-1)][n.old > n[1:(L-1)]] <- n.old[n.old > n[1:(L-1)]]
    
    # steps (6) and (7): update design and model
    for(l in 1:L){
      if(l==L) {
        if(n[l] == n.init) next
      }else{
        if(n[l] == n.old[l]) next
      }
      
      X[[l]] <- sobol(n[l], d, scrambling = 2)
      if(d == 1) X[[l]] <- matrix(X[[l]], ncol=1)
      if(l==L) y[[l]] <- c(y[[l]], f(X[[l]][(n.init+1):n[l],,drop=FALSE], l, MM, tt))
      else y[[l]] <- c(y[[l]], f(X[[l]][(n.old[l]+1):n[l],,drop=FALSE], l, MM, tt))
      
      if(l==1){
        nl.fit <- n[l]
        z[[l]] <- y[[l]]
      }else{
        nl.fit <- min(n[l-1], n[l])
        z[[l]] <- y[[l]][1:nl.fit] - y[[l-1]][1:nl.fit]
      }
      
      if(is.null(k) & n[l] < 10*d) {
        fit[[l]] <- GP.fit(X[[l]], z[[l]], k=3.5, ...)
      }else{
        fit[[l]] <- GP.fit(X[[l]], z[[l]], k=k, ...)
      }
    }

    # step (8): see if bias is converged
    if(L < 3) {
      tol.est <- Inf
    }else{
      if(is.null(alpha.set)){
        # df_alpha <- data.frame(mesh_l=rep(MM/tt^(2:L), n[2:L]),
        #                        z=abs(unlist(z[2:L])))
        # lm.model <- lm(data=df_alpha,log(z)~log(mesh_l))
        # alpha <- lm.model$coefficients[2]
        alpha <- rep(0, L-2)
        for(l in 3:L) {
          alpha[l-2] <- mean(log(abs(z[[l-1]][1:n[l]]/z[[l]]))/log(tt))
        }
        alpha <- mean(alpha)
        alpha.save[L] <- alpha
      }
      if(norm == "L2"){
        tol.est <- sqrt(mean(pred.GP(fit[[L]], sobol(n.sobol,d, scrambling =3))$mu^2)) /(tt^(alpha)-1)
      }else{ 
        tol.est <- max(abs(pred.GP(fit[[L]], sobol(n.sobol,d, scrambling =3))$mu)) /(tt^(alpha)-1)
      }
    }
    
    n.old <- n
    b[L] <- tol.est
    if(log.fg) cat("L=",L," => n =",n, ", bias =", tol.est, ";\n")
  }
  
  # step (9): making predictions
  for(l in 1:L){
    out <- pred.GP(fit[[l]], xnew)
    P[,l] <- out$mu
    V[,l] <- out$sig2
    B[,l] <- abs(out$mu /(tt^(alpha)-1))
  }
  
  # output 
  if(!model.save){
    return(list(mean=rowSums(P), z=z, P=P, V=V, B=B, b=b, s2=rowSums(V), bias=B[,L], n=n, X=X, Y=y, L=L, cost=cost[1:L], n.save=n.save, alpha.save=alpha.save[1:L], var.save=var.save[1:L]))
  }else{
    return(list(mean=rowSums(P), z=z, P=P, V=V, B=B, b=b, s2=rowSums(V), bias=B[,L], n=n, X=X, Y=y, L=L, model=fit, cost=cost[1:L], n.save=n.save, alpha.save=alpha.save[1:L], var.save=var.save[1:L]))
  }
}

