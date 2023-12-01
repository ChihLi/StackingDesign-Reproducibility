library(plgp)
eps <- 1e-8

# fit a Gaussian process model
GP.fit <- function(X, y, k, g=eps, 
                   init=rep(10, ncol(X)),
                   lower=0.1, upper=1000){
  # X is the input
  # y is output
  # k is the smoothness parameter in matern kernel
  # g is the nugget effect
  # init is the initial value of lengthscale parameter when performing optimization
  # lower is the lower bound of lengthscale parameter when performing optimization
  # upper is the upper bound of lengthscale parameter when performing optimization
  if(is.null(dim(X))) X <- matrix(X, ncol = 1)
  
  y <- scale(y, scale=FALSE)
  
  n <- length(y)
  if(!is.null(k)){ # if the smoothness is given
    obj <- function(par, X, Y) 
    {
      theta <- par
      R <- sqrt(distance(t(t(X)*theta)))
      K <- matern.kernel(R, k=k)
      Ki <- solve(K+diag(g,n))

      loocv <- mean((diag(1/diag(Ki)) %*% Ki %*% Y)^2)
      return(loocv)
    }
    
    tic <- proc.time()[3]
    
    outg <- optim(init, obj, 
                  method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y)
    toc <- proc.time()[3]
    
    theta <- outg$par
  }else{ # if the smoothness is not given, then it needs to be estimated
    k.candidate <- rev(c(3,5,7)/2) 
    theta.save <- matrix(0,nrow=length(k.candidate),ncol=ncol(X))
    loocv.save <- rep(0,length(k.candidate))
    
    for(i in 1:length(k.candidate)){
      obj <- function(par, X, Y) 
      {
        theta <- par
        R <- sqrt(distance(t(t(X)*theta)))
        K <- matern.kernel(R, k=k.candidate[i])
        Ki <- solve(K+diag(g,n))
        loocv <- mean((diag(1/diag(Ki)) %*% Ki %*% Y)^2)
        return(loocv)
      }
      
      tic <- proc.time()[3]
      
      outg <- optim(init, obj,
                    method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y)
      toc <- proc.time()[3]
      
      theta.save[i,] <- outg$par
      
      loocv.save[i] <- outg$value
    }
    k <- k.candidate[which.min(loocv.save)]
    theta <- theta.save[which.min(loocv.save),]
  }
  
  R <- sqrt(distance(t(t(X)*theta)))
  K <- matern.kernel(R, k=k)
  Ki <- solve(K+diag(g,n))
  
  return(list(theta = theta, k=k, Ki=Ki,
              g = g, X = X, y = y))
}

# predict on xnew based on a fitted GP model

pred.GP <- function(fit, xnew){
  
  xnew <- as.matrix(xnew)
  
  Ki <- fit$Ki
  theta <- fit$theta
  k <- fit$k
  g <- fit$g
  X <- fit$X
  y <- fit$y
  
  tau2hat <- drop(t(y) %*% Ki %*% y / nrow(X))
  
  RXX <- sqrt(distance(t(t(xnew)*theta)))
  RX <- sqrt(distance(t(t(xnew)*theta), t(t(X)*theta)))
  KXX <- matern.kernel(RXX, k=k)
  KX <- matern.kernel(RX, k=k)
  
  mup2 <- KX %*% Ki %*% y +  attr(y, "scaled:center") 
  Sigmap2 <- pmax(0,diag(tau2hat*(KXX + diag(g,nrow(xnew)) - KX %*% Ki %*% t(KX))))
  
  return(list(mu=mup2, sig2=Sigmap2))
}
