generate_data <- function(pb=1, ntrain=100, ntest=30, p=120, draw=FALSE, 
                          tslice=c(0,.2,.5, 1), relevant.slices=1) {
  # tmax is the largest time value
  # tslice gives the relevant time slices, the form should be: c(lb1,ub1,lb2,ub2, ...)
  
  nsim    <- ntrain + ntest
  tsteps  <- seq(min(tslice), max(tslice), length=p)
 
  ## Time slices
  J <- findInterval(tsteps, tslice[-c(1,length(tslice))])+1
  m <- table(J)
  nslices <- length(m)
  
  ## Index of data in relevant slices
  if (is.list(relevant.slices)){
    II <- list()
    for (i in 1:length(relevant.slices)) {
      II[[i]] <- J %in% relevant.slices[[i]]
    }
    relevant.slices <- unique(unlist(relevant.slices))
    I <- findInterval(tsteps, tslice) %in% relevant.slices
  } else {
    I <- J %in% relevant.slices
    II <- list(I)
  }
  
  color <- rep(c("black", "grey"), ceiling(length(m)/2))[1:length(m)]
  color <- rep(color, m)
  lwd   <- rep(1, p)
  
  nicecolor <- rainbow(length(relevant.slices))[sample.int(length(relevant.slices))]
  color[I]  <- rep(nicecolor, m[relevant.slices])
  lwd[I]    <- 3
  
  if (pb == 1){
    ## PROBLEM 1 - BASED ON GP TRAJECTORIES - single beta
    sigma <- .1
    sigma_eps <- 0.01
    model <- km(formula=~x+I(x^2), design=data.frame(x=tsteps), response=rep(0,p), 
                covtype="matern3_2", coef.trend=c(intercept=-5, beta1=4, beta2=-4), 
                coef.cov=.2/sqrt(3), coef.var=sigma)
    Sigma_x <- chol2inv(model@T)
    x <- simulate(model, nsim=nsim, newdata=NULL) + matrix(rnorm(nsim*p, sd=sigma_eps), nsim, p)
    
    # generate y
    beta <- beta0 <- matrix(0, p, length(II))
    y <- 0
    for (i in 1:length(II)) {
      beta[,i] <- beta0[,i] <- sin(tsteps*(2+i)*pi/2 - (i-1)*pi/3)
      
      y <- y + log(abs(as.vector(x[,II[[i]]]%*%matrix(beta[II[[i]],i], ncol=1) )))
      beta0[!II[[i]],i] <- 0
    }
    # y <- log(abs(y))
    
  } else if (pb == 2) {
    ## PROBLEM 2 - BASED ON BROWNIAN MOTION TRAJECTORIES
    sigma <- 1
    sigma_eps <- 0.1
    simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=sigma))))
    x <- t(replicate(nsim, simulate_bm()))
    Sigma_x <- NULL
    
    # generate beta_1 and beta_2
    beta <- beta0 <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2))
    
    if(length(relevant.slices)>2) II <- list(I, I)
    
    y <- as.vector(log(exp(beta[II[[1]],1] %*% t(x[ ,II[[1]]])) +
                         exp(abs(beta[II[[2]],2] %*% t(x[ ,II[[2]]])))) +
                     rnorm(nsim, sd=sigma_eps))
    beta0[!II[[1]],1] <- 0
    beta0[!II[[2]],2] <- 0
    
  } 
  
  ## Split data between training and testing
  X     <- x[1:ntrain,]
  Y     <- y[1:ntrain]
  Xtest <- x[ntrain + 1:ntest,]
  Ytest <- y[ntrain + 1:ntest]
  
  ## Plots
  if (draw){
    par(mfrow=c(2,2))
    # Draw X
    plot(tsteps, X[1,], type="l", ylim=range(X), col=color, xlab="t", ylab="X")
    for (i in 1:ntrain) 
      segments(x0=tsteps[-p], y0=X[i,-p], x1=tsteps[-1], y1=X[i,-1], col=color)
    
    # Draw beta
    if (ncol(beta)==1) {
      plot(tsteps, beta[,1], type="l", xlab="t", ylab="beta", main="Beta", 
           col="grey")
      segments(x0=tsteps[-p], y0=beta0[-p,1], x1=tsteps[-1], y1=beta0[-1,1], 
               col=color, lwd=lwd)
    } else {
      plot(NULL, xlab="t", ylab="beta", main="Beta", xlim=range(tsteps), 
           ylim=range(beta))
      for (i in 1:ncol(beta)) {
        colori <- rep("black", p);  colori[II[[i]]]  <- nicecolor[i]
        lwdi <- rep(1, p); lwdi[II[[i]]] <- 3
        lines(tsteps, beta[,i], col="grey")
        segments(x0=tsteps[-p], y0=beta0[-p,i], x1=tsteps[-1], y1=beta0[-1,i], 
                 col=colori, lwd=lwdi)
      }
    }
    
    # Draw X | y
    H <- 10
    slices_bounds <- quantile(Y, probs=seq(0, 1, length=H+1))
    Y_class <- cut(Y, breaks=slices_bounds, labels=FALSE, include.lowest=TRUE)
    colH <- rainbow(length(slices_bounds))
    hist(Y, breaks=slices_bounds, col=colH, main="Y slices", xlab="y")
    
    plot(tsteps, X[1,], type="l", ylim=range(X), col=colH[Y_class[1]], xlab="t",
         ylab="X", main="X | y (by slice)", lwd=lwd*2/3)
    for (i in 1:ntrain) segments(x0=tsteps[-p], y0=X[i,-p], x1=tsteps[-1], 
                                 y1=X[i,-1],col=colH[Y_class[i]], lwd=lwd*2/3)
    
  }
  return(list(tsteps=tsteps, X=X, Y=Y, Xtest=Xtest, Ytest=Ytest, m=m, 
              color=color, lwd=lwd, beta0=beta0, beta=beta, Sigma_x=Sigma_x,
              D=nslices))
}