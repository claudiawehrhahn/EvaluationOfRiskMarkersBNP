
# Functions to be used

# selecting scaling factors:
#   arguments:
#   X:            design matrix
#   xbeta0:       value defining an interval where the linear predictor lies with high probability
#   var0:         initial value for the scaling constant
#   varIncrease:  how much to increase the scaling constant in each iteration
tunningBetafun <- nimbleFunction(
  run = function(X = double(2), xbeta0 = double(0), var0 = double(0), varIncrease = double(0))
  {
    returnType(double(0))
    
    # sample size and (X^t%*%X)^-1
    nrec <- dim(X)[1]
    xtx1 <- inverse(t(X)%*%X)
    
    # objects to save scaling factors
    aux <- numeric(nrec)
    ans <- 0
    
    for(ix in 1:nrec) { #  select a scaling constant for each value of the predictor, x_i
      xi <- X[ix, ]  
      varAux <- t(xi) %*% xtx1 %*% xi 
      xVar <- varAux[1,1] # variance of x_i^t %*% beta
      
      tau2 <- var0 # initializing scaling constant. This value is increased until P(-xbeta0 < x^t*beta < xbeta0) >= 0.95
      probCond <- pnorm(xbeta0, 0, sd = sqrt(tau2 * xVar)) - pnorm(-xbeta0, 0, sqrt(tau2 * xVar))
      while(probCond >= 0.95) {
        tau2 <- tau2 + varIncrease
        probCond <- pnorm(xbeta0, 0, sqrt(tau2 * xVar)) - pnorm(-xbeta0, 0, sqrt(tau2 * xVar))
      }
      
      aux[ix] <- tau2
    }
    ans <- min(aux)
    return(ans)
  }
)
tunningBetafunC <- compileNimble(tunningBetafun)

# computing the dependent correlation: rho_j(x)
#   arguments:
#   beta: current value of beta^{rho}_{j}
#   xi:   value of the predictor
hZfun <- function(beta, xi){
  aux <- abs(xi%*%beta)
  return(2/(aux+1)-1)
}


# computing the dependent stick-variables: v_j(x)
#   arguments:
#   beta: current value of beta^{rho}_{j}
#   xi:   value of the predictor
hEtafun <- function(beta, xi){
  aux <- xi%*%beta 
  return(expit(aux))
}

# computing likelihood in log scale
#   arguments:
#   y:    data matrix
#   wj:   matrix of weights
#   rhoj: matrix of correlations
nimlik=nimbleFunction(
  run=function(y=double(2), wj=double(2), rhoj=double(2))
  {
    returnType(double(0))
    n <- dim(y)[1]
    L <- dim(wj)[2]
    ans <- 0
    for( i in 1:n){
      y1 <- y[i,1]
      y2 <- y[i,2]
      suml <- 0
      for(l in 1:L) {
        rho2 <- rhoj[i,l]*rhoj[i,l]
        det0 <- 1 - rho2
        suml <- suml + wj[i, l]*abs(det0)^(-0.5)*exp(-0.5*(rho2*(y1^2+y2^2)-2*y1*y2*rhoj[i,l])/det0)
      }
      ans = ans + log(suml)
    }
    return(ans)
  }
)
cnimlik <- compileNimble(nimlik)

# function transforming mcmc list output of the coefficients defining the stick breaking variables
# into a matrix
#   arguments:
#   x:      mcmc output as a vector
#   p:      number of predictors
#   nSB:    truncation level of the stick breaking representation
#   nsave:  number of saved mcmc iterations
betaVMatrix <- nimbleFunction(
  run = function(x = double(1), p = integer(0), nSB= integer(0), nsave= integer(0)) {
    returnType(double(2))
    
    nSB2 <- nSB - 1
    out <- matrix(0, ncol=nSB2, nrow=p*nsave)
    
    ii <- 1
    for(i in 1:nsave) {
      for(j in 1:nSB2) {
        for(k in 1:p) {
          out[(i-1)*p+k, j] <- x[ii]
          ii <- ii + 1  
        }
      }
    }
    return(out)
  }
) 
cbetaVMatrix <- compileNimble(betaVMatrix)


# function transforming mcmc list output of the coefficients defining the correlations into a matrix
#   arguments:
#   x:      mcmc output as a vector
#   p:      number of predictors
#   nSB:    truncation level of the stick breaking representation
#   nsave:  number of saved mcmc iterations
betaRMatrix <- nimbleFunction(
  run = function(x = double(1), p = integer(0), nSB= integer(0), nsave= integer(0)) {
    returnType(double(2))
    
    out <- matrix(0, ncol=nSB, nrow=p*nsave)
    
    ii <- 1
    for(i in 1:nsave) {
      for(j in 1:nSB) {
        for(k in 1:p) {
          out[(i-1)*p+k, j] <- x[ii]
          ii <- ii + 1  
        }
      }
    }
    return(out)
  }
) 
cbetaRMatrix <- compileNimble(betaRMatrix)


# function that computes posterior samples of copula density and copula cdf predictives for a given predictor, 
# for every mcmc iteration
#   arguments:
#   betaVSamAux:  matrix version of mcmc output of coefficients defining the stick-breaking variables
#   betaRSamAux:  matrix version of mcmc output of coefficients defining the correaltions
#   x0:           predictor
#   invNorm:      quantiles of the standard normal distribution on a grid of (0,1)
#   nsave:        number of saved mcmc iterations
compCondCopulaDenCDF <- nimbleFunction(
  run = function(betaVSamAux = double(2), betaRSamAux = double(2), x0 = double(1), invNorm = double(2),
                 nsave = integer(0)) {
    returnType(double(2))
    
    # size of the grid, number of predictors, and truncation level of stick-breaking representation
    ngrid <- dim(invNorm)[1]
    p <- length(x0)
    nSB = dim(betaRSamAux)[2]
    nSB2 <- nSB - 1
    
    # objects where stick-breaking variables con correlations are computed at each iteration of mcmc
    v <- numeric(nSB2)
    rho <- numeric(nSB)
    
    # output matrix. First ngrid columns save the copula density and colums (ngrid+1):(2*ngrid) save
    # the copula cdf. These computationare save per row block for each ietration of mcmc.
    out <- matrix(0, ncol=(2*ngrid), nrow=(nsave*ngrid))
    
    
    for(iter in 1:nsave) { # for each ietartion of mcmc
      for(l in 1:nSB2) { # compting stick-breaking variables for predictor x0
        v[l] <- 0
        for(k in 1:p) {
          v[l] <- v[l] + betaVSamAux[(iter-1)*p+k , l] * x0[k]
        }
      }
      v <- exp(v) / (1 + exp(v))
      w <- stick_breaking(v)
      
      for(l in 1:nSB) { # computing correlations at predictor x0
        rho[l] <- 0
        for(k in 1:p) {
          rho[l] <- rho[l] + betaRSamAux[(iter-1)*p+k , l] * x0[k]
        }
      }
      rho <- abs(rho)
      rho <- 2/(rho+1)-1
      
      for(i in 1:ngrid){ # computing posterior sample of copula density at values of a grid for predictor x0
        for(j in 1:ngrid){
          y1 <- invNorm[i,1]
          y2 <- invNorm[j,2]
          rho2 <- rho*rho
          det0 <- 1-rho2
          out[(iter-1)*ngrid+i , j]  <- sum(w*abs(det0)^(-0.5)*exp(-0.5*(rho2*(y1^2+y2^2)-2*y1*y2*rho)/det0))
        }
      }
      # normalizing copula density computed on grid and saved in colums 1:ngrid
      out[((iter-1)*ngrid+1):(iter*ngrid) , 1:ngrid]  <- out[((iter-1)*ngrid+1):(iter*ngrid) , 1:ngrid]  / sum(out[((iter-1)*ngrid+1):(iter*ngrid) , 1:ngrid] )
      
      # computing copula cdf based on computations of copula density and saved in colums (ngrid+1):(2*ngrid)
      sumCol <- matrix(0, ncol=ngrid, nrow=ngrid)
      for(j in 1:ngrid) {
        sumCol[1, j] <- out[(iter-1)*ngrid+1 , j]
        for(i in 2:ngrid) {
          sumCol[i, j] <- sumCol[i-1,j] + out[(iter-1)*ngrid+i , j]
        }
      }
      sumRow <- matrix(0, ncol=ngrid, nrow=ngrid)
      for(i in 1:ngrid) {
        sumRow[i, 1] <- sumCol[i, 1]
        for(j in 2:ngrid) {
          sumRow[i, j] <- sumRow[i, j-1] +  sumCol[i, j]
        }
      }
      out[((iter-1)*ngrid+1):(iter*ngrid) , (ngrid+1):(2*ngrid)] <- sumRow
      
      
    }
    return(out)
  }
)
ccompCondCopulaDenCDF <- compileNimble(compCondCopulaDenCDF)


# function that computes conditional copula density estimate for given predictor
#   arguments:
#   aux:          posterior samples of conditional copula density
#   nsave:        number of saved mcmc iterations
compCondDenHat <- nimbleFunction(
  run = function(aux = double(2), nsave = integer(0)) {
    returnType(double(2))
    ngrid <- dim(aux)[2]/2
    
    out <- matrix(0, ncol=ngrid, nrow=ngrid)
    
    # averaging over mcmc iterations
    for(iter in 1:nsave) {
      out <- out + aux[((iter-1)*ngrid+1):(iter*ngrid) , 1:ngrid]
    }
    out <- out / nsave
    
    return(out)
  }
)
ccompCondDenHat <- compileNimble(compCondDenHat)


# computing conditional Kendall's tau estimate for given predictor
#   arguments:
#   auxCopulaDenCDF:  posterior samples of conditional copula density and cdf
#   nsave:            number of saved mcmc iterations
compCondTau <- nimbleFunction(
  run = function(auxCopulaDenCDF = double(2), nsave = integer(0)) {
    returnType(double(2))
    ngrid <- dim(auxCopulaDenCDF)[2]/2
    out <- matrix(0, ncol=1, nrow=nsave)
    
    for(i in 1:nsave) {
      out[i, 1] <- 4 * sum( auxCopulaDenCDF[((i-1)*ngrid+1):(i*ngrid) , (ngrid+1):(2*ngrid)] * auxCopulaDenCDF[((i-1)*ngrid+1):(i*ngrid) , 1:ngrid] ) - 1
    }
    return(out)
  }
)
ccompCondTau <- compileNimble(compCondTau)



#----------------------------------------------------------------------------------


# function that samples the model
#   arguments:
#   U:          matrix with pseudo-data
#   X:          design matrix
#   nburn:      burn in period of mcmc iterations
#   njump:      thinning in mcmc iterations
#   nsave:      number of saved mcmc iterations
#   SimData:    whether simulated or real pseudo-data is used
condCopulaDensitySampling <- function(U, X, nburn, njump, nsave, SimData) {
  
  # sample size, number of predictors, dimension of U
  nrec <- nrow(U)
  p <- ncol(X)
  d <- ncol(U)
  
  #---- hyperparameters:
  muEta <- rep(0,p)
  muZ <- muEta
  if(SimData == 1) {      # hyperprior for simulated data
    SEta <- diag(2.25, p) 
    SZ <- diag(2.25, p) 
  } else {                # hyperprior for real data
    SEta <- matrix(0, ncol=p, nrow=p)
    SZ <- SEta
    
    # scaling constants of priors of beta^V and beta^Z
    tunningBetaV1 <- tunningBetafunC(X[, 1:3], 6, 0.01, 0.01) 
    tunningBetaRho1 <- tunningBetafunC(X[, 1:3], 6, 1, 1)
    
    SEta[1:3, 1:3] <- tunningBetaV1 * solve(t(X[, 1:3])%*%X[, 1:3])
    SEta[4:7, 4:7] <- 2.25 * solve(t(X[, 4:7])%*%X[, 4:7]) 
    SZ[1:3, 1:3] <- tunningBetaRho1 * solve(t(X[, 1:3])%*%X[, 1:3])
    SZ[4:7, 4:7] <- SEta[4:7, 4:7]
  }
  SEtai <- solve(SEta)
  SZi <- solve(SZ)
  
  #---- initializing parameters beta^V and beta^Z:
  maxN <- 25
  set.seed(1)
  bEta <- t(rmvnorm(maxN-1, muEta, SEta ))
  bZ <- t(rmvnorm(maxN, muZ, SZ ))
  
  #---- slice sampler parameter:
  wsliceV <- rep(2, p)
  
  
  #---- MCMC specification:
  ntotal <- nburn + nsave*njump
  
  
  #--- Storagingn objects
  postbZ <- list()
  postbEta <- list()
  llik <- c()
  
  # standard normal quantiles
  invNorm <- qnorm(U) 
  set.seed(1)
  #---- Algorithm
  for(niter in 1:ntotal){
    
    
    #-- update beta^V with a multivariate slice sampler:
    
    R <- t(sapply(1:nrec, function(i)hZfun(bZ,X[i,]))) # correlation matrix is fixed when updating bEta
    for(j1 in 1:(maxN-1)){
      x0 <- bEta[,j1] # current value
      
      V <- t(sapply(1:nrec, function(i)hEtafun(bEta,X[i,])))                                            # stick-breaking variables
      Vaux <- t(apply(V, 1, function(i)cumprod(1-i)))
      W <- t(sapply(1:nrec, function(i)c(V[i,1], V[i,2:(maxN-1)]*Vaux[i,-(maxN-1)], Vaux[i,maxN-1] )))  # stick-breaking weights
      fx0 <- cnimlik(invNorm, W, R) - 0.5*t(x0-c(muEta))%*%SEtai%*%(x0-c(muEta))                        # posterior in log scale
      
      #-- step a)
      yslice <- fx0 - rexp(1, 1)            # yslice in log-scale
      #-- step b)
      Li <- bEta[,j1] - wsliceV * runif(p)  # lower limits of hyper rectangle
      Ri <- Li + wsliceV                    # upper limits of hyper rectangle
      #-- step c)
      x1 <- Li + runif(p)*(Ri -Li)          # candidate
      bEta[,j1] <- x1
      
      V <- t(sapply(1:nrec, function(i)hEtafun(bEta, X[i,])))                                           # updated stick-breaking variables
      Vaux <- t(apply(V, 1, function(i)cumprod(1-i)))
      W <- t(sapply(1:nrec, function(i)c(V[i,1], V[i,2:(maxN-1)]*Vaux[i,-(maxN-1)], Vaux[i,maxN-1] )))  # updated stick-breaking weights
      fx1 <- cnimlik(invNorm, W, R) - 0.5*t(x1-c(muEta))%*%SEtai%*%(x1-c(muEta))                        # updated posterior in log scale
      
      while(yslice > fx1){  # slice sampler loop
        for(j2 in 1:p) {              # updating hyper rectangle
          if( x1[j2] < x0[j2]) { Li[j2] <- x1[j2] }  else { Ri[j2] <- x1[j2] }  
        }
        
        x1 <- Li + runif(p)*(Ri -Li)  # updated candidate
        bEta[,j1] <- x1
        V <- t(sapply(1:nrec, function(i)hEtafun(bEta, X[i,])))                                           # updated stick-breaking variables
        Vaux <- t(apply(V, 1, function(i)cumprod(1-i)))
        W <- t(sapply(1:nrec, function(i)c(V[i,1], V[i,2:(maxN-1)]*Vaux[i,-(maxN-1)], Vaux[i,maxN-1] )))  # updated stick-breaking weights
        fx1 <- cnimlik(invNorm, W, R) - 0.5*t(x1-c(muEta))%*%SEtai%*%(x1-c(muEta))                        # updated posterior in log scale
      }
    }
    
    
    #-- update betaZ:
    
    V <- t(sapply(1:nrec, function(i)hEtafun(bEta, X[i,])))
    Vaux <- t(apply(V, 1, function(i)cumprod(1-i)))
    W <- t(sapply(1:nrec, function(i)c(V[i,1], V[i,2:(maxN-1)]*Vaux[i,-(maxN-1)], Vaux[i,maxN-1] ))) # weights matrix is fixed when updating bZ
    
    for(j1 in 1:maxN){
      x0 <- bZ[,j1]   # current value
      
      R <- t(sapply(1:nrec, function(i)hZfun(bZ,X[i,])))                    # correlations
      fx0 <- cnimlik(invNorm, W, R) - 0.5*t(x0-c(muZ))%*%SZi%*%(x0-c(muZ))  # posterior in log scale
      
      #-- step a)
      yslice <- fx0 - rexp(1, 1)          # yslice in log-scale
      #-- step b)
      Li <- bZ[,j1] - wsliceV * runif(p)  # lower limits of hyper rectangle
      Ri <- Li + wsliceV                  # upper limits of hyper rectangle
      #-- step c)
      x1 <- Li + runif(p)*(Ri -Li)        # candidate
      bZ[,j1] <- x1
      
      R <- t(sapply(1:nrec, function(i)hZfun(bZ,X[i,])))                    # updated correlations
      fx1 <- cnimlik(invNorm, W, R) - 0.5*t(x1-c(muZ))%*%SZi%*%(x1-c(muZ))  # updated posterior in log scale
      
      while(yslice>fx1){    # slice sampler loop
        for(j2 in 1:p) {                  # updating hyper rectangle
          if( x1[j2] < x0[j2]) { Li[j2] <- x1[j2] }  else { Ri[j2] <- x1[j2] }  
        }
        x1  <-  Li + runif(p)*(Ri -Li)    # updated candidate
        bZ[,j1] <- x1
        R <- t(sapply(1:nrec, function(i)hZfun(bZ,X[i,])))                    # updated correlations
        fx1 <- cnimlik(invNorm, W, R) - 0.5*t(x1-c(muZ))%*%SZi%*%(x1-c(muZ))  # updated posterior in log scale
      }
    }
    
    llik[niter] <- cnimlik(invNorm, W, R) # computing likelhood in log scale
    
    postbZ[[niter]] <- bZ       # saving updated parameters for beta^Z
    postbEta[[niter]] <- bEta   # saving updated parameters for beta^V
    
    if(trunc(niter/100)==niter/100){print(niter)}
  }
  
  itersave=nburn+(1:nsave)*njump
  Eta=list()
  for(i in 1:nsave){Eta[[i]]=postbEta[[itersave[i]]]}
  Z=list()
  for(i in 1:nsave){Z[[i]]=postbZ[[itersave[i]]]}
  
  list(betaV  = Eta, betaR = Z, llik=llik)
}


# function that computes conditional copula density estimates and conditional kendall's tau for real data set.
# Here specific values of the covariates are considered.
#   arguments:
#   ix3:        index specifying which value in the grid of (gender, age) to use
#   x1grid:     grid for tg4 (tryglicerides)
#   ix1:        index specifying which value in the grid of x1grid to use
#   x2grid:     grid for bmi4 (body mass index)
#   ix2:        index specifying which value in the grid of x2grid to use
#   nsave:      number of saved mcmc iterations
#   mcmc:       mcmc output
#   invNorm:    quantiles of the normal distribution on a grid of the (0,1) interval
condCopulaAndAssociationAppli <- function(ix3, x1grid, ix1, x2grid, ix2, nsave, mcmc, invNorm) {
  if(ix3 == 1) {x0 <- c(1, x1grid[ix1], x2grid[ix2], 0, 0, 0, 1)} # (intercept, tryglicerides, body mass index, male, age < 55)
  if(ix3 == 2) {x0 <- c(1, x1grid[ix1], x2grid[ix2], 0, 0, 0, 0)} # (intercept, tryglicerides, body mass index, female, age < 55)
  if(ix3 == 3) {x0 <- c(1, x1grid[ix1], x2grid[ix2], 1, 0, 0, 1)} # (intercept, tryglicerides, body mass index, male, 55 <= age < 60)
  if(ix3 == 4) {x0 <- c(1, x1grid[ix1], x2grid[ix2], 1, 0, 0, 0)} # (intercept, tryglicerides, body mass index, female, 55 <= age < 60)
  if(ix3 == 5) {x0 <- c(1, x1grid[ix1], x2grid[ix2], 0, 1, 0, 1)} # (intercept, tryglicerides, body mass index, male, 65 <= age < 65)
  if(ix3 == 6) {x0 <- c(1, x1grid[ix1], x2grid[ix2], 0, 1, 0, 0)} # (intercept, tryglicerides, body mass index, female, 65 <= age < 65)
  if(ix3 == 7) {x0 <- c(1, x1grid[ix1], x2grid[ix2], 0, 0, 1, 1)} # (intercept, tryglicerides, body mass index, male, age >= 65)
  if(ix3 == 8) {x0 <- c(1, x1grid[ix1], x2grid[ix2], 0, 0, 1, 0)} # (intercept, tryglicerides, body mass index, female, age >= 65)
  
  # number of predictors and truncation level of random measure
  p <- length(x0)
  nSB <- ncol(mcmc$betaR[[1]])
  
  # mcmc outputs transformed into matrices
  betaVSamAux <- cbetaVMatrix(unlist(mcmc$betaV), p, nSB, nsave)
  betaRSamAux <- cbetaRMatrix(unlist(mcmc$betaR), p, nSB, nsave)
  
  auxCopulaDenCDF <- ccompCondCopulaDenCDF(betaVSamAux, betaRSamAux, x0, invNorm, nsave)  # posterior samples of coopula density and cds for given prodictor
  copulaDenHat <- ccompCondDenHat(auxCopulaDenCDF, nsave)                                 # computing copula density estimate for given predictor
  tauRho <- ccompCondTau(auxCopulaDenCDF, nsave)                                          # posterior samples of kendall's tau for given predictor
  
  print(paste("E=", E, "ix1 =", ix1, "ix2 =", ix2, "ix3 =", ix3))
  
  copula <- list(condCopulaDenHat = copulaDenHat) # copula estimate output for given predictor
  association <- list(condTau = tauRho[, 1])      # kendalls tau output for given preditor
  # saving output 
  save(copula, file=paste("./Results/Application/condCopula_E=", E, "_ix1=", ix1, "_ix2=", ix2, "_ix3=", ix3, ".RData", sep=""))
  save(association, file=paste("./Results/Application/condAssociation_E=", E, "_ix1=", ix1, "_ix2=", ix2, "_ix3=", ix3, ".RData", sep=""))
  
}

# function that computes conditional copula density estimates and conditional kendall's tau for sunthetic data sets.
# Here specific values of the covariates are considered.
#   arguments:
#   ix1:        index specifying which value in the grid of the (0,1) interval to use
#   x1grid:     grid of the (0,1) interval
#   nsave:      number of saved mcmc iterations
#   mcmc:       mcmc output
#   invNorm:    quantiles of the normal distribution on a grid of the (0,1) interval
condCopulaAndAssociationSim <- function(ix1, x1grid, nsave, mcmc, invNorm) {
  x0 <- c(1, x1grid[ix1]) # interpcept and value for predictor
  
  # number of predictors and truncation level of random measure
  p <- length(x0)
  nSB <- ncol(mcmc$betaR[[1]])
  
  # mcmc outputs transformed into matrices
  betaVSamAux <- cbetaVMatrix(unlist(mcmc$betaV), p, nSB, nsave)
  betaRSamAux <- cbetaRMatrix(unlist(mcmc$betaR), p, nSB, nsave)
  
  auxCopulaDenCDF <- ccompCondCopulaDenCDF(betaVSamAux, betaRSamAux, x0, invNorm, nsave)  # posterior samples of coopula density and cds for given prodictor
  copulaDenHat <- ccompCondDenHat(auxCopulaDenCDF, nsave)                                 # computing copula density estimate for given predictor
  tauRho <- ccompCondTau(auxCopulaDenCDF, nsave)                                          # posterior samples of kendall's tau for given predictor
  
  print(paste("E=", E, "ix1 =", ix1))
  
  copula <- list(condCopulaDenHat = copulaDenHat)   # copula estimate output for given predictor
  association <- list(condTau = tauRho[, 1])        # kendalls tau output for given preditor
  # saving output 
  save(copula, file=paste("./Results/Simulation/condCopula_E=", E, "_ix1=", ix1, ".RData", sep=""))
  save(association, file=paste("./Results/Simulation/condAssociation_E=", E, "_ix1=", ix1, ".RData", sep=""))
  
}
