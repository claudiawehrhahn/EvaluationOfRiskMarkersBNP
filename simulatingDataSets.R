
rm(list=ls())

library(copula)
library(parallel)

#-- simulating data sets:
# Simulation of data sets for Scenario I and Scenario II.
# Data sets are simulated from the joint distribution defined by means of a copula function and
# marginal distributions. In both cases marginals are univariate standard normal distributions.

# Scenario I: the copula function is given by a t-copula with correlation function equal to x, where x denotes
#               the predictor.
# Scenario II: the copula function is given by a mixture of copula functions. With probability x 
#               a t-copula with correlation -x(1-x)^2 and with  prob (1-x) a gumbel copula with
#               parameter x^2(1-x)+1, where x denotes the predictor.


simulatedData <- function(E) {
  
  set.seed(1)
  #-- data simulation:
  # sample size
  n <- 250
  # predictor
  x <- runif(n, 0,1)
  
  # object where data from the joint is saved
  uaux <- c()
  # weights in the mixture for Scenario II
  unif <- runif(n)
  
  if(E==1) { # Scenario I
    for(i in 1:n) {
      cop <- tCopula(x[i], dim=2)  # copula function
      joint  <- mvdc(cop, c("norm", "norm"),
                     param = list( list(mean = 0, sd = 1), list(mean = 0, sd = 1))) # joint distribution
      uaux <- rbind(uaux, rMvdc(1, joint))  # simulation
    }
  }
  if(E == 2) { # Scenario II
    for(i in 1:n) {
      if(unif[i]  <= x[i]) { 
        cop <- tCopula(-x[i]*(1-x[i])^2 , dim=2) # copula function
        joint  <- mvdc(cop, c("norm", "norm"),
                       param = list( list(mean = 0, sd = 1), list(mean = 0, sd = 1))) # joint distribution
        uaux <- rbind(uaux, rMvdc(1, joint))  # simulation
      } else {
        cop <- gumbelCopula(x[i]^2*(1-x[i]) + 1, dim=2) # copula function
        joint  <- mvdc(cop, c("norm", "norm"),
                       param = list( list(mean = 0, sd = 1), list(mean = 0, sd = 1))) # joint distribution
        uaux <- rbind(uaux, rMvdc(1, joint)) # simulation
      }
    }
  }
  
  # pseudo observations and predictors
  u <- cbind(rank(uaux[,1])/(n+1), rank(uaux[,2])/(n+1))  
  X <- x
  
  # equally spaced grids for computing copula functions and predictors 
  uGrid <- seq(0.01,0.99,by=0.01)
  nGrid <- length(uGrid)
  predGrid <- seq(0.01,0.99,by=0.01)

  # objects where density of copula, CDF of copula, and Kendall's tau are saved
  dTrueCopulas <- list()
  pTrueCopulas <- list()
  kendallsTau <- c()
  
  # computing copula densities, copula CDF, and Kendall's tau for each value of the predictor:
  for(ix in seq_along(predGrid)) {
    dTrueCopulas[[ix]] <- matrix(0, ncol = nGrid, nrow = nGrid)
    pTrueCopulas[[ix]] <- matrix(0, ncol = nGrid, nrow = nGrid)

    if(E == 1) {
      info <- "one t-copula with correlation function equal to x"
      cop <- tCopula(predGrid[ix], dim=2) 
      for(i in seq_along(uGrid)) {
        for( j in seq_along(uGrid) ) {
          dTrueCopulas[[ix]][i,j] <- dCopula(c(uGrid[i], uGrid[j]), cop) # copula density
          pTrueCopulas[[ix]][i,j] <- pCopula(c(uGrid[i], uGrid[j]), cop) # copula cdf
        }
      }
    }
    if(E == 2) { # mixture of 2 copulas.
      info <- "mixture of 2 copulas: one  t-copula with probability x with correlation
      -x(1-x)^2, one gumbel copula with prob (1-x) with parameter x^2(1-x)+1"
      
      cop1 <- tCopula(-predGrid[ix]*(1-predGrid[ix])^2 , dim=2)
      cop2 <- gumbelCopula(predGrid[ix]^2*(1-predGrid[ix]) + 1, dim=2)
      for(i in 1:nGrid) {
        for( j in 1:nGrid ) {
          dTrueCopulas[[ix]][i,j] <- predGrid[ix] * dCopula(c(uGrid[i], uGrid[j]), cop1) + 
            (1-predGrid[ix]) * dCopula(c(uGrid[i], uGrid[j]), cop2) # copula density
          pTrueCopulas[[ix]][i,j] <- predGrid[ix] * pCopula(c(uGrid[i], uGrid[j]), cop1)+
            (1-predGrid[ix])*pCopula(c(uGrid[i], uGrid[j]), cop2) # copula cdf
        }
      }
    } 
    
    dTrueCopulas[[ix]] <- dTrueCopulas[[ix]] / sum(dTrueCopulas[[ix]]) # normalizing copula density
    
    kendallsTau[ix] <- 4 * sum( pTrueCopulas[[ix]] * dTrueCopulas[[ix]] ) - 1 # true kendall's tau
    print(ix)
  }
  
  # output
  results <- list(uGrid = uGrid, nGrid = nGrid, predGrid = predGrid, info = info,
                  dTrueCopulas = dTrueCopulas,
                  pTrueCopulas = pTrueCopulas,
                  kendallsTau = kendallsTau,
                  u = u,
                  X = X)
  # saving output
  save(results, file=paste("./Data/simulationScenario", E, ".RData", sep=""))
}
  
simulatedData(1)

  