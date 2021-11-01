rm(list=ls())

library(copula)
library(nimble)
library(mvtnorm)
library(parallel)
library(ggplot2)

source("./Functions.R")


# Test Code: runns a very short MCMC for both simulated and real data sets

#------------------------------
#-- Simulation

#------------------------------------------------------------------------------------------------
#-- grids for computing conditional copula and Kendall tau's estimates

# grid for conditional copula densities
grid01 <- seq(0.01, 0.99, by = 0.1) 
grid01 <- cbind(grid01, grid01)
ngrid01 <- nrow(grid01)

# grid for covariate:
xgrid <- seq(0.01, 0.99, by = 0.1)
ngrid <- length(xgrid)

# normal quantiles of the (0,1)-grid used to compute conditional copula densities
invNorm <- qnorm(grid01)

#------------------------------------------------------------------------------------------------
#-- runnig MCMC

# mcmc specification
nburn <- 10    # burn-in period
njump <- 1       # size of thinning
nsave <- 10    # number of saved iterations
nCores <- 20      # number of cores to be used in posterior inference
E <- 1

#-- loading data set:
load(file=paste("./Data/simulationScenario", E, ".RData", sep=""))

#-- definig variables, design matrix, and sample size
U <- results$u
nrec <- nrow(U)
X <- cbind(rep(1, nrec), results$X)

#-- mcmc run
mcmc <- condCopulaDensitySampling(U, X, nburn, njump, nsave, SimData = 1)  # 10 hours per data set
save(mcmc, file=paste("./Results/Simulation/mcmcResults_E", E, ".RData", sep=""))

#-- computing conditional copula estimates and conditional Kendall's tau
mclapply(1:ngrid, condCopulaAndAssociationSim, x1grid = xgrid, 
         nsave = nsave, mcmc = mcmc, invNorm = invNorm, mc.cores =  nCores)



#--------------------------------------------------------------------------------
#-- plots for simulated data sets:

# grid of covariate:
xgrid <- seq(0.01, 0.99, by = 0.1) # seq(0.01, 0.99, by = 0.01)
ngrid <- length(xgrid)


#-- true kandalls's tau for each simulated scenario:
load(file=paste("./Data/simulationScenario", 1, ".RData", sep=""))
trueKendallsTauE1 <- results$kendallsTau



#-- estimating conditional Kendall's tau and 95% credible band for each simulated scenario:

# load the results for Kendall's tau 
nsave <- 10 # 10000 
kendallsTauE1 <- matrix(0, nrow = ngrid, ncol = nsave) # Kendall's tau for E=1
kendallsTauE2 <- matrix(0, nrow = ngrid, ncol = nsave) # Kendall's tau for E=2
for(ix1 in 1:ngrid) {
  load(file=paste("./Results/Simulation/condAssociation_E=", 1, "_ix1=", ix1, ".RData", sep=""))
  kendallsTauE1[ix1, ] <- association$condTau
}

# computing quantiles (estimate and 95% credible band)
kendallsTauE1Quantiles <- matrix(0, ncol=ngrid, nrow=3)
kendallsTauE1Quantiles[1, ] <- apply(kendallsTauE1, 1, function(x)quantile(x, 0.025))
kendallsTauE1Quantiles[2, ] <- apply(kendallsTauE1, 1, function(x)quantile(x, 0.5))
kendallsTauE1Quantiles[3, ] <- apply(kendallsTauE1, 1, function(x)quantile(x, 0.975))


#-- plots of kendall's tau for  scenario I:
# limits of plots
yKendallsTau <- c(-0.3, 0.9)

# accomodating outputs for using ggplot 
y <- c(kendallsTauE1Quantiles[1, ], rev(kendallsTauE1Quantiles[3, ]))
x <- c(xgrid, rev(xgrid))
hat <- c()
true <- c()
x2 <- c()
for(i in 1:ngrid) {
  hat <- c(hat, rep(kendallsTauE1Quantiles[2, i], 2))
  true <- c(true, rep(trueKendallsTauE1[i], 2))
  x2 <- c(x2, rep(xgrid[i], 2))
}

tests <-  data.frame(y, x, hat, true, x2)
colnames(tests) <- c("Correlation", "x", "Hat", "True", "x2")

pdf(file=paste("./Results/Simulation/curve_kendallsTau_E=",E, ".pdf",sep=""))
par(cex=2, mar=c(2.3, 2.3, 1, 1), mgp=c(2.2,1,0))
ggplot(tests, aes(x = x, y = Correlation))+ 
  geom_polygon(data = tests, alpha = 0.3) +
  geom_line(aes(x = x2, y = Hat), lwd=2, lty=2) +
  geom_line(aes(x = x2, y = True), lwd=2, lty=1) +
  xlim(0,1) +
  ylim(-0.3, 1) +
  theme_grey(base_size = 30) +
  theme(legend.position = "none")+ 
  scale_color_manual(values = "grey20")
dev.off()


#-- computing L1 distance:

kendallsTauE1Estimate <- apply(kendallsTauE1, 1, mean)
print(mean(abs(kendallsTauE1Estimate - trueKendallsTauE1))) # 0.09344631


