Software Documentation for “A Copula-based Fully Bayesian Nonparametric
Evaluation of Cardiovascular Risk Markers for normoglycemic patients in
the Mexico City Diabetes Study”

In what follows we give a brief description of the files and folders 
contained in this Software Documentation together with instructions
for reproducibility of the results for Simulation Scenario I.

#———————————————————————————————————————————————————————————————————————

Folders:

1. Data: where simulated data set is saved. We already included the 
    simulated data set for Scenario I.
2. Results: where the outputs and graphics are saved.
3. DataApplication: contains data from the MCDS analysis.

#———————————————————————————————————————————————————————————————————————

Files:

1. Functions.R: contains several functions used when running the Markov 
    chain Monte Carlo (MCMC) algorithm and when computing Kendall’s tau.
2. runCodeSimulation.R: runs the MCMC sampler for simulation scenario I
    using the proposed model,  computes L1 distance, and reproduces plot  
    reported in the paper for simulation Scenario I.
3. runCodeTest.R: runs a very short MCMC sampler for simulation Scenario I
    using the proposed model. 
4. simulatingDataSets.R: computes true quantities and (pseudo) randomly 
    generates the data set for simulation Scenario I.

#———————————————————————————————————————————————————————————————————————

Reproducibility:

1. Execute the runCodeTest.R file to verify that everything is in order 
    before running the long MCMC.
2. Execute runCodeSimulation.R file.
3. To reproduce the simulated data set execute the simulatingDataSets.R 
    file.

#———————————————————————————————————————————————————————————————————————

Timings:

- running the mcmc takes ~10 hours.
- postprocessing take about ~14 hours using 20 cores.
