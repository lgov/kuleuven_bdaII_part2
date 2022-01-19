# Code for the second project in BDA II.
# For nimble, see resources: https://github.com/nimble-training/AHMnimble
# and the book 
# Applied Hierarchical Modeling in Ecology: Analysis of distribution, abundance and species richness in R and BUGS
library(nimble)
library(coda)

setwd("D:/asus_documents/ku_leuven/courses/bayesian_II/project2/")
set.seed(13)

# Variable name Description
# id Patient id
# day Time (in days) in the Intensive Care Unit since the study entry
# icustaymv Number of days in mechanical ventilation
# fail Type of failure (1: dead in the ICU; 2: alive discharge; 0: censored at
#                       day 30)
# sofa Sequential Organ Failure Assessment (score with range 0-24)
# di Asynchrony index is defined as the proportion of asynchronous events
# among the total number of ventilator cycles (range 0-1)
# age Age (yrs)
ai <- readRDS("AIdataset.Rds")
nburnin <- 1000
niter <- 10000
nchains <- 3

# First fit a linear mixed model 
# We also have to account for missing data in the sofa variable.
lmm_code <- nimbleCode({
    for(i in 1:N) {
        y[i] ~ dnorm(mu_y[i], tau_y)
        mu_y[i] <- b_unidentifiable[id[i]] + betas[1] + betas[2]*age[i] + betas[3]*day[i]
    }
    # Set up the priors.
    # For y params
    tau_y <- pow(sd_y, -2)
    sd_y ~ dunif(0, 10)
    # For main effects
    for(k in 1:3){
        betas[k] ~ dnorm(0, 0.001)
    }
    # For the random effects
    # Post sweeping
    for(j in 1:J){
        b_unidentifiable[j] ~ dnorm(mu_b, tau_b)
        # Monitor this set of random effects
        b[j] <- b_unidentifiable[j] - avg_b
    }
    # Intercept to monitor
    beta_star <- betas[1] + avg_b
    avg_b <- mean(b_unidentifiable[])
    tau_b <- pow(sd_b, -2)
    sd_b ~ dunif(0, 10)
    mu_b ~ dnorm(0, .0001)
    
})
params <- c("beta_star", "betas", "b") 
N <- nrow(ai)
J <- length(unique(ai$id))
# Apparently using id as a constant is more computationally efficient.
constants <- list(N=N, J=J, id=ai$id)
data <- list(y=ai$sofa, age=ai$age-mean(ai$age), day=ai$day)
inits <- list(betas=c(0,0,0))
dimensions <- list(b=J, b_unidentifiable=J)
lmm <- nimbleMCMC(code = lmm_code, constants = constants, inits = inits, 
                  data = data, nburnin=nburnin, niter=niter, nchains=nchains, monitors=params, dimensions=dimensions,
                  summary = TRUE, WAIC = TRUE)

lmm$summary
lmm$WAIC
# If summary and WAIC are used, then need to do lmm$samples$chain1
# otherwise do $lmm$chain1
plot(lmm$samples$chain1[,1], type="l")
# intercept
# ignore 141
plot(lmm$samples$chain1[,140], type="l")
# age
plot(lmm$samples$chain1[,142], type="l")
# day
plot(lmm$samples$chain1[,143], type="l")
# Random effects seem to ave converged.
plot(lmm$samples$chain1[,1], type="l")

