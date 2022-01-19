# Part 2 initial exploration
library(dplyr)
library(ggplot2)
library(GGally) 
library(nimble)

############### Load and preprocess data ######################
AIdata = readRDS('AIdataset.Rds') 
# Data where icustm plays a role
sum(AIdata$day>AIdata$icustaymv)
# As only a few specific cases due to censoring we remove it
AIdata = AIdata %>% select(-c(icustaymv))
# Make fail into factor
AIdata$fail = as.factor(AIdata$fail)

# Add small amount to di to avoid zero values which cause errors
AIdata$di = AIdata$di+0.0001

# Plot pair and relation between data
ggpairs(AIdata %>% select(-c(id)),   
        columns = 1:4,     
        aes(color = fail,  
            alpha = 0.5))    


################## Analysis 1 analyse sofa score ##################

# Overview analysis
ggplot(AIdata, aes(x=day, y=age, color=sofa)) +
  geom_point(size=5)+
  geom_smooth(method=lm,formula = y~x,se = TRUE)+
  theme_classic()

# Pairplot
ggpairs(AIdata %>% select(c(age,day,sofa)) %>% filter(!is.na(sofa)),   
        columns = 1:3)
# Should motivate choice of smoothing technique

##### Train linear mixed model with age and day covariates
# This can be improved by smoothing variables and other type of regression
sof_model <- nimbleCode({
  for (i in 1:N) {
    mu <- beta0+inprod(beta[1:nbetas],x[i,1:nbetas])
    y[i] ~ dnorm(mu,sigma)
  }
  beta0 ~ dnorm(0,0.001)
  for(i in 1:nbetas){
    beta[i]~dnorm(0,0.001)
  }
  sigma ~ dunif(0,100)
})

# Make small adjusted data smoothing with log
X <- AIdata %>% select(-c(id,di,fail)) %>% filter(!is.na(sofa))
X$day = log(X$day)
sof_data <- list(
  x = data.matrix(X %>% select(-c(sofa))),
  y = X$sofa
)

## constants, data, and initial values
sof_const<- list(N = nrow(X),nbetas=2)
sof_inits <- list(beta0=0,beta =rep(0,2), sigma = 1)

## create the model object
sofModel <- nimbleModel(code = sof_model, constants = sof_const, data = sof_data, 
                         inits = sof_inits, check = FALSE)

# Use default MCMC algorithm
sof_glmmMCMC <- buildMCMC(sofModel)

# Compile MCMC algorithm
C_sof_glmmModel <- compileNimble(sofModel)
C_sof_glmmMCMC <- compileNimble(sof_glmmMCMC, project = sofModel)

# Run MCMC and sample
burnin = 10000
n_samples = 1000
n_chains=1
samples <- runMCMC(C_sof_glmmMCMC, niter = burnin+n_samples,nchains=n_chains)

outpt = samples[burnin:(n_samples+burnin),]
par(mfrow = c(2, 2))
coda::traceplot(coda::as.mcmc(outpt),col=1:4)

# Extract and predict sofa
betas = colMeans(outpt)[1:2]
beta0 = colMeans(outpt)[3]
sig = sqrt(1/colMeans(outpt)[4])

# Predict actual values
x = AIdata %>% select(age,day)
x$day = log(x$day)
x = data.matrix(x)
mu = beta0 +x %*% betas[2:1]
sofa_pred = mu


# Add into AIdata
AIdata$sofa_pred = sofa_pred 
AIdata$sofa_adj = AIdata$sofa-sofa_pred 

# Adjust NA values to be prediction
AIdata[is.na(AIdata$sofa_adj),]$sofa_adj= AIdata[is.na(AIdata$sofa_adj),]$sofa_pred

AIdata
############### Analysis 2 ####################################

############ Plotting and explroation analysis
ggplot(AIdata, aes(x=day, y=di, color=age)) +
  geom_point(size=5)+
  geom_smooth(method=lm,formula = y~x,se = TRUE)+
  theme_classic()

# Explore data
ggpairs(AIdata %>% select(c(di,age,day)),columns = 1:3)

########### Model
di_model <- nimbleCode({
  for (i in 1:N) {
    eta[i] <- inprod(betas[1:nbetas],x[i,1:nbetas]) 
    mu[i] <- exp(eta[i])/(1+exp(eta[i]))
    alpha[i] <- phi*mu[i]
    beta[i] <- phi*(1-mu[i])
    y[i] ~ dbeta(alpha[i],beta[i])
  }
  for(i in 1:nbetas){
    betas[i]~dnorm(0,0.001)
  }
  phi~dnorm(0,0.001)
})

X <- AIdata %>% select(c(day,age,di))
X$intercept = 1
di_data <- list(
  x = data.matrix(X %>% select(-c(di))),
  y = X$di
)

n = nrow(data$x)
nbetas=3

## constants, data, and initial values
di_const <- list(N = n,nbetas=nbetas)

# Calculate initial values
# Use regular regression for beta intiation
x = di_data$x
z = exp(logit(AIdata$di))
betas =  solve(t(x) %*% x) %*% (t(x) %*% z)

# Use suggested intiial esstimate for phi
y_pred = x %*% betas
muhat = exp(y_pred)/(1+exp(y_pred))
ehat = z-y_pred
sighat  = (t(ehat) %*% ehat / (nrow(x)-nbetas))[1] * (1/(muhat*(1-muhat)))^2
phi_init = mean(muhat*(1-muhat)/sighat)-1

# Insert as inits
di_inits <- list(betas =c(betas), phi = phi_init)

## create the model object
di_glmmModel <- nimbleModel(code = di_model, constants = di_const, 
                            data = di_data, inits = di_inits, check = FALSE)

# Use default MCMC algorithm
di_glmmMCMC <- buildMCMC(di_glmmModel)

# Compile MCMC algorithm
di_CglmmModel <- compileNimble(di_glmmModel)
di_CglmmMCMC <- compileNimble(glmmMCMC, project = di_glmmModel)

# Run MCMC and sample
burnin = 10000
n_samples = 1000
n_chains=1
samples <- runMCMC(di_CglmmMCMC, niter = burnin+n_samples,nchains=n_chains)

# Analyse output
outpt = samples[burnin:(n_samples+burnin),]
coda::traceplot(coda::as.mcmc(outpt),col=1:4)

# Extract and predict sofa
betas = colMeans(outpt)[1:2]
beta0 = colMeans(outpt)[3]
phi = colMeans(outpt)[4]

# Predict actual values
x = AIdata %>% select(age,day)
x$day = log(x$day)
x = data.matrix(x)
mu = beta0 +x %*% betas
eta = exp(mu)/(1+exp(mu))
di_pred = eta


# Add into AIdata
AIdata$di_pred = di_pred 
AIdata$di_adj = AIdata$di-di_pred 

AIdata


################# Analysis 3 ##############################

# Adjust fails
AIdata$fail_adj = AIdata$fail
AIdata[AIdata$fail_adj==0,]$fail_adj = 1
AIdata[AIdata$fail_adj==2,]$fail_adj = 0

# Set up model
weiProHazModel <- nimbleCode({
  for (i in 1:N) {
    # Linear predictor
    mu[i] <- inprod(beta[1:nbetas],x[i,1:nbetas])
    # Baseline Hazard function
    h0 <- lambda^(-k)*k*t[i]^(k-1)
    # Hazard esstimate
    y[i] <- h0*exp(mu[i])
  }
  for(i in 1:nbetas){
    beta[i]~dnorm(0,0.001)
  }
  lambda ~ dgamma(1, .0001)
  k ~ dgamma(1, .0001)
})

X <- AIdata %>% select(c(fail_adj,di_pred,sofa_pred,day))
wei_data <- list(
  x = data.matrix(X %>% select(-c(fail_adj,day))),
  y = as.numeric(X$fail_adj),
  t = X$day
)

## constants, data, and initial values
wei_const <- list(N = nrow(X),nbetas=2)
wei_inits <- list(beta =rep(0,2), lambda=1,k=1)

## create the model object
wei_glmmModel <- nimbleModel(code = weiProHazModel, constants = wei_const, 
                         data = wei_data, inits = wei_inits, check = FALSE)

# Use default MCMC algorithm
wei_glmmMCMC <- buildMCMC(wei_glmmModel)

# Compile MCMC algorithm
wei_CglmmModel <- compileNimble(wei_glmmModel)
wei_CglmmMCMC <- compileNimble(wei_glmmMCMC, project = wei_glmmModel)

# Run MCM and 
burnin = 100000
n_samples = 1000
n_chains=1
samples <- runMCMC(wei_CglmmMCMC, niter = burnin+n_samples,nchains=n_chains)

outpt = samples[burnin:(n_samples+burnin),]
coda::traceplot(coda::as.mcmc(outpt),col=1:4)

# Extract and predict sofa
betas = colMeans(outpt)[1:2]
k = colMeans(outpt)[3]
lambda = colMeans(outpt)[4]

betas

################# Analysis 4 #################################
bls_Model <- nimbleCode({
  for (i in 1:N) {
    eps[i] ~ dnorm(0,sd=sigma)
    logit(p[i]) <- inprod(beta[1:nbetas],x[i,1:nbetas])+eps[i]
    y[i] ~ dbern(p[i])
  }
  for(i in 1:nbetas){
    beta[i]~dnorm(0,0.001)
  }
  sigma ~ dunif(0,1000)
  })

X <- AIdata %>% select(di_adj,sofa_adj,day,age,fail_adj)
X$intercept = 1
bls_data <- list(
  x = data.matrix(X %>% select(-c(fail_adj))),
  y = as.integer(X$fail==2)
)

## constants, data, and initial values
bls_const <- list(N = nrow(X),nbetas=5)
bls_inits <- list(beta =rep(0,5), sigma = 1)

## create the model object
bls_glmmModel <- nimbleModel(code = bls_Model, constants = bls_const, 
                             data = bls_data,inits = bls_inits, check = FALSE)

# Use default MCMC algorithm
bls_glmmMCMC <- buildMCMC(bls_glmmModel)

# Compile MCMC algorithm
bls_CglmmModel <- compileNimble(bls_glmmModel)
bls_CglmmMCMC <- compileNimble(bls_glmmMCMC, project = bls_glmmModel)

# Run MCM and 
burnin = 100000
n_samples = 1000
n_chains=1
samples <- runMCMC(bls_CglmmMCMC, niter = burnin+n_samples,nchains=n_chains)

outpt = samples[1:(n_samples+burnin),]
coda::traceplot(coda::as.mcmc(outpt),col=1:4)

# Extract and predict sofa
betas = colMeans(outpt)[1:5]
sigma = colMeans(outpt)[6]

betas
sigma
