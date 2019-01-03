set.seed(32)

# data plotting
ids <- unique(NIC$id)
plot(NIC$age[NIC$id==1], NIC$bw[NIC$id==1], "l")
for (idi in NIC$id) lines(NIC$age[NIC$id==idi], NIC$bw[NIC$id==idi]);

# model for questions 1 and 2
cat("
    model {
    for (i in 1:N) {
    body_weight[i] ~ dnorm(body.mu[i], body.tau)
    body.mu[i] <- beta.const + beta.breed*breed[i] + beta.age*age[i] + b[id[i],1] + b[id[i],2]*age[i]
    }
    
    for (i in 1:M) {
    b[i, 1:2] ~ dmnorm(b_mean[1:2], b.tau[1:2,1:2])
    }
    
    # error SD
    body.tau ~ dgamma(0.00001,0.00001)
    body.sigma <- 1 / sqrt(body.tau)
    
    # fixed effects
    beta.const ~ dnorm(0, 1.0E-6)
    beta.breed ~ dnorm(0, 1.0E-6)
    beta.age ~ dnorm(0, 1.0E-6)
    
    # random effects
    b.sigma[1] ~ dunif(0,10)
    b.sigma[2] ~ dunif(0,10)
    b.corr ~ dunif(-1,1)
    
    for (i in 1:2) {
    b.sigma2[i,i] <- pow(b.sigma[i], 2)
    }
    b.sigma2[1,2] <- b.sigma[1] * b.sigma[2] * b.corr
    b.sigma2[2,1] <- b.sigma[1] * b.sigma[2] * b.corr
    
    b.tau <- inverse(b.sigma2[1:2,1:2])
    
    error <- mean(body_weight - body.mu)
    b.intercept.means <- mean(b[,1])
    b.slope.means <- mean(b[,2])
    
    }", file="./project1/code/jags-programs/NIC-norm-uni-model.jag")

NIC_data <- list(
  body_weight=normalize(NIC$bw),
  breed=normalize(NIC$breed),
  age=normalize(NIC$age),
  id=NIC$id,
  N=nrow(NIC),
  b_mean=c(0,0),
  b_R=diag(0.1,2),
  M=length(unique(NIC$id))
)

NIC_inits <- list(
  list(
    body.tau=25,
    beta.const=-0.07,
    beta.breed=-0.05,
    beta.age=0.8,
    b.tau=matrix(c(32,0,0,58),2,2)
  ),
  list(
    body.tau=25,
    beta.const=-0.02,
    beta.breed=-0.051,
    beta.age=0.82,
    b.tau=matrix(c(32,0,0,58),2,2)
  ),
  list(
    body.tau=24,
    beta.const=-0.1,
    beta.breed=-0.41,
    beta.age=0.84,
    b.tau=matrix(c(32,0,0,58),2,2)
  )
)

# model fitting
NIC_model_question1 <- jags(
  NIC_data,
  parameters.to.save=c("beta.const", "beta.breed", "beta.age", "body.sigma", "b.sigma", "b.corr"),
  inits=NIC_inits,
  n.iter=50000,
  n.chains=3,
  n.burnin=25000,
  model.file="./project1/code/jags-programs/NIC-norm-wish-model.jag",
  n.thin=1,
  DIC=T
)

# convergence checks
NIC_model_MCMC <- as.mcmc(NIC_model_question1)
NIC_model_gg <- ggs(NIC_model_MCMC)
print(NIC_model_question1)
summary(NIC_model_MCMC)
ggmcmc(NIC_model_gg, "./project1/question1_plots-norm-wish-distribution.pdf")
superdiag(NIC_model_MCMC, burnin = 1000)


# question 2

# model fitting
NIC_model_question2 <- jags(
  NIC_data,
  parameters.to.save=c("error", "b.intercept.means", "b.slope.means"),
  inits=NIC_inits,
  n.iter=50000,
  n.chains=3,
  n.burnin=25000,
  model.file="./project1/code/jags-programs/NIC-norm-wish-model.jag",
  n.thin=1,
  DIC=T
)

NIC_model_MCMC <- as.mcmc(NIC_model_question2)
NIC_model_gg <- ggs(NIC_model_MCMC)

# histograms
ggs_histogram(NIC_model_gg)

# Q-Q plots
data <- NIC_model_gg
error <- data$value[data$Parameter == "error"]
b.intercept.means <- data$value[data$Parameter == "b.intercept.means"]
b.slope.means <- data$value[data$Parameter == "b.slope.means"]

qqnorm(error, col="darkgreen", main="Normal Q-Q Plot of Measurement Error")
qqline(error)

qqnorm(b.intercept.means, col="darkgreen", main="Normal Q-Q Plot of Random Intercept Means")
qqline(b.intercept.means)

qqnorm(b.slope.means, col="darkgreen", main="Normal Q-Q Plot of Random Slope Means")
qqline(b.slope.means)


