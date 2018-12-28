set.seed(32)

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
      body.tau ~ dgamma(0.001,0.001)
      body.sigma <- 1 / sqrt(body.tau)

      # fixed effects
      beta.const ~ dnorm(0, 1.0E-3)
      beta.breed ~ dnorm(0, 1.0E-3)
      beta.age ~ dnorm(0, 1.0E-3)

      # random effects
      b.tau ~ dwish(b_R, 2)
      b.sigma2[1:2,1:2]  <- inverse(b.tau[1:2,1:2])
      for (i in 1:2) {
        b.sigma[i] <- sqrt(b.sigma2[i,i])
      }
      b.corr <- b.sigma2[1,2]/(b.sigma[1]*b.sigma[2])
      
  
    }", file="./project1/code/jags-programs/NIC-model1.jag")

NIC_data <- list(
  body_weight=NIC$bw,
  breed=NIC$breed,
  age=NIC$age,
  id=NIC$id,
  N=nrow(NIC),
  b_mean=c(0,0),
  b_R=diag(0.1,2),
  M=length(unique(NIC$id))
)

# NIC_inits <- list()

NIC_model <- jags(
  NIC_data,
  parameters.to.save=c("beta.const", "beta.breed", "beta.age", "body.sigma", "b.sigma", "b.corr"),
  n.iter=10000,
  n.chains=3,
  n.burnin=5000,
  model.file="./project1/code/jags-programs/NIC-model1.jag",
  n.thin=1,
  DIC=T
)

NIC_model_MCMC <- as.mcmc(NIC_model)
print(NIC_model)
plot(NIC_model)
traceplot(NIC_model)
superdiag(NIC_model_MCMC)



