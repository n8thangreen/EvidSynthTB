# make plots using jags output
# LTBI evidence synthesis model


library(MCMCvis)
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

MCMCtrace(out,
          pdf = FALSE,
          ind = TRUE)


n_sample <- (n_iter-n_burnin)/n_thin

## with priors ----

sens <- rbeta(n = n_sample, shape1 = 100, shape2 = 5)
spec <- rbeta(n = n_sample, shape1 = 100, shape2 = 5)
lambda <- rbeta(n = n_sample, shape1 = 5, shape2 = 100)
logitp_latent <- rnorm(n = n_sample, mean = 0, sd = 1/0.001)
p_latent <- exp(logitp_latent) / (1 + exp(logitp_latent))
p_latent[is.nan(p_latent)] <- 0 #remove NaN

MCMCtrace(out,
          pdf = FALSE,
          ind = TRUE,
          iter = n_sample,
          params = c("lambda", "sens", "spec", "p_latent"),
          priors = as.matrix(data.frame(lambda, sens, spec, p_latent, p_latent)))

# forest plot
MCMCplot(out,
         params = c("sens", "spec", "lambda", "p_latent"))

## two model plots
# MCMCplot(object = MCMC_data,
#          object2 = MCMC_data2,
#          params = c("sens", "spec", "lambda"),
#          offset = 0.1)


## BGR statistics
coda::gelman.plot(mcmcplots::as.mcmc.rjags(out))
