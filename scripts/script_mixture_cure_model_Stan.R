
# Evidence synthesis LTBI screening: Stan
# with real data


library(rstan)
library(shinystan)
library(purrr)
library(readr)
library(dplyr)

# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


load("data input/cleaned_migrant_predict_data.RData")

dat_m <- dat_m[!is.na(dat_m$age), ]
dat_m <- dat_m[dat_m$ethnicity != "", ]

dat_input <-
  list(
    N = nrow(dat_m),
    t = as.numeric(dat_m$time),
    d = as.numeric(dat_m$status),
    pos = dat_m$pos,
    age = as.numeric(dat_m$age) - mean(dat_m$age, na.rm = TRUE),
    ethnicity =  as.numeric(droplevels(dat_m$ethnicity)),
    n_eth = length(unique(dat_m$ethnicity)),
    t_lim = 20,
    scale_alpha = 10,
    scale_beta = 10,
    a_eth = 0.01,
    b_eth = 0.01,
    ## normal priors
    mu_shape = -0.8481,
    sigma_shape = 1,
    ## gamma priors
    a_lambda = 0.01, # rate
    b_lambda = 0.01)

params <- c(
  "S_pred",
  "cf", "alpha",
  "beta_age", "beta_eth",
  "lambda", "shape")

n_iter <- 10e3
n_burnin <- 3e1
n_thin <- 2e1 #floor((n_iter - n_burnin)/500)

##TODO: try with artificial clean data?


###########
# run MCMC
###########

out <- stan(data = dat_input,
            pars = params,
            file = here::here("BUGS", "Stan_code_mixture_cure_model.stan"),
            chains = 1,
            iter = n_iter,
            warmup = n_burnin,
            thin = n_thin,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 20))

stan_output <- extract(out)

save(stan_output, file = here::here("data output", "stan_output_mcm.RData"))

mean(stan_output$lambda)
mean(stan_output$shape)
mean(stan_output$cf)


#######
# plot
#######

library(reshape2)
library(ggplot2)

# base R
stan_output$S_pred %>%
  as.data.frame %>%
  summarise(across(.fns = ~ mean(.x, na.rm = TRUE ))) %>%
  # summarise(across(.fns = ~ median(.x, na.rm = TRUE ))) %>%
  unlist() %>%
  plot(ylim = c(0.5, 1), type = "l")


plot_dat <-
  stan_output$S_pred %>%
  as.data.frame() %>%
  mutate(sim = 1:n()) %>%
  melt(id.vars = "sim",
       variable.name = "time") %>%
  mutate(time = as.numeric(gsub("V", "", time))) %>%
  group_by(time) %>%
  summarise(median = median(value),
            mean = mean(value),
            lower50 = quantile(value, probs = 0.25),
            upper50 = quantile(value, probs = 0.75),
            lower95 = quantile(value, probs = 0.025),
            upper95 = quantile(value, probs = 0.975))

ggplot(plot_dat, aes(time, median)) +
  geom_line() +
  geom_line(aes(y = mean), linetype = 2) +
  ylab("Survival") +
  geom_ribbon(aes(x = time, ymin = lower95, ymax = upper95),
              linetype = 0,
              alpha = 0.2) +
  geom_ribbon(aes(x = time, ymin = lower50, ymax = upper50),
              linetype = 0,
              alpha = 0.2) +
  ylim(0, 1)


# LTBI prevalence

for (i in 1:ncol(stan_output$beta_eth)) {

  lrg <- stan_output$alpha + stan_output$beta_eth[, i]
  print(summary(1/(1 + exp(-lrg))))
  hist(1/(1 + exp(-lrg)), breaks = 40)
}


