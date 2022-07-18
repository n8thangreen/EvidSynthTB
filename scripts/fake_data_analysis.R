
# simulates artificial ltbi and tb data

library(rstan)
library(shinystan)
library(purrr)
library(readr)
library(dplyr)


N <- 100

# covariates
ids <- sample(1:3, size = N, replace = TRUE)

rdat <-
  data.frame(
    age =  rep(50, N),
    ethnicity = c("White", "China", "Other")[ids],
    t = round(flexsurv::rgompertz(N, shape = 1, rate = 0.1), 3),
    probs = c(0.1, 0.5, 0.8)[ids]) |>
    mutate(x = rbinom(n = N, size = 1, prob = probs),
           d = sample(c(0,1), size = N, replace = TRUE),
           d = ifelse(x == 0, 0, d))

rdat

# for (i in prob) {
#   x[i] <- purrr::rbernoulli(1, p = i)
# }

# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


dat_input <-
  list(
    N = nrow(rdat),
    t = as.numeric(rdat$t),
    d = as.numeric(rdat$d),
    pos = rdat$x,
    age = as.numeric(rdat$age) - mean(rdat$age, na.rm = TRUE),
    ethnicity =  as.numeric(as.factor(rdat$ethnicity)),
    n_eth = length(unique(rdat$ethnicity)),
    t_lim = 20,
    scale_alpha = 10,
    scale_beta = 10,
    a_eth = 0.01,
    b_eth = 0.01,
    ## normal priors
    mu_shape = 1,
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

mean(stan_output$lambda)
mean(stan_output$shape)


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


