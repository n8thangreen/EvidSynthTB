
# Evidence synthesis LTBI screening: Stan
# with real data


library(rstan)
library(shinystan)
library(purrr)
library(readr)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


edetecttb <- read_csv("../../data/clean_edetecttb_11oct_withregionMAR.csv")
load("data input/cleaned_migrant_predict_data.RData")

predicttb <- filter(dat_m, pos == TRUE)

dat_input <-
  list(
    N = nrow(predicttb),
    t = as.numeric(predicttb$time),
    d = as.numeric(predicttb$status),
    mu_lambda = -4,
    sigma_lambda = 1,
    mu_gamma = -1,
    sigma_gamma = 1,
    t_lim = 20)

params <- c("ppred",
            "lambda",
            "gamma")

n_iter <- 2e3
n_burnin <- 1e1
n_thin <- 1e1 #floor((n_iter - n_burnin)/500)


##############
## run MCMC ##
##############

out <- stan(data = dat_input,
            pars = params,
            file = here::here("BUGS", "Stan_code_predicttb.stan"),
            chains = 2,
            iter = n_iter,
            warmup = n_burnin,
            thin = n_thin,
            control = list(adapt_delta = 0.9,
                           max_treedepth = 20))

stan_output <- extract(out)

save(stan_output, file = here::here("data output", "stan_output.RData"))

# flexsurv:
# shape       rate
# 0.42823033 0.02314312

mean(stan_output$lambda)
mean(stan_output$gamma)


###########
# plot

library(reshape2)
library(ggplot2)

# base R
stan_output$ppred %>%
  as.data.frame %>%
  summarise_all(mean) %>%
  unlist() %>%
  plot(ylim = c(0.8, 1), type = "l")


plot_dat <-
  stan_output$ppred %>%
  as.data.frame() %>%
  mutate(sim = 1:n()) %>%
  melt(id.vars = "sim",
       variable.name = "time") %>%
  mutate(time = as.numeric(gsub("V", "", time))) %>%
  group_by(time) %>%
  summarise(mean = mean(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975))

ggplot(plot_dat, aes(time, mean)) +
  geom_line() +
  ylab("Survival") +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
              linetype = 0,
              alpha = 0.2) +
  ylim(0, 1)

