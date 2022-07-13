
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

dat_m <- dat_m[!is.na(dat_m$age), ]

dat_input <-
  list(
    N = nrow(dat_m),
    t = as.numeric(dat_m$time),
    d = as.numeric(dat_m$status),
    K = 2,
    age = as.numeric(dat_m$age) - mean(dat_m$age, na.rm = TRUE),
    ethnicity =  as.numeric(dat_m$ethnicity),
    t_lim = 20,
    scale_alpha = 1,
    scale_beta = 1,
    ## normal priors
    mu_shape = -0.8481,
    sigma_shape = 0.1,
    ## gamma priors
    a_lambda = 0.1, # rate
    b_lambda = 1)

params <- c(
  # "ppred",
  "cf", "alpha",
  "lambda",
  "shape")

n_iter <- 10e3
n_burnin <- 3e1
n_thin <- 2e1 #floor((n_iter - n_burnin)/500)


###########
# run MCMC
###########

out <- stan(data = dat_input,
            pars = params,
            file = here::here("BUGS", "Stan_code_predicttb.stan"),
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


#######
# plot
#######

library(reshape2)
library(ggplot2)

# # base R
# stan_output$ppred %>%
#   as.data.frame %>%
#   summarise_all(mean) %>%
#   unlist() %>%
#   plot(ylim = c(0.8, 1), type = "l")
#
#
# plot_dat <-
#   stan_output$ppred %>%
#   as.data.frame() %>%
#   mutate(sim = 1:n()) %>%
#   melt(id.vars = "sim",
#        variable.name = "time") %>%
#   mutate(time = as.numeric(gsub("V", "", time))) %>%
#   group_by(time) %>%
#   summarise(mean = mean(value),
#             lower = quantile(value, probs = 0.025),
#             upper = quantile(value, probs = 0.975))
#
# ggplot(plot_dat, aes(time, mean)) +
#   geom_line() +
#   ylab("Survival") +
#   geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
#               linetype = 0,
#               alpha = 0.2) +
#   ylim(0, 1)


