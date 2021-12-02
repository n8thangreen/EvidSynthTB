
library(brms)

fit_brms <- brm(formula = pos ~ age*ethnicity,  
                data = pos_dat, 
                family = bernoulli(link = "logit"),
                warmup = 500, 
                iter = 2000, 
                chains = 2, 
                inits = "0", 
                cores = 2,
                seed = 123)

save(fit_brms, file = "data output/rbms.RData")

stan_output <- brms::as_draws_df(fit_brms)

summary(fit_brms)

stanplot(fit_brms,
         type = "areas",
         prob = 0.95) +
  xlim(-1, 1)


# frequentist
fit <-
  dat_m %>%
  glm(pos ~ age*ethnicity,
      data = .,
      family = "binomial")
