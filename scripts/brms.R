
library(brms)
library(ggplot2)

load("data input/cleaned_migrant_predict_data.RData")

pos_dat <-
  select(dat_m, pos, age, ethnicity) %>%
  mutate(age = as.numeric(age)) %>%
  na.omit()

fit_brms <- brm(formula = pos ~ age*ethnicity,
                data = pos_dat,
                family = bernoulli(link = "logit"),
                warmup = 500,
                iter = 2000,
                thin = 10,
                chains = 2,
                inits = "0",
                cores = 3,
                seed = 123)

save(fit_brms, file = "data output/brms_pltbi.RData")

Xnew <-
  expand.grid(age = 1:89,
              ethnicity = levels(pos_dat$ethnicity))

pred_brms <-
  fitted(fit_brms, newdata = Xnew)

p_ltbi_brms <- cbind(Xnew, prob_ltbi = pred_brms)

ggplot(p_ltbi_brms, aes(x = age, y = prob_ltbi.Estimate, group = ethnicity, col = ethnicity)) +
  geom_line()

save(p_ltbi_brms, file = "data output/p_ltbi_brms.RData")

## frequentist

fit <-
  dat_m %>%
  glm(pos ~ age*ethnicity,
      data = .,
      family = "binomial")

pred <-
  predict(fit, newdata = Xnew, type = "response")

prob_ltbi <- cbind(Xnew, prob_ltbi = pred)

ggplot(prob_ltbi, aes(x = age, y = prob_ltbi, group = ethnicity, col = ethnicity)) +
  geom_line()

