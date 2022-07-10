
# combine tb data sets
#
# want to take progression rate from predict-tb
# and denominators from e-detect
# to then estimate number of expected cases by
# groups for each year after entry
# then we can compare with ETS
#
# ETS/migrant journey: migrant numbers by (age), ethnicity
# PREDICT: prevalence, progression rate
#          by age, contact/migrant, ethnicity
#
#
# 5 year age groups
# start at different years for each country because PES different
#   Bangladeshi: 2011
#   Black African: 2014
#   Indian: 2015
#   Pakistani: 2011


##TODO: need to filter each ethnic group by
##      range of year; different for different ethnicity
##      matched for ETS and migrant journey

library(readr)
library(dplyr)
library(stringr)
library(reshape2)

ETS <- read_csv("../../data/ETS_for_PREDICT_comparison_new.csv")

load("../migrantJourney/data/migrant_pop.RData")  # entry populations
load("../migrantJourney/data/posterior_S.RData")  # survival posteriors

migrant_S <- melt(res)

# predict-tb model fit stats
load("~/R/EvidSynthTB/data output/brms_pltbi.RData")   # Stan output
load("~/R/EvidSynthTB/data output/p_ltbi_brms.RData")
# active tb survival
load("~/R/EvidSynthTB/data output/stan_output.RData")  # gompertz

# aggregate ethnicity into 5 yr age groups
pltbi_agegrp <-
  p_ltbi_brms %>%
  mutate(age_grp = cut(age, breaks = seq(0, 95, 5), right = FALSE)) %>%
  group_by(age_grp, ethnicity) %>%
  summarise(mean_pltbi = mean(prob_ltbi.Estimate),
            low_pltbi = mean(prob_ltbi.Q2.5),
            upp_pltbi = mean(prob_ltbi.Q97.5))  ##TODO: how to best aggregate this statistic?

# combined ages
pltbi_ethgrp <-
  p_ltbi_brms %>%
  group_by(ethnicity) %>%
  summarise(mean_pltbi = mean(prob_ltbi.Estimate),
            low_pltbi = mean(prob_ltbi.Q2.5),
            upp_pltbi = mean(prob_ltbi.Q97.5))

## clean migrant journey data

#
entry_pop <-
  migrant_pop |>
  mutate(nationality = as.character(nationality)) |>
  mutate(nationality = replace(nationality, nationality == "Sub_Saharan_Africa", "Black African"),
         nationality = replace(nationality, nationality == "Bangladesh", "Bangladeshi"),
         nationality = replace(nationality, nationality == "India", "Indian"),
         nationality = replace(nationality, nationality == "Pakistan", "Pakistani")) %>%
  rename(ethnicity = nationality,
         year_issue = year_issued) |>
  group_by(year_issue, ethnicity) |>
  summarise(pop = sum(pop))

# combine migrant counts and S probabilities
# combine estimate number ltbi
dat <-
  migrant_S %>%
  mutate(nat = as.character(nat)) |>
  mutate(nat = replace(nat, nat == "Sub_Saharan_Africa", "Black African"),
         nat = replace(nat, nat == "Bangladesh", "Bangladeshi"),
         nat = replace(nat, nat == "India", "Indian"),
         nat = replace(nat, nat == "Pakistan", "Pakistani")) |>
  rename(ethnicity = nat,
         S_expire = value) |>
  inner_join(entry_pop) |>
  mutate(year_pop = S_expire*pop) |>
  left_join(pltbi_ethgrp, by = c("ethnicity")) |>
  # left_join(pltbi_agegrp, by = c("age_grp", "ethnicity")) |>
  mutate(n_ltbi = year_pop*mean_pltbi,
         n_ltbi_low = year_pop*low_pltbi,
         n_ltbi_upp = year_pop*upp_pltbi) |>
  filter(year_issue <= year_expire) |>
  mutate(time = year_expire - year_issue + 1)
  # full_join(data.frame(time = 1:10), by = character()) |>
  # na.omit()

##TODO: how to do this for posterior samples?
##      assume p tb is small so point estimate?

## tb progression
# fitted values
ppred <- stan_output$ppred

S_m <-
  data.frame(time = 0:ncol(ppred),
             S_tb = c(1, colMeans(ppred))) %>%
  mutate(p_tb = lag(S_tb) - S_tb)

dat_with_tb <-
  left_join(dat, S_m, by = "time") %>%
  mutate(n_tb = round(n_ltbi*p_tb, 2),
         n_tb_low = round(n_ltbi_low*p_tb, 2),
         n_tb_upp = round(n_ltbi_upp*p_tb, 2))

save(dat_with_tb, file = "data output/dat_with_tb.RData")
write.csv(dat_with_tb, file = "data output/dat_with_tb.csv")


###################
# compare with ETS
###################

library(readr)
library(reshape2)

ETS <- read_csv("../../data/ETS_for_PREDICT_comparison_new.csv")

# clean ETS data
# to long format
ETSm <-
  melt(ETS,
       id.vars = c("age_grp_at_entry", "ethgrp", "entry_year"),
       variable.name = "time",
       value.name = "n_tb_ETS") %>%
  mutate(time = str_remove(str_remove(time, "(\\[|\\().,"), "\\]")) %>%
  rename(age_grp = "age_grp_at_entry",
         year = "entry_year",
         ethnicity = "ethgrp") %>%
  filter(year != "< 2000") %>%
  mutate(year = as.numeric(year),
         time = as.numeric(time),
         ethnicity = replace(ethnicity, ethnicity == "Black-African", "Black African"),
         ethnicity = replace(ethnicity, ethnicity == "White", "Europe")) %>%
  ## filter each ethnicity by different calendar years
  # filter(year %in% 2015:2018) %>%
  # aggregate over calendar time
  group_by(ethnicity, time, year) %>%
  summarise(n_tb_ETS = sum(n_tb_ETS, na.rm = TRUE)) |>
  rename(year_issue = year)


combined_dat <-
  full_join(dat_with_tb, ETSm,
             by = c("time", "ethnicity", "year_issue")) |>
  mutate(year_issue = as.factor(year_issue))

save(combined_dat, file = "data output/combined_dat.RData")

## aggregate year of issue

combined_a <-
  combined_dat |>
  group_by(ethnicity, time) |>
  summarise(n_tb = sum(n_tb, na.rm = TRUE),
            n_tb_ETS = sum(n_tb_ETS, na.rm = TRUE))


########
# plots
########

library(ggplot2)

ggplot(combined_dat, aes(time, n_tb, col = year_issue)) +
  geom_line() +
  xlim(0,8) +
  theme_bw() +
  facet_wrap(vars(ethnicity), scales = "free_y")
  geom_ribbon(aes(x = time, ymin = n_tb_low, ymax = n_tb_upp, group = ethnicity),
              inherit.aes = FALSE,
              linetype = 0,
              alpha = 0.1) +


ggdat <-
  melt(combined_agg, measure.vars = c("n_tb", "n_tb_ETS")) %>%
  mutate(value = ifelse(is.na(value), 0, value))

ggplot(ggdat, aes(time, value, col = variable, linetype = ethnicity)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(data = combined_agg, aes(x = time, ymin = n_tb_low, ymax = n_tb_upp, group = ethnicity),
              inherit.aes = FALSE,
              linetype = 0,
              alpha = 0.1) +
facet_wrap(vars(ethnicity), scales = "free_y")

