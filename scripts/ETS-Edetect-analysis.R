
# combine tb data sets
#
# want to take progression rate from predict-tb
# and denominators from e-detect
# to then estimate number of expected cases by
# groups for each year after entry
# then we can compare with ETS
#
# ETS/e-detect: migrant numbers by age, ethnicity
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
##      matched for ETS and edetecttb

library(readr)
library(dplyr)
library(stringr)

ETS <- read_csv("../../data/ETS_for_PREDICT_comparison_new.csv")
edetecttb <- read_csv("../../data/clean_edetecttb_11oct_withregionMAR.csv")

# predict-tb model fit stats
load("~/R/EvidSynthTB/data output/brms_pltbi.RData")
load("~/R/EvidSynthTB/data output/p_ltbi_brms.RData")

# aggregate into 5 yr age groups
pltbi_agegrp <-
  p_ltbi_brms %>%
  mutate(age_grp = cut(age, breaks = seq(0, 95, 5), right = FALSE)) %>%
  group_by(age_grp, ethnicity) %>%
  summarise(mean_pltbi = mean(prob_ltbi.Estimate),
            low_pltbi = mean(prob_ltbi.Q2.5),
            upp_pltbi = mean(prob_ltbi.Q97.5))  ##TODO: how to best aggregate this statistic?

## clean EDETECT data
# total numbers screened per age group, ethnicity
edetecttb_a <-
  edetecttb %>%
  mutate(region = replace(region, region == "Sub_Saharan_Africa", "Black African"),
         region = replace(region, region == "Bangladesh", "Bangladeshi"),
         region = replace(region, region == "India", "Indian"),
         region = replace(region, region == "Pakistan", "Pakistani")) %>%
  ## filter each ethnicity by different calendar years
  filter(year %in% 2015:2018) %>%
  ## different aggregations
  # group_by(sex, age_grp, year, visa_type, region) %>%  # non-aggregated
  # group_by(age_grp, year, region) %>%                  # with calendar time
  group_by(age_grp, region) %>%                         # without calendar time
  summarise(n_screened = sum(n_screened)) %>%
  rename(ethnicity = region)

# combine data and estimate number ltbi
edetecttb_a <-
  edetecttb_a %>%
  left_join(pltbi_agegrp, by = c("age_grp", "ethnicity")) %>%
  mutate(n_ltbi = n_screened*mean_pltbi,
         n_ltbi_low = n_screened*low_pltbi,
         n_ltbi_upp = n_screened*upp_pltbi) %>%
  full_join(data.frame(time = 1:10), by = character()) %>%
  na.omit()

##TODO: how to do this for posterior samples?
##      assume p tb is small so point estimate?

# fitted values
ppred <- stan_output$ppred

S_m <-
  summary(fs_m, t = 0:20) %>%
  ##TODO: is this the right stat to use?
  mutate(diff = lag(est) - est)

edetecttb_with_est <-
  full_join(edetecttb_a, S_m, by = "time") %>%
  mutate(n_tb = round(n_ltbi*diff, 2))

save(edetecttb_with_est, file = "edetecttb_with_estimates.RData")
write.csv(edetecttb_with_est, file = "edetecttb_with_estimates.csv")


###################
# compare with ETS

library(readr)
library(reshape2)

ETS <- read_csv("../data/ETS_for_PREDICT_comparison_new.csv")

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
  filter(year %in% 2015:2018) %>%
  # aggregate over calendar time
  group_by(age_grp, ethnicity, time) %>%
  summarise(n_tb_ETS = sum(n_tb_ETS, na.rm = TRUE))


combined_dat <-
  inner_join(edetecttb_with_est, ETSm,
             by = c("time", "age_grp", "ethnicity"))     # aggregated over calendar year
             # by = c("time", "age_grp", "year", "ethnicity"))

save(combined_dat, file = "combined_dat.RData")


##########
# plots

library(ggplot2)

ggdat <-
  melt(combined_dat, measure.vars = c("n_tb", "n_tb_ETS")) %>%
  mutate(value = ifelse(is.na(value), 0, value)) #%>%
  ## single age easier to read
  # filter(age_grp == "[40,45)")

ggplot(ggdat, aes(time, value, col = variable, linetype = ethnicity)) +
  geom_line() +
  theme_bw() +
  # facet_grid(age_grp ~ year, scales = "free_y")
  facet_wrap(vars(age_grp), scales = "free_y")


