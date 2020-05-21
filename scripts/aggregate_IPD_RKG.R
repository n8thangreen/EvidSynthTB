#
# aggregate_IPT.R
# N Green
#
# Group together individuals in the PREDICT-TB dataset
# to use in LTBI screening DES model

setwd("~/Rishi Gupta/Peter White NIHR HTA")

library(dplyr)

## use this to test the code ----
# 
# dummy data
# DATA_MAIN <-
#   data.frame(
#     age = c(10,20,30,40,10,20,30,40,60,70),
#     sex = c(0,1,0,1,0,1,0,1,0,1),
#     ethnicity = c(1,2,3,4,5,1,2,3,4,5),
#     countryofbirth = c("UK", "India", "Other","UK", "India",
#                        "Other", "UK", "India", "Other","UK"),
#     countryofbirthTBincidence = c(1, 2, 1, 3, 3, 4, 5, 4, 4, 5),
#     yearssinceentry = c(NA,1,2,NA,10,4,NA,2,10,NA),
#     prevbcg = c(1,1,0,1,1,0,1,1,0,1),
#     activetb_site = c(1,1,1,1,NA,0,0,0,0,NA),
#     reasonforscreening = c(1,1,2,2,3,1,1,2,2,3),
#     activetb = c(1,1,1,1,0,1,1,1,1,0),
#     qfn_plus_result = c(1,1,1,2,2,2,0,1,1,1)
#   )

#Read in data

DATA_MAIN<- read.csv("Abubakar2018.csv")
DATA_MAIN <- DATA_MAIN %>% mutate_all(na_if,"")

#Change composite IGRA variable (pos_either) to numeric

DATA_MAIN$pos_either <- as.character(DATA_MAIN$pos_either)
DATA_MAIN$pos_either[DATA_MAIN$pos_either=="Positive"] <- 1
DATA_MAIN$pos_either[DATA_MAIN$pos_either=="Negative"] <- 0
DATA_MAIN$pos_either <- as.numeric(DATA_MAIN$pos_either)

#Drop patients with missing IGRA (otherwise we were getting NAs in the LTBI column in the final res dataframe)

DATA_MAIN <- DATA_MAIN[complete.cases(DATA_MAIN[ , "pos_either",]),]

# create new fields with discrete groups
DATA_MAIN$age_grp <- cut(DATA_MAIN$age, c(0,15,35,55,100))
DATA_MAIN$yearssinceentry_grp <- cut(DATA_MAIN$yearssinceentry, c(0,5,10,100))

## run this on the real data ----

# variables to aggregate over
fields_groups <- c(
  "age_grp",
  "sex",
  "ethnicity",
  # "countryofbirth",
  "inc_cob_participant2",
  "yearssinceentry_grp",
  "prevbcg",
  #"activetb_site", This is the site of disease in the study participant (which is an outcome rather than exposure - I have therefore removed it)
  "reasonforscreening")

# variables to obtain totals with categories
fields_out <- c(
  "activetb",
  "pos_either")


## replace DATA_MAIN with whatever the data are actually called

dat_out <- DATA_MAIN[, fields_out]
dat_groups <- DATA_MAIN[, fields_groups]


# don't know how missing data is coded
# have assumed its an NA
#
# include NAs as a seperate category (I have added some more variables here where there is a significant amount of missing data)
dat_groups$yearssinceentry_grp <- factor(dat_groups$yearssinceentry_grp, exclude = "")
#dat_groups$activetb_site  <- factor(dat_groups$activetb_site, exclude = "")
dat_groups$prevbcg  <- factor(dat_groups$prevbcg, exclude = "")
dat_groups$inc_cob_participant2  <- factor(dat_groups$inc_cob_participant2, exclude = "")

# total category size
pop <- aggregate(rep(1, nrow(dat_out)),
                 by = as.list(dat_groups),
                 FUN = "sum") %>% 
  rename(pop = "x")

# active tb counts
n_tb <- aggregate(dat_out$activetb,
                  by = as.list(dat_groups),
                  FUN = function(x) sum(x == 1)) %>% 
  rename(tb = "x")

# ltbi positive counts
n_ltbi <- aggregate(dat_out$pos_either,
                    by = as.list(dat_groups),
                    FUN = function(x) sum(x == 1)) %>% 
  rename(ltbi = "x")

# combined array

res <- 
  merge(pop, n_tb,
        by = fields_groups) %>% 
  merge(n_ltbi,
        by = fields_groups)

save(res, file = "aggregated_data.RData")
write.csv(res, "aggregated_data.csv")
