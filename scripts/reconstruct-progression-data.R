
# Reconstructing data from Kaplan-Meier curve
# of TB prgression in total cohort
# from Abubakar (2018) UK-PREDICT-TB
# using:
# https://www.r-bloggers.com/2019/12/reconstructing-data-from-kaplan-meier-curves/


# Determine which rows the upper and lower values of each interval are
# need this for risk set calculation
#
find_interval_limits <- function(start_time,
                                 surv_time){

  if (max(surv_time) > max(start_time))
    stop("Intervals must span all survival times. Later interval missing.")

  interval <- 1
  end_time <- start_time[2]
  upper <- NULL

  for (i in seq_along(surv_time)){
    if (surv_time[i] >= end_time) {
      upper[interval] <- i - 1
      interval <- interval + 1
      end_time <- start_time[interval + 1]
    }
  }

  cbind(lower = c(1, upper + 1),
        upper = c(upper, length(surv_time)))
}


library(survHE)

digdata <- read.csv("raw-data/abubaker-2018_curve Dataset.csv")

# survival data
digdata$Dataset.y <- (100 - digdata$Dataset.y)/100
plot(digdata, ylim = c(0.95, 1), type = "l")

# sort in ascending time
digdata <- digdata[order(digdata$Dataset.x), ]

surv_inp <- data.frame(ID = 1:nrow(digdata),
                       time = digdata$Dataset.x,
                       survival = round(digdata$Dataset.y, 3))

# ensure monotonic
surv_inp$diff <- c(TRUE, diff(surv_inp$survival)<=0)

surv_inp <- surv_inp |>
  mutate(survival = ifelse(!diff, lag(survival), survival)) |>
  select(-diff)

# # check
# surv_inp$diff_mono <- c(TRUE, diff(surv_inp$mono)<=0)

surv_inp$time[1] <- 0

write.table(surv_inp, file = "raw-data/surv_inp.txt", row.names = FALSE)

# risk sets

interval_limits <-
  find_interval_limits(start_time = c(0, 1, 2, 3, 4, 5, 6),
                       surv_time = surv_inp$time)

# add extraction date specific values to numbers at risk data from table in the paper

atrisk <-
  data.frame(Interval = 1:6,
             time = c(0, 1, 2, 3, 4, 5),
             interval_limits,
             nrisk = c(6380, 6220, 5654, 3963, 1438, 36)) |>
  as.matrix()

atrisk

write.table(atrisk, file = "raw-data/nrisk_inp.txt", row.names = FALSE)

# digitize data

survHE::digitise(surv_inp = "raw-data/surv_inp.txt",
                 nrisk_inp = "raw-data/nrisk_inp.txt",
                 km_output = "raw-data/KMdata.txt",
                 ipd_output = "raw-data/KMdata.txt")

IPDdata <- read.table("raw-data/IPDdata.txt", header = TRUE)
KMdata <- read.table("raw-data/KMdata.txt", header = TRUE)

IPD <- as.data.frame(IPDdata)
KM.est <- survfit(Surv(time, event) ~ 1, data = IPD)
                  # type = "kaplan-meier")
plot(KM.est, ylim = c(0.95, 1))
lines(digdata, ylim = c(0.95, 1), type = "l", col = "red")
