
# Evidence synthesis LTBI screening: jags
# with real data


library(R2jags)
library(R2WinBUGS)
library(purrr)
library(readr)


edetecttb <- read_csv("../../data/clean_edetecttb_11oct_withregionMAR.csv")
predicttb <-
  read.csv("E:/DIDE-PC_2019/NIHR_HTA_LTBI_project/data/UK-PREDICT-TB/041120-Rishi/PREDICT_data_for_NG.csv")



jags_dat_input <-
  list(
    len_gp = nrow(dat),                   # number of groups
    Xm  = dat$n_migrant,                  # number of migrants
    Xp  = dat$n_pos,                      # number of positive test results
    # Xn  = dat$n_migrant - dat$n_pos,      # number of negative test results
    Xtb = dat$n_tb                        # number of observed active tb cases
  )

params <-
  c("sens", "spec",
    "lambda",
    "p_latent"
    # "delta",
    # "pred_p", "pred_tb", "pred_latent",
    # "OR"
    # "ppost",
    # "delta_c",
    # "thresh"
  )

# inits <- function(){
#   list(
#     sens = ,
#     spec = ,
#     delta =
#   )
# }


#test
# n_iter <- 10000
# n_burnin <- 100
# n_thin <- 10

n_iter <- 1e6
n_burnin <- 1e3
n_thin <- 1e2 #floor((n_iter - n_burnin)/500)


##############
## run MCMC ##
##############

out <- jags(jags_dat_input,
            parameters.to.save = params,
            model.file = here::here("scripts", "BUGS_code.txt"),
            n.chains = 2,
            n.iter = n_iter,
            n.burnin = n_burnin,
            n.thin = n_thin,
            DIC = TRUE,
            working.directory = here::here("scripts"),
            progress.bar = "text")

BUGSoutput <- out$BUGSoutput

folder_nm <- "BUGSoutput"
dir.create(here::here("data output", folder_nm), showWarnings = FALSE)

save(BUGSoutput, file = here::here(folder_nm, "BUGSoutput.RData"))
save(jags_dat_input, file = here::here(folder_nm, "jags_dat_input.RData"))
