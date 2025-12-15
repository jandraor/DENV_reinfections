library(boot)
library(flavipack)
library(future)
library(furrr)
library(dplyr)
library(lubridate)
library(nloptr)
library(pomp)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

source("./R_scripts/CPC.R")
source("./R_scripts/KFCS.R")
source("./R_scripts/inference.R")


KFCS_age_inf_df <- KFCS_get_prob_inf_at_age() |> filter(age_round >= 2)
final_age_1     <- KFCS_age_inf_df$age_round |> max()

KFCS_titre_df <- KFCS_get_titre_data() |> filter(age >= 2)

KFCS_avg_titre <- KFCS_titre_df |> group_by(age) |>
  summarise(mean = mean(log2_mean),
            n    = n())

CPC_age_inf_df <- CPC_get_prob_inf_at_age() |> filter(age >= 2)
final_age_2    <- CPC_age_inf_df$age |> max()

CPC_titre  <- CPC_get_tidy_data() |> filter(age >= 2)

CPC_avg_titre <- CPC_titre |> group_by(age) |>
  summarise(mean = mean(log2_mean),
            n    = n())


age_inf_data_list <- list(KFCS_age_inf_df, CPC_age_inf_df)
final_age_vctr    <- c(final_age_1, final_age_2)
titre_data_list   <- list(KFCS_avg_titre, CPC_avg_titre)

set.seed(111020225)

start_df <- sobol_design(
  lower = c("lambda_1" = 0.01,
            "lambda_2" = 0.01,
            "rho"      = 0.001,
            "r1"       = 0.001,
            "log_A0"   = 0.1,
            "phi"      = 0,
            "sd_1"     = 0.01,
            "sd_2"     = 0.01),
  upper = c("lambda_1" = 0.25,
            "lambda_2" = 0.25,
            "rho"      = 0.25,
            "r1"       = 0.99,
            "log_A0"   = 2,
            "phi"      = 5,
            "sd_1"     = 10,
            "sd_2"     = 10),
  nseq =  200)

start_list <- transpose(start_df)

plan(multisession, workers = availableCores())

res_list <- find_MLE(start_list, titre_data_list,
                     age_inf_data_list, final_age_vctr)

plan(sequential)
