library(boot)
library(future)
library(dplyr)
library(flavipack)
library(furrr)
library(lubridate)
library(progressr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

source("./R_scripts/CPC.R")
source("./R_scripts/inference.R")
source("./R_scripts/KFCS.R")


rho_vals <- seq(0.002, 0.012, by = 0.000025)

MLE_df <- readRDS("./saved_objects/inference/MLE/MLE_estimate.rds")

lambda_1_MLE <- MLE_df |> pull(lambda_1_MLE)
lambda_2_MLE <- MLE_df |> pull(lambda_2_MLE)
rho_MLE      <- MLE_df |> pull(rho_MLE)

MLE_vals <- c(lambda_1_MLE, lambda_2_MLE, rho_MLE)

CPC_age_inf_df <- CPC_get_prob_inf_at_age() |> filter(age >= 2)
final_age_1    <- CPC_age_inf_df$age |> max()

KFCS_age_inf_df <- KFCS_get_prob_inf_at_age() |> filter(age_round >= 2)
final_age_2     <- KFCS_age_inf_df$age_round |> max()

profile_rho_df <- maximise_over_profile(rho_vals,
                                        par_name    = "rho",
                                        MLE_vals    = MLE_vals,
                                        fixed_pos   = 3,
                                        data_df1    = CPC_age_inf_df,
                                        data_df2    = KFCS_age_inf_df,
                                        final_age_1 = final_age_1,
                                        final_age_2 = final_age_2) |>
  select(rho_MLE, ll) |>
  rename(value = rho_MLE) |>
  mutate(par_name = "rho")
