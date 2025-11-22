library(boot)
library(dplyr)
library(flavipack)
library(furrr)
library(future)
library(lubridate)
library(pomp)
library(purrr)
library(progressr)
library(readr)
library(stringr)
library(tidyr)

source("./R_scripts/CPC.R")
source("./R_scripts/inference.R")
source("./R_scripts/KFCS.R")

CPC_age_inf_df <- CPC_get_prob_inf_at_age() |> filter(age >= 2)
final_age_1    <- CPC_age_inf_df$age |> max()

KFCS_age_inf_df <- KFCS_get_prob_inf_at_age() |> filter(age_round >= 2)
final_age_2     <- KFCS_age_inf_df$age_round |> max()

loglik_df <- find_MLE(CPC_age_inf_df, KFCS_age_inf_df,
                      final_age_1, final_age_2)
