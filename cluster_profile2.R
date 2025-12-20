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

MLE_vals    <- get_MLE_2()
MLE_vals$ll <- NULL

plan(multisession, workers = availableCores())
# lambda 1----------------------------------------------------------------------

lambda_1_vals <- seq(0.04, 0.065, by = 0.0001)

prof_objs <- construct_profile(fixed_vals        = lambda_1_vals,
                               fixed_pos         = 1,
                               MLE_vals          = MLE_vals,
                               titre_data_list   = titre_data_list,
                               age_inf_data_list = age_inf_data_list,
                               final_age_vctr    = final_age_vctr,
                               n_indiv           = 1e5)

# lambda 2----------------------------------------------------------------------

lambda_2_vals <- seq(0.07, 0.10, by = 0.0001)

prof2_objs <- construct_profile(fixed_vals        = lambda_2_vals,
                                fixed_pos         = 2,
                                MLE_vals          = MLE_vals,
                                titre_data_list   = titre_data_list,
                                age_inf_data_list = age_inf_data_list,
                                final_age_vctr    = final_age_vctr,
                                n_indiv           = 1e5)

# rho --------------------------------------------------------------------------

rho_vals <- seq(0.005, 0.011, by = 0.000025)

prof3_objs <- construct_profile(fixed_vals        = rho_vals,
                                fixed_pos         = 3,
                                MLE_vals          = MLE_vals,
                                titre_data_list   = titre_data_list,
                                age_inf_data_list = age_inf_data_list,
                                final_age_vctr    = final_age_vctr,
                                n_indiv           = 1e5)

# log_A0------------------------------------------------------------------------

log_A0_vals <- seq(1, 2, by = 0.005)

prof4_objs <- construct_profile(fixed_vals        = log_A0_vals,
                                fixed_pos         = 4,
                                MLE_vals          = MLE_vals,
                                titre_data_list   = titre_data_list,
                                age_inf_data_list = age_inf_data_list,
                                final_age_vctr    = final_age_vctr,
                                n_indiv           = 1e5)

# phi---------------------------------------------------------------------------

phi_vals <- seq(5.3, 5.8, by = 0.002)

prof5_objs <- construct_profile(fixed_vals        = phi_vals,
                                fixed_pos         = 5,
                                MLE_vals          = MLE_vals,
                                titre_data_list   = titre_data_list,
                                age_inf_data_list = age_inf_data_list,
                                final_age_vctr    = final_age_vctr,
                                n_indiv           = 1e5)

# sd_1--------------------------------------------------------------------------

sd_1_vals <- seq(2, 5, by = 0.01)

prof6_objs <- construct_profile(fixed_vals        = sd_1_vals,
                                fixed_pos         = 6,
                                MLE_vals          = MLE_vals,
                                titre_data_list   = titre_data_list,
                                age_inf_data_list = age_inf_data_list,
                                final_age_vctr    = final_age_vctr,
                                n_indiv           = 1e5)

# sd_2--------------------------------------------------------------------------

sd_2_vals <- seq(1.5, 4.5, by = 0.01)

prof7_objs <- construct_profile(fixed_vals        = sd_2_vals,
                                fixed_pos         = 7,
                                MLE_vals          = MLE_vals,
                                titre_data_list   = titre_data_list,
                                age_inf_data_list = age_inf_data_list,
                                final_age_vctr    = final_age_vctr,
                                n_indiv           = 1e5)
