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

source("./R_scripts/inference.R")

data_obj <- get_inference_data()

age_inf_data_list <- data_obj$age_inf_data_list
titre_data_list   <- data_obj$titre_data_list
final_age_vctr    <- data_obj$final_age_vctr
#-------------------------------------------------------------------------------
n_indiv <- 1e4

ll_df   <- get_loglik_2()

box <- ll_df |> filter(ll > max(ll)- 20) |> sapply(range)

plan(multisession, workers = availableCores())

colnames <- c("lambda_1", "lambda_2", "rho", "log_A0", "phi", "sd")

#-------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript profile_ci.R <param>")
}

param <- args[1]


arg_list <- list(titre_data_list   = titre_data_list,
                 age_inf_data_list = age_inf_data_list,
                 final_age_vctr    = final_age_vctr,
                 n_indiv           = n_indiv,
                 ds                = "two_datasets")


# lambda 1----------------------------------------------------------------------

if(param == "lambda_1")
{
  lambda_1_vals <- seq(0.05, 0.06, by = 0.0002)

  set.seed(2174)

  profile_design(
    lambda_1  = lambda_1_vals,
    lower = box[1, setdiff(colnames, param)],
    upper = box[2, setdiff(colnames, param)],
    nprof = 20, type = "sobol"
  ) -> guesses_df

  guesses_df <- guesses_df[, colnames]

  # starting points(sp)
  sp_list <- transpose(guesses_df)

  arg_list$fixed_pos       <- 1
  arg_list$starting_points <- sp_list

  prof_objs <- do.call(construct_profile, arg_list)
}

# lambda 2----------------------------------------------------------------------
if(param == "lambda_2")
{
  lambda_2_vals <- seq(0.07, 0.10, by = 0.001)

  set.seed(0954)

  guesses_df <- profile_design(
    lambda_2  = lambda_2_vals,
    lower = box[1, setdiff(colnames, param)],
    upper = box[2, setdiff(colnames, param)],
    nprof = 20, type = "sobol")

  guesses_df <- guesses_df[, colnames]

  # starting points(sp)
  sp_list <- transpose(guesses_df)

  prof2_objs <- construct_profile(starting_points   = sp_list,
                                  fixed_pos         = 2,
                                  titre_data_list   = titre_data_list,
                                  age_inf_data_list = age_inf_data_list,
                                  final_age_vctr    = final_age_vctr,
                                  n_indiv           = n_indiv)
}

# rho --------------------------------------------------------------------------
if(param == "rho")
{
  rho_vals <- seq(0.005, 0.011, by = 0.0002)

  set.seed(1002)

  guesses_df <- profile_design(
    rho   = rho_vals,
    lower = box[1, setdiff(colnames, param)],
    upper = box[2, setdiff(colnames, param)],
    nprof = 20, type = "sobol")

  guesses_df <- guesses_df[, colnames]

  # starting points(sp)
  sp_list <- transpose(guesses_df)

  prof3_objs <- construct_profile(starting_points   = sp_list,
                                  fixed_pos         = 3,
                                  titre_data_list   = titre_data_list,
                                  age_inf_data_list = age_inf_data_list,
                                  final_age_vctr    = final_age_vctr,
                                  n_indiv           = n_indiv)
}

# log_A0------------------------------------------------------------------------
if(param == "log_A0")
{
  log_A0_vals <- seq(1.1, 1.8, by = 0.01)

  set.seed(0836)

  guesses_df <- profile_design(
    log_A0 = log_A0_vals,
    lower  = box[1, setdiff(colnames, param)],
    upper  = box[2, setdiff(colnames, param)],
    nprof  = 20, type = "sobol")

  guesses_df <- guesses_df[, colnames]

  # starting points(sp)
  sp_list <- transpose(guesses_df)

  prof4_objs <- construct_profile(starting_points   = sp_list,
                                  fixed_pos         = 4,
                                  titre_data_list   = titre_data_list,
                                  age_inf_data_list = age_inf_data_list,
                                  final_age_vctr    = final_age_vctr,
                                  n_indiv           = n_indiv)
}

# phi---------------------------------------------------------------------------
if(param == "phi")
{
  phi_vals <- seq(5.9, 6.5, by = 0.01)

  set.seed(2109)

  guesses_df <- profile_design(
    phi    = phi_vals,
    lower  = box[1, setdiff(colnames, param)],
    upper  = box[2, setdiff(colnames, param)],
    nprof  = 20, type = "sobol")

  guesses_df <- guesses_df[, colnames]

  # starting points(sp)
  sp_list <- transpose(guesses_df)

  arg_list$fixed_pos       <- 5
  arg_list$starting_points <- sp_list

  prof5_objs <- do.call(construct_profile, arg_list)
}

# sd----------------------------------------------------------------------------
if(param == "sd")
{
  sd_vals <- seq(2.5, 3.2, by = 0.01)

  set.seed(2118)

  guesses_df <- profile_design(
    sd     = sd_vals,
    lower  = box[1, setdiff(colnames, param)],
    upper  = box[2, setdiff(colnames, param)],
    nprof  = 20, type = "sobol")

  guesses_df <- guesses_df[, colnames]

  # starting points(sp)
  sp_list <- transpose(guesses_df)

  arg_list$fixed_pos       <- 6
  arg_list$starting_points <- sp_list

  prof6_objs <- do.call(construct_profile, arg_list)
}
