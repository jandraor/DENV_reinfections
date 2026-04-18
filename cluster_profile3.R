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

source("./R_scripts/DGP.R")
source("./R_scripts/inference.R")

data_obj <- get_inference_data()

age_inf_data_list <- data_obj$age_inf_data_list
titre_data_list   <- data_obj$titre_data_list
final_age_vctr    <- data_obj$final_age_vctr
#-------------------------------------------------------------------------------
n_indiv <- 1e4

ll_df   <- get_loglik_2("alternative")

box <- ll_df |> filter(ll > max(ll)- 20) |> sapply(range)

plan(multisession, workers = availableCores())
#-------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript profile_ci.R <param>")
}

param <- args[1]

alternative_profile(param,
                    titre_data_list,
                    age_inf_data_list,
                    final_age_vctr,
                    n_indiv,
                    box,
                    computation = "parallel")
