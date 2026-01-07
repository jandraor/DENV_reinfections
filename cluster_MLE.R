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

start_list <- get_starting_points("one_dataset")

plan(multisession, workers = availableCores())

res_list <- find_MLE(start_list, age_inf_data_list,
                     final_age_vctr, n_indiv = 1e4)

plan(sequential)
