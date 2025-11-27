source("./R_scripts/titre_simulators.R")

library(dplyr)
library(flavipack)
library(furrr)
library(future)
library(purrr)
library(stringr)
library(tidyr)

max_n_inf <- 20
damp_df <- data.frame(seq_inf = 1:max_n_inf,
                      decay_rate = pmax(seq(0.2, by = -0.05,
                                            length.out = max_n_inf), 0.05),
                      dynamics = "Dampening")
sim_titres_list <- simulate_titres_scenarios(damp_df)
